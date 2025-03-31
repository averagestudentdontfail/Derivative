#include "mfem.hpp"
#include "dpg-integrators.hpp" // Include our custom integrators
#include <iostream>
#include <fstream>
#include <memory> // For unique_ptr

using namespace mfem;
using namespace std;

// Function signature matches the call in heston-slv-main.cpp
bool TestPoissonDPG(int order, int ref_levels, bool visualization)
{
    int myid = Mpi::WorldRank();
    int num_procs = Mpi::WorldSize();

    if (myid == 0) {
        cout << "\n--- Testing Poisson DPG ---" << endl;
        cout << "  Order: " << order << ", Ref Levels: " << ref_levels << endl;
    }

    // 1. Mesh
    Mesh serial_mesh = Mesh::MakeCartesian2D(4, 4, Element::QUADRILATERAL, true); // Small base mesh
    for (int i = 0; i < ref_levels; ++i) {
        serial_mesh.UniformRefinement();
    }
    ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
    serial_mesh.Clear(); // Free memory
    int dim = mesh.Dimension();

    if (myid == 0) {
        cout << "  Mesh NE: " << mesh.GetNE() << endl;
    }

    // 2. Define FE Collections and Spaces
    enum Var { U = 0, SIGMA = 1, NUM_VARS = 2 }; // Trial variables
    enum TestVar { V = 0, W = 1, NUM_TEST_VARS = 2 }; // Test variables

    int test_order = order + 1; // Enriched test space

    // Trial spaces
    FiniteElementCollection *u_fec = new DG_FECollection(order, dim);
    FiniteElementCollection *sigma_fec = new DG_FECollection(order, dim, BasisType::GaussLobatto, FiniteElement::VECTOR);
    ParFiniteElementSpace u_fes(&mesh, u_fec);
    ParFiniteElementSpace sigma_fes(&mesh, sigma_fec);

    Array<ParFiniteElementSpace *> trial_fes(NUM_VARS);
    trial_fes[U] = &u_fes;
    trial_fes[SIGMA] = &sigma_fes;

    Array<int> trial_block_offsets(NUM_VARS + 1);
    trial_block_offsets[0] = 0;
    trial_block_offsets[1] = u_fes.GetTrueVSize();
    trial_block_offsets[2] = sigma_fes.GetTrueVSize();
    trial_block_offsets.PartialSum();

    // Test spaces
    FiniteElementCollection *v_fec = new DG_FECollection(test_order, dim); // Need H(div) later? Start with DG.
    FiniteElementCollection *w_fec = new DG_FECollection(test_order, dim, BasisType::GaussLobatto, FiniteElement::VECTOR); // Need H(curl) later? Start with DG.
    ParFiniteElementSpace v_fes(&mesh, v_fec);
    ParFiniteElementSpace w_fes(&mesh, w_fec);

    Array<ParFiniteElementSpace *> test_fes(NUM_TEST_VARS);
    test_fes[V] = &v_fes;
    test_fes[W] = &w_fes;

    Array<int> test_block_offsets(NUM_TEST_VARS + 1);
    test_block_offsets[0] = 0;
    test_block_offsets[1] = v_fes.GetTrueVSize();
    test_block_offsets[2] = w_fes.GetTrueVSize();
    test_block_offsets.PartialSum();


    if (myid == 0) {
        cout << "  Trial DOFs: u=" << trial_block_offsets[1] << ", sigma=" << trial_block_offsets[2] - trial_block_offsets[1] << ", Total=" << trial_block_offsets[2] << endl;
        cout << "  Test DOFs:  v=" << test_block_offsets[1] << ", w=" << test_block_offsets[2] - test_block_offsets[1] << ", Total=" << test_block_offsets[2] << endl;
    }

    // 3. Define Bilinear/Linear Forms (DPG System Components)

    // --- RHS Linear Form F ---
    unique_ptr<BlockVector> F(new BlockVector(test_block_offsets));
    *F = 0.0;

    ConstantCoefficient one(1.0);
    unique_ptr<ParLinearForm> f_form(new ParLinearForm);
    f_form->Update(test_fes[V], F->GetBlock(V), 0); // Accumulate into V block
    f_form->AddDomainIntegrator(new DomainLFIntegrator(one)); // Source term f = 1 for -div(sigma) = f
    f_form->Assemble();


    // --- Test Inner Product G ---
    unique_ptr<BlockOperator> G(new BlockOperator(test_block_offsets));
    G->owns_blocks = true; // G will delete the blocks

    unique_ptr<ParBilinearForm> g00(new ParBilinearForm(test_fes[V]));
    g00->AddDomainIntegrator(new TestInnerProductIntegrator(1.0)); // (v, v')_L2
    g00->Assemble();
    g00->Finalize();
    G->SetBlock(V, V, g00.release()); // G takes ownership

    unique_ptr<ParBilinearForm> g11(new ParBilinearForm(test_fes[W]));
    g11->AddDomainIntegrator(new TestInnerProductIntegrator(1.0)); // (w, w')_L2
    g11->Assemble();
    g11->Finalize();
    G->SetBlock(W, W, g11.release()); // G takes ownership


    // --- Coupling Form B ---
    // Equation 1: sigma + grad u = 0 -> (sigma, w) - (u, div w) + <{u}, [w.n]> = 0 (sign flip for grad)
    // Equation 2: -div sigma = f -> (-div sigma, v) = (f, v)
    unique_ptr<BlockOperator> B(new BlockOperator(test_block_offsets, trial_block_offsets));
    B->owns_blocks = true; // B will delete the blocks

    // B(0,0) = (v, u): Zero term
    // B(0,1) = (v, sigma): (-div sigma, v) -> using MixedVectorWeakDivergenceIntegrator
    unique_ptr<ParMixedBilinearForm> b01(new ParMixedBilinearForm(trial_fes[SIGMA], test_fes[V]));
    b01->AddDomainIntegrator(new MixedVectorWeakDivergenceIntegrator(one)); // Default sign is (div u, v), we want (-div sigma, v)
    // Need to handle boundary terms? -(sigma.n, v)_bdr. Assume natural BC for now (sigma.n = 0)
    b01->Assemble();
    b01->Finalize();
    B->SetBlock(V, SIGMA, b01.release());


    // B(1,0) = (w, u): -(u, div w) + <{u}, [w.n]>_int + <u, w.n>_bdr
    unique_ptr<ParMixedBilinearForm> b10(new ParMixedBilinearForm(trial_fes[U], test_fes[W]));
    ConstantCoefficient neg_one(-1.0);
    b10->AddDomainIntegrator(new MixedScalarWeakGradientIntegrator()); // -(u, div w) term handled by integrator sign? Check impl. YES, need neg_one coeff.
    // Correction: The integrator computes (u, div w). We need -(u, div w).
    // Let's redefine MixedScalarWeakGradientIntegrator to take a coefficient.
    // Or add the neg_one coefficient here. Let's try coefficient.
    // Reverting: No, the standard form is (grad u, w) -> -(u, div w) + bdr. Our form is (u, div w). So sign is correct.
    b10->AddInteriorFaceIntegrator(new UltraweakJumpIntegrator()); // Handles <{u}, [w.n]> and <u, w.n>_bdr
    // Boundary condition u=0 on bdr: The jump term <u, w.n>_bdr becomes <0, w.n> = 0.
    b10->Assemble();
    b10->Finalize();
    B->SetBlock(W, U, b10.release());

    // B(1,1) = (w, sigma): (sigma, w)
    unique_ptr<ParMixedBilinearForm> b11(new ParMixedBilinearForm(trial_fes[SIGMA], test_fes[W]));
    b11->AddDomainIntegrator(new VectorFEMassIntegrator(one)); // (sigma, w)
    b11->Assemble();
    b11->Finalize();
    B->SetBlock(W, SIGMA, b11.release());


    // 4. Set up Solvers
    HypreBoomerAMG G_prec(*G); // Simple AMG for the whole block G for now
    G_prec.SetPrintLevel(0);
    HyprePCG G_solver(G->GetComm());
    G_solver.SetTol(1e-12);
    G_solver.SetMaxIter(200);
    G_solver.SetPrintLevel(0);
    G_solver.SetOperator(*G);
    G_solver.SetPreconditioner(G_prec);

    SolverOperator G_inv_op(G_solver); // Wrap solver as Operator
    TransposeOperator B_T_op(B.get());

    // A = B^T G^{-1} B
    ProductOperator GinvB_op(&G_inv_op, B.get(), false, false);
    unique_ptr<ProductOperator> A(new ProductOperator(&B_T_op, &GinvB_op, false, false));

    // RHS_eff = B^T G^{-1} F
    BlockVector G_inv_F(test_block_offsets);
    G_inv_F = 0.0;
    G_solver.Mult(*F, G_inv_F); // Apply G^{-1} to F

    BlockVector RHS_eff(trial_block_offsets);
    RHS_eff = 0.0;
    B_T_op.Mult(G_inv_F, RHS_eff); // Apply B^T

    // Solver for A * U = RHS_eff
    HypreGMRES A_solver(A->GetComm());
    A_solver.SetTol(1e-8);
    A_solver.SetMaxIter(500);
    A_solver.SetPrintLevel(1); // Print GMRES iterations
    A_solver.SetOperator(*A);
    // No preconditioner for A yet

    // 5. Solve the system
    BlockVector U_sol(trial_block_offsets);
    U_sol = 0.0;
    A_solver.Mult(RHS_eff, U_sol);

    bool converged = (A_solver.GetFinalRelResidual() < 1e-8);
    if (myid == 0) {
        cout << "  GMRES Converged: " << (converged ? "YES" : "NO") << endl;
        cout << "  GMRES Iterations: " << A_solver.GetNumIterations() << endl;
        cout << "  Final Rel Residual: " << A_solver.GetFinalRelResidual() << endl;
    }

    // 6. Visualize Solution
    if (visualization)
    {
        char vishost[] = "localhost";
        int visport = 19916;
        socketstream vis_u, vis_s;

        ParGridFunction u_gf(&u_fes);
        u_gf.Distribute(U_sol.GetBlock(U));

        vis_u.open(vishost, visport);
        vis_u.precision(8);
        vis_u << "parallel " << num_procs << " " << myid << "\n";
        vis_u << "solution\n" << mesh << u_gf
              << "window_title 'Solution u (Poisson DPG)'" << flush;


        ParGridFunction sigma_gf(&sigma_fes);
        sigma_gf.Distribute(U_sol.GetBlock(SIGMA));

        vis_s.open(vishost, visport); // Use different port or close previous
        vis_s.precision(8);
        vis_s << "parallel " << num_procs << " " << myid << "\n";
        vis_s << "solution\n" << mesh << sigma_gf
              << "window_title 'Solution sigma (Poisson DPG)'" << flush;
    }


    // 7. Cleanup (unique_ptr handles B, G, F blocks, A)
    delete w_fec;
    delete v_fec;
    delete sigma_fec;
    delete u_fec;

    return converged; // Return true if solver converged
}