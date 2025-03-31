#include "heston-slv-multigrid.hpp"

#define MFEM_MULTIGRID_HAS_PRINT_LEVEL

namespace mfem {

HestonSLVMultigrid::HestonSLVMultigrid(FiniteElementSpaceHierarchy &fespaces,
                                       Array<int> &ess_bdr, double r_,
                                       double q_, double kappa_, double theta_,
                                       double xi_, double rho_)
    : GeometricMultigrid(fespaces, ess_bdr), r(r_), q(q_), kappa(kappa_),
      theta(theta_), xi(xi_), rho(rho_) {
  // Create operators and solvers for each level of the hierarchy
  ConstructCoarseOperatorAndSolver(fespaces.GetFESpaceAtLevel(0));

  for (int level = 1; level < fespaces.GetNumLevels(); ++level) {
    ConstructOperatorAndSmoother(fespaces.GetFESpaceAtLevel(level), level);
  }
}

HestonSLVMultigrid::~HestonSLVMultigrid() {
  // Free all bilinear forms
  for (int i = 0; i < bfs.Size(); i++) {
    delete bfs[i];
  }
}

void HestonSLVMultigrid::ConstructBilinearForm(FiniteElementSpace &fespace) {
  // Create coefficient functions for the Heston-SLV PDE
  ConstantCoefficient r_coef(r);
  ConstantCoefficient neg_r_coef(-r);
  ConstantCoefficient one(1.0);

  // S-diffusion: 0.5*V*S^2
  VectorFunctionCoefficient s_drift_vec(2, [this](const Vector &x, Vector &v) {
    // Drift in S-direction: (r-q)*S
    v(0) = (r - q) * x(0); // S-component
    v(1) = 0.0;            // V-component
  });

  VectorFunctionCoefficient v_drift_vec(2, [this](const Vector &x, Vector &v) {
    // Drift in V-direction: kappa*(theta-V)
    v(0) = 0.0;                    // S-component
    v(1) = kappa * (theta - x(1)); // V-component
  });

  // Create the coefficients for diffusion terms
  FunctionCoefficient s_diff([this](const Vector &x) {
    return 0.5 * x(1) * x(0) * x(0); // 0.5*V*S^2
  });

  FunctionCoefficient v_diff([this](const Vector &x) {
    return 0.5 * xi * xi * x(1); // 0.5*xi^2*V
  });

  FunctionCoefficient mixed_diff([this](const Vector &x) {
    return rho * xi * x(1) * x(0); // rho*xi*V*S
  });

  // Create and set up the bilinear form
  BilinearForm *form = new BilinearForm(&fespace);

  // Set partial assembly for better performance when available
  bool use_partial = true;
  try {
    form->SetAssemblyLevel(AssemblyLevel::PARTIAL);
  } catch (...) {
    use_partial = false;
  }

  // Add integrators for all terms in the Heston-SLV PDE

  // Reaction term: -r*u
  form->AddDomainIntegrator(new MassIntegrator(neg_r_coef));

  // Convection terms: (r-q)*S * du/dS and kappa*(theta-V) * du/dV
  form->AddDomainIntegrator(new ConvectionIntegrator(s_drift_vec));
  form->AddDomainIntegrator(new ConvectionIntegrator(v_drift_vec));

  // Diffusion terms
  form->AddDomainIntegrator(
      new DiffusionIntegrator(s_diff)); // 0.5*V*S^2 * d^2u/dS^2
  form->AddDomainIntegrator(
      new DiffusionIntegrator(v_diff)); // 0.5*xi^2*V * d^2u/dV^2

  // Mixed derivative term: rho*xi*V*S * d^2u/dSdV
  form->AddDomainIntegrator(new MixedCrossDerivativeIntegrator(mixed_diff));

  // Add stabilization/interface terms if using DG spaces
  const FiniteElementCollection *fec = fespace.FEColl();
  if (dynamic_cast<const DG_FECollection *>(fec)) {
    form->AddInteriorFaceIntegrator(
        new DGDiffusionIntegrator(s_diff, -1.0, 10.0));
    form->AddInteriorFaceIntegrator(
        new DGDiffusionIntegrator(v_diff, -1.0, 10.0));

    // For mixed derivatives, we would need specialized interface terms
    // This is a simplification for the multigrid operator
  }

  form->Assemble();

  if (!use_partial) {
    form->Finalize();
  }

  bfs.Append(form);
}

void HestonSLVMultigrid::ConstructCoarseOperatorAndSolver(
    FiniteElementSpace &coarse_fespace) {
  ConstructBilinearForm(coarse_fespace);

  OperatorPtr opr;
  opr.SetType(Operator::ANY_TYPE);
  bfs[0]->FormSystemMatrix(*essentialTrueDofs[0], opr);
  opr.SetOperatorOwner(false);

  // Use direct solver on the coarsest level
  Solver *direct_solver = nullptr;

// Choose appropriate direct solver based on availability
#ifdef MFEM_USE_SUITESPARSE
  direct_solver = new UMFPackSolver();
#elif defined(MFEM_USE_MUMPS)
  direct_solver = new MUMPSSolver();
#else
  // Fallback to dense solver for small systems, or GS iterations
  if (coarse_fespace.GetTrueVSize() < 1000) {
    direct_solver = new DenseMatrixInverseOperator();
  } else {
    GSSmoother *smoother = new GSSmoother;
    smoother->SetOperator(*opr.Ptr());
    smoother->SetPositiveDiagonal(true);
    direct_solver = smoother;
  }
#endif

  if (direct_solver) {
    direct_solver->SetOperator(*opr.Ptr());
  } else {
    direct_solver = new GSSmoother();
    direct_solver->SetOperator(*opr.Ptr());
  }

  AddLevel(opr.Ptr(), direct_solver, true, true);
}

void HestonSLVMultigrid::ConstructOperatorAndSmoother(
    FiniteElementSpace &fespace, int level) {
  const Array<int> &ess_tdof_list = *essentialTrueDofs[level];
  ConstructBilinearForm(fespace);

  OperatorPtr opr;
  opr.SetType(Operator::ANY_TYPE);
  bfs[level]->FormSystemMatrix(ess_tdof_list, opr);
  opr.SetOperatorOwner(false);

  // Create appropriate smoother based on the problem characteristics
  Solver *smoother = nullptr;

  // Extract diagonal for Chebyshev smoother or Jacobi relaxation
  Vector diag(fespace.GetTrueVSize());
  bfs[level]->AssembleDiagonal(diag);

  // Check if the operator is a HypreParMatrix
  bool is_parallel = dynamic_cast<const HypreParMatrix *>(opr.Ptr()) != nullptr;

  if (is_parallel) {
    // For parallel case, use Chebyshev smoother
    smoother =
        new OperatorChebyshevSmoother(*opr.Ptr(), diag, ess_tdof_list, 2);
  } else {
    // For serial case, use SOR or Chebyshev
    // Determine if PDE is convection-dominated or diffusion-dominated
    double convection_dominance =
        std::min(std::abs(r - q), kappa) / (0.5 * xi * xi * theta + 0.25);

    if (convection_dominance > 10.0) {
      // Convection-dominated: use Gauss-Seidel with appropriate ordering
      smoother = new GSSmoother();
    } else {
      // Diffusion-dominated: use Chebyshev
      smoother =
          new OperatorChebyshevSmoother(*opr.Ptr(), diag, ess_tdof_list, 2);
    }
  }

  // Set the operator for the smoother
  smoother->SetOperator(*opr.Ptr());

  // Add the level to the hierarchy
  AddLevel(opr.Ptr(), smoother, true, true);
}

// Implementation of MixedCrossDerivativeIntegrator

MixedCrossDerivativeIntegrator::MixedCrossDerivativeIntegrator(Coefficient &q)
    : Q(q) {}

void MixedCrossDerivativeIntegrator::AssembleElementMatrix(
    const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat) {
  int dof = el.GetDof();
  int dim = el.GetDim();

  MFEM_ASSERT(dim == 2, "MixedCrossDerivativeIntegrator requires 2D elements");

  DenseMatrix dshape(dof, dim);
  Vector shape(dof);

  elmat.SetSize(dof);
  elmat = 0.0;

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = 2 * el.GetOrder() + Trans.OrderW();
    ir = &IntRules.Get(el.GetGeomType(), order);
  }

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    Trans.SetIntPoint(&ip);

    el.CalcShape(ip, shape);
    el.CalcDShape(ip, dshape);

    double coef = Q.Eval(Trans, ip);
    double w = ip.weight * coef * Trans.Weight();

    // Construct the mixed derivative matrix: d^2/(dS*dV) term
    for (int j = 0; j < dof; j++) {
      for (int k = 0; k < dof; k++) {
        // Mixed derivatives: du/dS * dv/dV + du/dV * dv/dS
        elmat(j, k) +=
            w * (dshape(j, 0) * dshape(k, 1) + dshape(j, 1) * dshape(k, 0));
      }
    }
  }
}

// Implementation of DPGSystemMultigrid

DPGSystemMultigrid::DPGSystemMultigrid(BlockOperator *block_operator,
                                       Array<int> &block_offsets,
                                       Array<Mesh *> &meshes,
                                       Array<double> &heston_params,
                                       int max_levels)
    : Solver(block_operator->Height()), A(block_operator),
      offsets(block_offsets), num_vars(block_offsets.Size() - 1),
      r(block_offsets), z(block_offsets) {
  MFEM_ASSERT(heston_params.Size() >= 6,
              "Heston parameters array must have at least 6 elements");
  MFEM_ASSERT(meshes.Size() == num_vars,
              "Number of meshes must match number of variables");

  double r_value = heston_params[0];
  double q_value = heston_params[1];
  double kappa_value = heston_params[2];
  double theta_value = heston_params[3];
  double xi_value = heston_params[4];
  double rho_value = heston_params[5];

  // Initialize space hierarchies and multigrid solvers for each variable
  for (int var = 0; var < num_vars; var++) {
    // Create FE space hierarchy
    FiniteElementSpaceHierarchy *hierarchy = CreateHestonHPHierarchy(
        meshes[var],
        4, // max_order
        2, // h_levels - adjust based on memory constraints
        1  // p_levels - typically 1 for DPG
    );

    // Create empty array for essential BCs
    Array<int> ess_bdr;

    // Create multigrid solver for this variable
    mg_solvers.push_back(std::make_unique<HestonSLVMultigrid>(
        *hierarchy, ess_bdr, r_value, q_value, kappa_value, theta_value,
        xi_value, rho_value));

    // Set multigrid parameters through the base class methods
    MultigridBase *mg_base = mg_solvers[var].get();
    mg_base->SetCycleType(MultigridBase::CycleType::VCYCLE, 1, 1);
  }
}

DPGSystemMultigrid::~DPGSystemMultigrid() {
  // mg_solvers are managed by unique_ptr and cleaned up automatically
}

void DPGSystemMultigrid::Mult(const Vector &x, Vector &y) const {
  // Apply the multigrid preconditioner block by block
  r = 0.0;
  z = 0.0;

  // Copy x into the block vector r
  for (int var = 0; var < num_vars; var++) {
    // Create a temporary vector to hold the block
    Vector x_block(offsets[var + 1] - offsets[var]);
    for (int i = 0; i < x_block.Size(); i++) {
      x_block(i) = x(offsets[var] + i);
    }
    r.GetBlock(var) = x_block;
  }

  // Apply multigrid preconditioner to each block
  for (int var = 0; var < num_vars; var++) {
    mg_solvers[var]->Mult(r.GetBlock(var), z.GetBlock(var));
  }

  // Copy z into the output vector y
  for (int var = 0; var < num_vars; var++) {
    Vector &z_block = z.GetBlock(var);
    for (int i = 0; i < z_block.Size(); i++) {
      y(offsets[var] + i) = z_block(i);
    }
  }
}

void DPGSystemMultigrid::SetPrintLevel(int print_lvl) {
  // Your MFEM version doesn't have SetPrintLevel, so we just store the value
  // but don't try to apply it to the multigrid solvers
  // This method is kept as a placeholder for compatibility
}

// Implementation of CreateHestonHPHierarchy

FiniteElementSpaceHierarchy *
CreateHestonHPHierarchy(Mesh *mesh, int max_order, int h_levels, int p_levels) {
  // Create base collection for the coarsest level
  FiniteElementCollection *base_fec = new DG_FECollection(1, mesh->Dimension());

  // Create space for the coarsest level
  FiniteElementSpace *coarse_space = new FiniteElementSpace(mesh, base_fec);

  // Create the hierarchy
  FiniteElementSpaceHierarchy *hierarchy =
      new FiniteElementSpaceHierarchy(mesh, coarse_space, true, true);

  // Add h-refinement levels
  for (int level = 0; level < h_levels; ++level) {
    hierarchy->AddUniformlyRefinedLevel();
  }

  // Add p-refinement levels
  Array<FiniteElementCollection *> collections;
  collections.Append(base_fec);

  for (int level = 0; level < p_levels; ++level) {
    // Increase polynomial order for each p-level
    int order = std::min(2 << level, max_order);
    collections.Append(new DG_FECollection(order, mesh->Dimension()));
    hierarchy->AddOrderRefinedLevel(collections.Last());
  }

  return hierarchy;
}

} // namespace mfem