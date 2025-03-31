#include "mfem.hpp"
#include "heston-slv-dpg.hpp" // Include your solver class header
// #include "dpg-error-estimator.hpp" // Needed later for test function
#include <iostream>
#include <fstream>

using namespace mfem;
using namespace std;

// Forward declaration for the test function we will create
// bool TestDPGErrorEstimator(int order, int dim, int ref_levels);
bool TestPoissonDPG(int order, int ref_levels, bool visualization);


int main(int argc, char *argv[])
{
    // Initialize MPI.
    Mpi::Init(argc, argv);
    int myid = Mpi::WorldRank();
    int num_procs = Mpi::WorldSize();

    // Default options
    int order = 1;
    int ref_levels = 1;
    bool visualization = true;
    bool run_poisson_test = true; // Flag to run the Poisson test

    OptionsParser args(argc, argv);
    args.AddOption(&order, "-o", "--order", "Finite element order (polynomial degree).");
    args.AddOption(&ref_levels, "-r", "--ref-levels", "Number of mesh refinement levels.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis", "--no-visualization", "Enable or disable GLVis visualization.");
    args.AddOption(&run_poisson_test, "-test-poisson", "--test-poisson", "-no-test-poisson", "--no-test-poisson", "Run the Poisson DPG test case.");

    args.Parse();
    if (!args.Good())
    {
        if (myid == 0)
        {
            args.PrintUsage(cout);
        }
        Mpi::Finalize();
        return 1;
    }
    if (myid == 0)
    {
        args.PrintOptions(cout);
    }

    try
    {
        if (run_poisson_test)
        {
             if (myid == 0) {
                 cout << "\nRunning Poisson DPG Test Case..." << endl;
             }
             bool poisson_passed = TestPoissonDPG(order, ref_levels, visualization);
             if (myid == 0) {
                 cout << "Poisson DPG Test " << (poisson_passed ? "PASSED" : "FAILED") << endl;
             }
        }
        else
        {
            if (myid == 0) {
                cout << "\nSkipping Poisson DPG Test Case." << endl;
                cout << "Heston Solver execution is currently a stub." << endl;
            }
            // Placeholder for Heston solver execution
            // UltraweakDPGHestonSLV heston_solver(order, 1); // Example
            // heston_solver.Initialize();
            // heston_solver.Solve();
        }

        // Placeholder for Error Estimator test (Phase 1)
        // if (myid == 0) {
        //     cout << "\nRunning Error Estimator Test..." << endl;
        // }
        // bool est_passed = TestDPGErrorEstimator(order, 2, ref_levels); // 2D test
        // if (myid == 0) {
        //      cout << "Error Estimator Test " << (est_passed ? "PASSED" : "FAILED") << endl;
        // }

    }
    catch (exception &e)
    {
        if (myid == 0)
        {
            cout << "Standard Exception Caught: " << e.what() << endl;
        }
        Mpi::Finalize();
        return 1;
    }

    Mpi::Finalize();
    return 0;
}

// Define TestPoissonDPG in test_poisson_dpg.cpp