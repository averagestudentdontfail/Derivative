// File: /src/dpg-error-estimator.hpp

#ifndef DPG_ERROR_ESTIMATOR_HPP
#define DPG_ERROR_ESTIMATOR_HPP

#include "mfem.hpp"
#include <memory>
#include <stdexcept>
#include <vector>

namespace mfem {

/**
 * @brief Computes residual-based error estimates for ultraweak DPG formulations.
 *
 * This class calculates local error indicators eta_K = || F - B*U ||_{V_K'}
 * by solving local Riesz representation problems on each element K.
 * V_K' is the dual of the local test space V_K, and the norm is induced
 * by the test space inner product G.
 */
class UltraweakDPGErrorEstimator {
private:
    // Pointers to the global FE spaces (not owned)
    ParFiniteElementSpace *trial_fes_; // Block structure assumed later
    ParFiniteElementSpace *test_fes_;  // Block structure assumed later

    // Pointers to the global DPG operators (not owned)
    // These represent the full system B(u,v) = F(v)
    BlockOperator *B_op_; // Maps trial space U to test space dual V'
    Operator *G_op_;      // Represents test space inner product G(v,v)
    BlockVector *F_vec_;  // Represents the RHS linear form F(v)

    // Pointer to the global solution vector (not owned)
    const BlockVector *U_sol_;

    // Cached local error estimates
    Vector local_errors_;
    bool errors_computed_;

    // Parallel communication
    MPI_Comm comm_;
    int myid_;
    int nranks_;

    // Verbosity
    int verbose_ = 0;

    // Helper method to assemble local matrices/vectors needed for estimation
    // This is a key part that needs careful implementation.
    void GetElementData(int elem_idx, DenseMatrix &B_elem, DenseMatrix &G_elem,
                        Vector &F_elem, Vector &U_elem);

public:
    /**
     * @brief Constructor.
     * @param trial_fes Pointer to the parallel trial finite element space (block structure expected).
     * @param test_fes Pointer to the parallel test finite element space (block structure expected).
     * @param B Pointer to the global DPG operator B.
     * @param G Pointer to the global test space inner product operator G.
     * @param F Pointer to the global RHS vector F.
     * @param U Pointer to the global solution vector U.
     */
    UltraweakDPGErrorEstimator(ParFiniteElementSpace *trial_fes,
                               ParFiniteElementSpace *test_fes,
                               BlockOperator *B, Operator *G, BlockVector *F,
                               const BlockVector *U);

    virtual ~UltraweakDPGErrorEstimator() = default; // Use default destructor

    /**
     * @brief Compute the local error estimate for a single element.
     * @param elem_idx The index of the element.
     * @return The computed error estimate eta_K = sqrt(|r_elem * err_repr|).
     * @throws std::out_of_range If elem_idx is invalid.
     */
    virtual double GetLocalError(int elem_idx);

    /**
     * @brief Get the vector of local error estimates for all elements owned by this processor.
     * Recomputes errors if they haven't been computed since the last Reset().
     * @return A const reference to the vector of local error estimates.
     */
    virtual const Vector &GetLocalErrors();

    /**
     * @brief Clears cached error estimates, forcing recomputation on the next call to GetLocalErrors().
     * Should be called after the solution U or the mesh changes.
     */
    virtual void Reset();

    /**
     * @brief Set the verbosity level.
     * @param level 0 (silent), 1 (basic info), 2 (detailed info).
     */
    void SetVerbosity(int level) { verbose_ = level; }

    // Provide access to underlying components if needed by refiners
    ParFiniteElementSpace *GetTrialFESpace() const { return trial_fes_; }
    ParFiniteElementSpace *GetTestFESpace() const { return test_fes_; }
    const BlockVector* GetSolutionVector() const { return U_sol_; }
};


/**
 * @brief Performs adaptive mesh refinement based on DPG error estimates.
 */
class DPGAdaptiveRefiner {
private:
    UltraweakDPGErrorEstimator &estimator_; // Reference to the estimator

    // Refinement strategy parameters
    double total_error_fraction_; // For bulk (Doerfler) marking [0,1]
    double threshold_factor_;     // For max marking (fraction of max error) [0,1]
    int strategy_;                // 0 = Max strategy, 1 = Bulk/Doerfler strategy

    // Non-conforming refinement control
    int nc_limit_ = 0; // Max level of hanging nodes (0 means conforming)

    // Verbosity
    int verbose_ = 0;

public:
    /**
     * @brief Constructor.
     * @param est Reference to the UltraweakDPGErrorEstimator.
     * @param err_fraction Fraction of total squared error for bulk marking (default 0.7).
     * @param threshold Threshold factor relative to max error for max strategy (default 0.7).
     * @param strategy_type 0 for max strategy, 1 for bulk/Doerfler (default 1).
     */
    DPGAdaptiveRefiner(UltraweakDPGErrorEstimator &est, double err_fraction = 0.7,
                       double threshold = 0.7, int strategy_type = 1);

    /**
     * @brief Determines which elements should be refined based on the chosen strategy.
     * @param errors Vector of local error estimates (typically from estimator_.GetLocalErrors()).
     * @param marked_elements Output array storing the indices of elements marked for refinement.
     * @return The number of elements marked for refinement.
     */
    int MarkElementsForRefinement(const Vector &errors, Array<int> &marked_elements);

    /**
     * @brief Applies refinement to the mesh and updates associated FE spaces and solution vectors.
     *
     * NOTE: This method assumes that the FE spaces and solution vectors associated
     * with the estimator are the ones that need updating after refinement.
     * The user needs to handle the re-assembly of operators after calling this.
     *
     * @param mesh The ParMesh to refine.
     * @param trial_fes The trial ParFiniteElementSpace to update.
     * @param test_fes The test ParFiniteElementSpace to update.
     * @param solution The solution BlockVector to update (interpolated to the new mesh).
     * @return The number of elements that were actually refined.
     */
    int Refine(ParMesh &mesh, ParFiniteElementSpace &trial_fes,
               ParFiniteElementSpace &test_fes, BlockVector &solution);


    /**
     * @brief Set the maximum level of hanging nodes allowed during refinement.
     * @param nc_limit 0 for conforming refinement, >0 for non-conforming.
     */
    void SetNCLimit(int nc_limit) { nc_limit_ = nc_limit; }

     /**
      * @brief Set the verbosity level.
      * @param level 0 (silent), 1 (basic info).
      */
    void SetVerbosity(int level) { verbose_ = level; }
};


/**
 * @brief Performs adaptive mesh derefinement (coarsening) based on DPG error estimates.
 */
class DPGAdaptiveDerefiner {
private:
    UltraweakDPGErrorEstimator &estimator_; // Reference to the estimator

    // Derefinement parameters
    double error_threshold_; // Elements with error below this are candidates
    int op_ = 1;             // Comparison operator (1 for '<', see DerefineByError)
    int nc_limit_ = 0;       // Max level of hanging nodes after derefinement

    // Verbosity
    int verbose_ = 0;

public:
    /**
     * @brief Constructor.
     * @param est Reference to the UltraweakDPGErrorEstimator.
     * @param error_threshold Error threshold below which elements are candidates for derefinement.
     * @param max_nc_level Maximum level of hanging nodes allowed after derefinement (default 0).
     */
    DPGAdaptiveDerefiner(UltraweakDPGErrorEstimator &est, double error_threshold,
                         int max_nc_level = 0);

    /**
     * @brief Applies derefinement to the mesh and updates associated FE spaces and solution vectors.
     *
     * NOTE: This method assumes that the FE spaces and solution vectors associated
     * with the estimator are the ones that need updating after derefinement.
     * The user needs to handle the re-assembly of operators after calling this.
     *
     * @param mesh The ParMesh to derefine.
     * @param trial_fes The trial ParFiniteElementSpace to update.
     * @param test_fes The test ParFiniteElementSpace to update.
     * @param solution The solution BlockVector to update (interpolated to the new mesh).
     * @return True if the mesh was changed, false otherwise.
     */
    bool Derefine(ParMesh &mesh, ParFiniteElementSpace &trial_fes,
                  ParFiniteElementSpace &test_fes, BlockVector &solution);


    /**
     * @brief Set the maximum level of hanging nodes allowed after derefinement.
     * @param nc_limit 0 for conforming refinement, >0 for non-conforming.
     */
    void SetNCLimit(int nc_limit) { nc_limit_ = nc_limit; }

    /**
     * @brief Set the verbosity level.
     * @param level 0 (silent), 1 (basic info).
     */
    void SetVerbosity(int level) { verbose_ = level; }
};


/**
 * @brief Standalone function to test the DPG error estimator on a simple problem (e.g., Poisson).
 *
 * This function will set up a basic DPG system for a known problem,
 * solve it, compute the error estimates, and perform basic validation checks.
 *
 * @param order Polynomial order for trial space.
 * @param dim Dimension of the problem (e.g., 2 for 2D).
 * @param ref_levels Number of initial uniform refinements for the mesh.
 * @param comm MPI communicator.
 * @return True if the basic tests pass, false otherwise.
 */
bool TestDPGErrorEstimator(int order, int dim, int ref_levels, MPI_Comm comm);


} // namespace mfem

#endif // DPG_ERROR_ESTIMATOR_HPP