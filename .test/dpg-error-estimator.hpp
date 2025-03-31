#ifndef DPG_ERROR_ESTIMATOR_HPP
#define DPG_ERROR_ESTIMATOR_HPP

#include "mfem.hpp"
#include <memory>
#include <stdexcept>
#include <vector>

namespace mfem {

// Forward declaration for OptimalTestFunctionComputer, which is now defined
// only in dpg-integrators.hpp
class OptimalTestFunctionComputer;

/**
 * @brief Error estimator for the ultraweak DPG method
 *
 * This class implements the residual-based error estimator for the ultraweak
 * formulation of the DPG method. It computes the norm of the residual in the
 * appropriate test space for each element.
 */
class UltraweakDPGErrorEstimator {
private:
  /// Finite element space for the trial variables
  FiniteElementSpace *trial_fes;

  /// Enriched finite element space for the test variables
  FiniteElementSpace *enriched_test_fes;

  /// Bilinear form for the test space inner product
  BilinearForm *test_inner_product;

  /// Bilinear form for the DPG formulation
  MixedBilinearForm *dpg;

  /// Linear form for the DPG formulation
  LinearForm *dpg_linear_form;

  /// Solution grid function
  GridFunction *solution;

  /// Maximum number of iterations for local solvers
  int max_it;

  /// Relative tolerance for local solvers
  double rel_tol;

  /// Vector of local error estimates
  Vector local_errors;

  /// Verbosity level (0=silent, 1=minimal, 2=detailed)
  int verbose;

public:
  /**
   * @brief Constructor for UltraweakDPGErrorEstimator
   *
   * @param trial_space Trial finite element space
   * @param test_space Enriched test finite element space
   * @param ip_form Bilinear form for test space inner product
   * @param dpg_form Bilinear form for the DPG formulation
   * @param l_form Linear form for the DPG formulation
   * @param sol Solution grid function
   * @param max_iterations Maximum iterations for local solvers (default: 100)
   * @param relative_tolerance Relative tolerance for local solvers (default:
   * 1e-10)
   * @param verbosity Verbosity level (default: 0)
   * @throws std::invalid_argument If any required pointer is null
   */
  UltraweakDPGErrorEstimator(FiniteElementSpace *trial_space,
                             FiniteElementSpace *test_space,
                             BilinearForm *ip_form, MixedBilinearForm *dpg_form,
                             LinearForm *l_form, GridFunction *sol,
                             int max_iterations = 100,
                             double relative_tolerance = 1e-10,
                             int verbosity = 0);

  /**
   * @brief Destructor
   */
  virtual ~UltraweakDPGErrorEstimator();

  /**
   * @brief Get the local error estimates for all elements
   *
   * Computes the error estimate for each element if not already computed.
   *
   * @return const Vector& Vector containing error estimates for each element
   */
  virtual const Vector &GetLocalErrors();

  /**
   * @brief Reset the error estimator
   *
   * Clears cached error estimates, forcing recomputation on next call.
   */
  virtual void Reset();

  /**
   * @brief Compute the error estimate for a single element
   *
   * @param elem_idx Index of the element
   * @return double Error estimate for the element
   * @throws std::out_of_range If elem_idx is out of range
   */
  double GetLocalError(int elem_idx);

  /**
   * @brief Compute error estimates for all elements
   *
   * @param errors Vector to store the error estimates
   */
  void ComputeEstimator(Vector &errors);

  /**
   * @brief Set verbosity level
   *
   * @param level Verbosity level (0=silent, 1=minimal, 2=detailed)
   */
  void SetVerboseLevel(int level) { verbose = level; }

  /**
   * @brief Get the trial finite element space
   *
   * @return FiniteElementSpace* Trial space pointer
   */
  FiniteElementSpace *GetTrialSpace() const { return trial_fes; }

  /**
   * @brief Get the test finite element space
   *
   * @return FiniteElementSpace* Test space pointer
   */
  FiniteElementSpace *GetTestSpace() const { return enriched_test_fes; }
};

/**
 * @brief Class for adaptive mesh refinement based on error estimates
 *
 * This class uses error estimates from an UltraweakDPGErrorEstimator to mark
 * elements for refinement based on various strategies.
 */
class DPGAdaptiveRefiner {
private:
  /// Reference to the error estimator
  UltraweakDPGErrorEstimator &estimator;

  /// Fraction of total error to be reduced (for bulk marking)
  double total_error_fraction;

  /// Threshold factor relative to maximum error (for maximum strategy)
  double threshold_factor;

  /// Refinement strategy (0=maximum, 1=bulk/Dörfler)
  int strategy;

public:
  /**
   * @brief Constructor for DPGAdaptiveRefiner
   *
   * @param est Reference to an UltraweakDPGErrorEstimator
   * @param err_fraction Fraction of total error for bulk marking (default: 0.7)
   * @param threshold Threshold factor for maximum strategy (default: 0.7)
   * @param strategy_type Refinement strategy (0=maximum, 1=bulk/Dörfler)
   * (default: 0)
   */
  DPGAdaptiveRefiner(UltraweakDPGErrorEstimator &est, double err_fraction = 0.7,
                     double threshold = 0.7, int strategy_type = 0);

  /**
   * @brief Mark elements for refinement
   *
   * @param errors Vector of error estimates for each element
   * @param marked_elements Array to store indices of marked elements
   * @return int Number of marked elements
   */
  int MarkElementsForRefinement(const Vector &errors,
                                Array<int> &marked_elements);

  /**
   * @brief Refine the mesh based on error estimates
   *
   * @param mesh Mesh to be refined
   * @param solution Grid function to be updated after refinement
   * @return int Number of elements marked for refinement
   */
  int Refine(Mesh &mesh, GridFunction &solution);

  /**
   * @brief Set the refinement strategy
   *
   * @param strategy_type Refinement strategy (0=maximum, 1=bulk/Dörfler)
   */
  void SetStrategy(int strategy_type) { strategy = strategy_type; }

  /**
   * @brief Reset the refiner
   *
   * Clears any cached data, forcing recomputation on next call.
   */
  void Reset() { estimator.Reset(); }
};

/**
 * @brief Class for adaptive mesh derefinement based on error estimates
 *
 * This class uses error estimates to identify elements that can be coarsened.
 */
class DPGAdaptiveDerefiner {
private:
  /// Reference to the error estimator
  UltraweakDPGErrorEstimator &estimator;

  /// Error threshold for derefinement
  double threshold;

  /// Maximum level of hanging nodes
  int nc_limit;

public:
  /**
   * @brief Constructor for DPGAdaptiveDerefiner
   *
   * @param est Reference to an UltraweakDPGErrorEstimator
   * @param error_threshold Error threshold for derefinement (default: 0.1)
   * @param max_nc_level Maximum level of hanging nodes (default: 3)
   */
  DPGAdaptiveDerefiner(UltraweakDPGErrorEstimator &est,
                       double error_threshold = 0.1, int max_nc_level = 3);

  /**
   * @brief Mark elements for derefinement
   *
   * @param errors Vector of error estimates for each element
   * @param marked_elements Array to store indices of marked elements
   * @return int Number of marked elements
   */
  int MarkElementsForDerefinement(const Vector &errors,
                                  Array<int> &marked_elements);

  /**
   * @brief Derefine the mesh based on error estimates
   *
   * @param mesh Mesh to be derefined
   * @param solution Grid function to be updated after derefinement
   * @return bool True if any elements were derefined
   */
  bool Apply(Mesh &mesh, GridFunction &solution);

  /**
   * @brief Set the error threshold
   *
   * @param error_threshold New threshold value
   */
  void SetThreshold(double error_threshold) { threshold = error_threshold; }

  /**
   * @brief Set the maximum level of hanging nodes
   *
   * @param max_nc_level Maximum level value
   */
  void SetNCLimit(int max_nc_level) { nc_limit = max_nc_level; }

  /**
   * @brief Reset the derefiner
   *
   * Clears any cached data, forcing recomputation on next call.
   */
  void Reset() { estimator.Reset(); }
};

/**
 * @brief Function to test the DPG error estimator functionality
 *
 * This function creates a simple test case and verifies the error estimator.
 *
 * @param order Polynomial order for the test
 * @param dim Dimension of the test problem (1, 2, or 3)
 * @param ref_levels Number of uniform refinement levels
 * @return bool True if all tests pass, false otherwise
 */
bool TestDPGErrorEstimator(int order = 2, int dim = 2, int ref_levels = 1);

} // namespace mfem

#endif // DPG_ERROR_ESTIMATOR_HPP