#include "dpg-error-estimator.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace mfem {

// Implementation of UltraweakDPGErrorEstimator

UltraweakDPGErrorEstimator::UltraweakDPGErrorEstimator(
    FiniteElementSpace *trial_space, FiniteElementSpace *test_space,
    BilinearForm *ip_form, MixedBilinearForm *dpg_form, LinearForm *l_form,
    GridFunction *sol, int max_iterations, double relative_tolerance,
    int verbosity)
    : trial_fes(trial_space), enriched_test_fes(test_space),
      test_inner_product(ip_form), dpg(dpg_form), dpg_linear_form(l_form),
      solution(sol), max_it(max_iterations), rel_tol(relative_tolerance),
      verbose(verbosity) {
  // Validate input parameters
  if (!trial_space) {
    throw std::invalid_argument(
        "UltraweakDPGErrorEstimator: trial_space is null");
  }
  if (!test_space) {
    throw std::invalid_argument(
        "UltraweakDPGErrorEstimator: test_space is null");
  }
  if (!sol) {
    throw std::invalid_argument("UltraweakDPGErrorEstimator: solution is null");
  }

  // No need to validate ip_form, dpg_form, l_form as they can be null in some
  // cases

  // Check mesh compatibility
  if (trial_space->GetMesh() != test_space->GetMesh()) {
    throw std::invalid_argument("UltraweakDPGErrorEstimator: trial and test "
                                "spaces must use the same mesh");
  }

  // Initialize local errors array if we have a valid mesh
  if (trial_space && trial_space->GetMesh()) {
    local_errors.SetSize(trial_space->GetMesh()->GetNE());
    local_errors = 0.0;
  }

  if (verbose > 0) {
    std::cout << "Created UltraweakDPGErrorEstimator with:" << std::endl
              << "  Trial space size: " << trial_fes->GetVSize() << std::endl
              << "  Test space size: " << enriched_test_fes->GetVSize()
              << std::endl
              << "  Number of elements: " << local_errors.Size() << std::endl;
  }
}

UltraweakDPGErrorEstimator::~UltraweakDPGErrorEstimator() {
  // We don't own any of the pointers, so no need to delete them
}

const Vector &UltraweakDPGErrorEstimator::GetLocalErrors() {
  // Check if we need to compute the errors
  if (local_errors.Size() == 0 && trial_fes && trial_fes->GetMesh()) {
    local_errors.SetSize(trial_fes->GetMesh()->GetNE());
    ComputeEstimator(local_errors);
  } else if (local_errors.Size() != trial_fes->GetMesh()->GetNE()) {
    // If the mesh changed, we need to resize and recompute
    local_errors.SetSize(trial_fes->GetMesh()->GetNE());
    ComputeEstimator(local_errors);
  }

  return local_errors;
}

void UltraweakDPGErrorEstimator::Reset() {
  if (trial_fes && trial_fes->GetMesh()) {
    local_errors.SetSize(trial_fes->GetMesh()->GetNE());
    local_errors = 0.0;
  }

  if (verbose > 1) {
    std::cout << "Reset UltraweakDPGErrorEstimator" << std::endl;
  }
}

double UltraweakDPGErrorEstimator::GetLocalError(int elem_idx) {
  // Check if the element index is valid
  Mesh *mesh = trial_fes->GetMesh();
  if (!mesh) {
    throw std::runtime_error("UltraweakDPGErrorEstimator: mesh is null");
  }

  if (elem_idx < 0 || elem_idx >= mesh->GetNE()) {
    throw std::out_of_range(
        "UltraweakDPGErrorEstimator: element index out of range");
  }

  // Get the DOFs for this element in the trial and test spaces
  Array<int> trial_dofs, test_dofs;
  trial_fes->GetElementDofs(elem_idx, trial_dofs);
  enriched_test_fes->GetElementDofs(elem_idx, test_dofs);

  // Extract the local solution
  Vector local_sol(trial_dofs.Size());
  if (!solution) {
    throw std::runtime_error("UltraweakDPGErrorEstimator: solution is null");
  }
  solution->GetSubVector(trial_dofs, local_sol);

  // Make sure bilinear forms are assembled
  if (dpg) {
    dpg->Assemble();
  }
  if (dpg_linear_form) {
    dpg_linear_form->Assemble();
  }
  if (test_inner_product) {
    test_inner_product->Assemble();
  }

  // Get the local right-hand side vector
  Vector l_elem(test_dofs.Size());
  l_elem = 0.0;

  if (dpg_linear_form) {
    dpg_linear_form->GetSubVector(test_dofs, l_elem);
  }

  // Compute B*u (the action of the operator on the solution)
  Vector Bu(test_dofs.Size());
  Bu = 0.0;

  if (dpg) {
    SparseMatrix &mat = dpg->SpMat();

    // Multiply the local submatrix by the local solution
    for (int i = 0; i < test_dofs.Size(); i++) {
      int row = test_dofs[i];
      if (row >= mat.Height())
        continue;

      for (int j = 0; j < trial_dofs.Size(); j++) {
        int col = trial_dofs[j];
        if (col >= mat.Width())
          continue;

        double val = mat.Elem(row, col);
        if (val != 0.0) {
          Bu(i) += val * local_sol(j);
        }
      }
    }
  }

  // Compute the residual: r = l - B*u
  Vector res(test_dofs.Size());
  res = l_elem;
  res -= Bu;

  // Extract the test space inner product matrix for this element
  DenseMatrix G_elem(test_dofs.Size(), test_dofs.Size());
  G_elem = 0.0;

  if (test_inner_product) {
    SparseMatrix &G_mat = test_inner_product->SpMat();

    for (int i = 0; i < test_dofs.Size(); i++) {
      int row = test_dofs[i];
      if (row >= G_mat.Height())
        continue;

      for (int j = 0; j < test_dofs.Size(); j++) {
        int col = test_dofs[j];
        if (col >= G_mat.Width())
          continue;

        double val = G_mat.Elem(row, col);
        if (val != 0.0) {
          G_elem(i, j) = val;
        }
      }
    }
  }

  // Add regularization to diagonal if necessary
  double min_diag = std::numeric_limits<double>::max();
  for (int i = 0; i < G_elem.Height(); i++) {
    if (G_elem(i, i) < min_diag) {
      min_diag = G_elem(i, i);
    }

    if (G_elem(i, i) <= 1e-14) {
      // Add regularization for numerical stability
      G_elem(i, i) = 1e-10;
    }
  }

  // If all diagonal entries are very small, scale the matrix
  if (min_diag < 1e-12 && min_diag > 0) {
    double scale = 1.0 / min_diag;
    G_elem *= scale;
  }

  // Solve G * err_repr = res for the error representation
  Vector err_repr(test_dofs.Size());

  // Use the dense matrix inverse to solve the system
  DenseMatrixInverse G_inv(G_elem);

  try {
    G_inv.Mult(res, err_repr);
  } catch (const std::exception &e) {
    std::cerr << "Warning: Exception in element " << elem_idx
              << " when computing error: " << e.what() << std::endl;

    // If inversion fails, use diagonal approximation as fallback
    for (int i = 0; i < G_elem.Height(); i++) {
      err_repr(i) = res(i) / G_elem(i, i);
    }
  }

  // Compute the error indicator: sqrt(res^T * err_repr)
  double error_indicator = sqrt(std::abs(res * err_repr));

  // Check for NaN or inf
  if (!std::isfinite(error_indicator)) {
    std::cerr << "Warning: Non-finite error indicator in element " << elem_idx
              << ": " << error_indicator << std::endl;
    return 0.0; // Return zero as fallback
  }

  if (verbose > 1) {
    std::cout << "Element " << elem_idx << " error: " << error_indicator
              << std::endl;
  }

  return error_indicator;
}

void UltraweakDPGErrorEstimator::ComputeEstimator(Vector &errors) {
  int ne = trial_fes->GetMesh()->GetNE();
  errors.SetSize(ne);

  if (verbose > 0) {
    std::cout << "Computing error estimates for " << ne << " elements..."
              << std::endl;
  }

  double total_error = 0.0;
  double max_error = 0.0;
  int max_error_elem = -1;

  // Compute error for each element
  for (int i = 0; i < ne; i++) {
    errors(i) = GetLocalError(i);

    total_error += errors(i) * errors(i);

    if (errors(i) > max_error) {
      max_error = errors(i);
      max_error_elem = i;
    }
  }

  total_error = sqrt(total_error);

  if (verbose > 0) {
    std::cout << "Error estimation complete:" << std::endl
              << "  Total error: " << total_error << std::endl
              << "  Maximum error: " << max_error << " (element "
              << max_error_elem << ")" << std::endl;
  }

  // Update the stored local errors
  local_errors = errors;
}

// Implementation of DPGAdaptiveRefiner

DPGAdaptiveRefiner::DPGAdaptiveRefiner(UltraweakDPGErrorEstimator &est,
                                       double err_fraction, double threshold,
                                       int strategy_type)
    : estimator(est), total_error_fraction(err_fraction),
      threshold_factor(threshold), strategy(strategy_type) {
  // Validate input parameters
  if (err_fraction <= 0.0 || err_fraction > 1.0) {
    throw std::invalid_argument(
        "DPGAdaptiveRefiner: err_fraction must be in (0,1]");
  }
  if (threshold <= 0.0 || threshold > 1.0) {
    throw std::invalid_argument(
        "DPGAdaptiveRefiner: threshold must be in (0,1]");
  }
  if (strategy_type < 0 || strategy_type > 1) {
    throw std::invalid_argument("DPGAdaptiveRefiner: strategy must be 0 or 1");
  }
}

int DPGAdaptiveRefiner::MarkElementsForRefinement(const Vector &errors,
                                                  Array<int> &marked_elements) {
  marked_elements.DeleteAll();

  // Check for empty errors vector
  if (errors.Size() == 0) {
    return 0;
  }

  double max_error = errors.Max();

  // If all errors are very close to zero, don't mark anything
  if (max_error < 1e-14) {
    return 0;
  }

  double threshold = threshold_factor * max_error;
  int size = errors.Size();

  // Maximum strategy (strategy 0):
  // Mark all elements with error > threshold*max_error
  if (strategy == 0) {
    for (int i = 0; i < size; i++) {
      if (errors(i) > threshold) {
        marked_elements.Append(i);
      }
    }
  }
  // Bulk/Dörfler marking strategy (strategy 1):
  // Mark a minimal subset of elements whose combined
  // squared error exceeds the specified fraction of total squared error
  else if (strategy == 1) {
    // Create pairs of (error, index) for sorting
    std::vector<std::pair<double, int>> error_indices(size);
    for (int i = 0; i < size; i++) {
      error_indices[i] = std::make_pair(errors(i), i);
    }

    // Sort in descending order of error
    std::sort(
        error_indices.begin(), error_indices.end(),
        [](const std::pair<double, int> &a, const std::pair<double, int> &b) {
          return a.first > b.first;
        });

    // Calculate total squared error
    double total_error_sq = 0.0;
    for (int i = 0; i < size; i++) {
      total_error_sq += errors(i) * errors(i);
    }

    // Target squared error to achieve
    double target_error_sq = total_error_fraction * total_error_sq;

    // Accumulate squared error until we reach the target
    double accumulated_error_sq = 0.0;
    for (int i = 0; i < size; i++) {
      double error_value = error_indices[i].first;
      int elem_idx = error_indices[i].second;

      // Skip elements with very small error
      if (error_value < threshold) {
        break;
      }

      marked_elements.Append(elem_idx);
      accumulated_error_sq += error_value * error_value;

      if (accumulated_error_sq >= target_error_sq) {
        break;
      }
    }
  }

  return marked_elements.Size();
}

int DPGAdaptiveRefiner::Refine(Mesh &mesh, GridFunction &solution) {
  // Get the error estimates from the estimator
  const Vector &errors = estimator.GetLocalErrors();

  // Create array to store marked elements
  Array<int> marked_elements;

  // Mark elements for refinement
  int num_marked = MarkElementsForRefinement(errors, marked_elements);

  // If any elements were marked, refine the mesh
  if (num_marked > 0) {
    // Check if the mesh supports nonconforming refinement
    if (!mesh.Nonconforming()) {
      mesh.EnsureNCMesh(true);
    }

    // Apply refinement
    mesh.GeneralRefinement(marked_elements);

    // Update the solution
    solution.Update();
  }

  return num_marked;
}

// Implementation of DPGAdaptiveDerefiner

DPGAdaptiveDerefiner::DPGAdaptiveDerefiner(UltraweakDPGErrorEstimator &est,
                                           double error_threshold,
                                           int max_nc_level)
    : estimator(est), threshold(error_threshold), nc_limit(max_nc_level) {
  // Validate input parameters
  if (error_threshold <= 0.0) {
    throw std::invalid_argument(
        "DPGAdaptiveDerefiner: error_threshold must be positive");
  }
  if (max_nc_level < 0) {
    throw std::invalid_argument(
        "DPGAdaptiveDerefiner: max_nc_level must be non-negative");
  }
}

int DPGAdaptiveDerefiner::MarkElementsForDerefinement(
    const Vector &errors, Array<int> &marked_elements) {
  marked_elements.DeleteAll();

  // Get the mesh from the estimator
  Mesh *mesh = estimator.GetTrialSpace()->GetMesh();
  if (!mesh) {
    return 0;
  }

  // Check if the mesh supports derefinement
  if (!mesh->Nonconforming()) {
    return 0;
  }

  // Loop through all elements and mark those with small errors
  for (int i = 0; i < errors.Size(); i++) {
    if (errors(i) < threshold) {
      // Mark element for possible derefinement
      marked_elements.Append(i);
    }
  }

  return marked_elements.Size();
}

bool DPGAdaptiveDerefiner::Apply(Mesh &mesh, GridFunction &solution) {
  // Get the error estimates from the estimator
  const Vector &errors = estimator.GetLocalErrors();

  // Create array to store marked elements
  Array<int> marked_elements;

  // Mark elements for derefinement
  int num_marked = MarkElementsForDerefinement(errors, marked_elements);

  if (num_marked == 0) {
    return false;
  }

  // Apply derefinement
  mesh.DerefineByError(errors, threshold, nc_limit);

  // Update the solution
  solution.Update();

  return true;
}

// Implementation of test function
bool TestDPGErrorEstimator(int order, int dim, int ref_levels) {
  std::cout << "\n=== Testing DPG Error Estimator ===" << std::endl;
  std::cout << "  Order: " << order << ", Dimension: " << dim
            << ", Refinement levels: " << ref_levels << std::endl;

  bool all_tests_passed = true;

  try {
    // 1. Create a simple mesh
    Mesh *mesh = nullptr;

    if (dim == 1) {
      // Create 1D mesh using the more modern approach
      mesh = new Mesh(Mesh::MakeCartesian1D(10, 1.0));
    } else if (dim == 2) {
      // Create 2D mesh using the more modern approach
      mesh = new Mesh(Mesh::MakeCartesian2D(10, 10, Element::QUADRILATERAL,
                                            true, 1.0, 1.0));
    } else if (dim == 3) {
      // Create 3D mesh using the more modern approach
      mesh = new Mesh(Mesh::MakeCartesian3D(4, 4, 4, Element::HEXAHEDRON, true,
                                            1.0, 1.0, 1.0));
    } else {
      std::cerr << "Invalid dimension: " << dim << std::endl;
      return false;
    }

    // Perform uniform refinement
    for (int l = 0; l < ref_levels; l++) {
      mesh->UniformRefinement();
    }

    // Make sure the mesh is nonconforming
    mesh->EnsureNCMesh(true);

    std::cout << "  Created mesh with " << mesh->GetNE() << " elements."
              << std::endl;

    // 2. Define finite element spaces for trial and test functions
    H1_FECollection trial_fec(order, dim);
    H1_FECollection test_fec(order + 1, dim); // Enriched test space

    FiniteElementSpace trial_space(mesh, &trial_fec);
    FiniteElementSpace test_space(mesh, &test_fec);

    std::cout << "  Created trial space with " << trial_space.GetVSize()
              << " DOFs and test space with " << test_space.GetVSize()
              << " DOFs." << std::endl;

    // 3. Create a simple source function and solution
    ConstantCoefficient one(1.0);
    FunctionCoefficient f([](const Vector &x) {
      double r2 = 0.0;
      for (int i = 0; i < x.Size(); i++) {
        r2 += (x(i) - 0.5) * (x(i) - 0.5);
      }
      return exp(-10.0 * r2);
    });

    // 4. Set up the bilinear and linear forms
    BilinearForm a(&trial_space);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));
    a.AddDomainIntegrator(new MassIntegrator(one));
    a.Assemble();
    a.Finalize();

    LinearForm b(&trial_space);
    b.AddDomainIntegrator(new DomainLFIntegrator(f));
    b.Assemble();

    // 5. Define and solve the linear system
    GridFunction x(&trial_space);
    x = 0.0;

    Array<int> ess_tdof_list;

    SparseMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

    // Use a simple solver
    GSSmoother smoother(A);
    PCG(A, smoother, B, X, 0, 1000, 1e-12, 0.0);

    a.RecoverFEMSolution(X, b, x);

    std::cout << "  Solved the linear system." << std::endl;

    // 6. Set up the test space inner product
    BilinearForm m_test(&test_space);
    m_test.AddDomainIntegrator(new DiffusionIntegrator(one));
    m_test.AddDomainIntegrator(new MassIntegrator(one));
    m_test.Assemble();
    m_test.Finalize();

    // 7. Set up the mixed DPG bilinear form
    MixedBilinearForm dpg(&trial_space, &test_space);
    dpg.AddDomainIntegrator(new DiffusionIntegrator(one));
    dpg.AddDomainIntegrator(new MassIntegrator(one));
    dpg.Assemble();
    dpg.Finalize();

    // 8. Create an error estimator
    UltraweakDPGErrorEstimator estimator(&trial_space, &test_space, &m_test,
                                         &dpg, &b, &x, 100, 1e-10, 1);

    // 9. Compute the error estimates
    const Vector &errors = estimator.GetLocalErrors();

    // Verify that we have the correct number of error estimates
    if (errors.Size() != mesh->GetNE()) {
      std::cerr << "  ERROR: Number of error estimates (" << errors.Size()
                << ") does not match number of elements (" << mesh->GetNE()
                << ")" << std::endl;
      all_tests_passed = false;
    } else {
      std::cout << "  Computed error estimates for " << errors.Size()
                << " elements." << std::endl;
    }

    // Check that all error estimates are finite and non-negative
    bool all_errors_valid = true;
    for (int i = 0; i < errors.Size(); i++) {
      if (!std::isfinite(errors(i)) || errors(i) < 0.0) {
        std::cerr << "  ERROR: Invalid error estimate at element " << i << ": "
                  << errors(i) << std::endl;
        all_errors_valid = false;
        all_tests_passed = false;
        break;
      }
    }

    if (all_errors_valid) {
      std::cout << "  All error estimates are valid." << std::endl;
    }

    // 10. Test the adaptive refiner
    DPGAdaptiveRefiner refiner(estimator, 0.7, 0.5, 0);
    Array<int> marked_elements;
    int num_marked = refiner.MarkElementsForRefinement(errors, marked_elements);

    std::cout << "  Marked " << num_marked << " elements for refinement."
              << std::endl;

    // Verify the number of marked elements is sensible
    if (num_marked == 0 || num_marked > mesh->GetNE() / 2) {
      std::cout << "  WARNING: Number of marked elements (" << num_marked
                << ") is unexpected." << std::endl;
    }

    // 11. Test the adaptive derefiner
    DPGAdaptiveDerefiner derefiner(estimator, 1e-4, 3);
    int num_derefinable =
        derefiner.MarkElementsForDerefinement(errors, marked_elements);

    std::cout << "  Found " << num_derefinable
              << " elements eligible for derefinement." << std::endl;
    delete mesh;

    std::cout << "=== Test " << (all_tests_passed ? "PASSED" : "FAILED")
              << " ===" << std::endl;

  } catch (const std::exception &e) {
    std::cerr << "ERROR: Exception in TestDPGErrorEstimator: " << e.what()
              << std::endl;
    all_tests_passed = false;
  }

  return all_tests_passed;
}

} // namespace mfem