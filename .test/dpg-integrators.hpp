#ifndef DPG_INTEGRATORS_HPP
#define DPG_INTEGRATORS_HPP

#include "mfem.hpp"

namespace mfem {

// Base class for Ultraweak DPG integrators
class UltraweakDPGIntegrator : public BilinearFormIntegrator {
protected:
  Coefficient *Q;

public:
  UltraweakDPGIntegrator(Coefficient *q = nullptr);
  virtual ~UltraweakDPGIntegrator();

  virtual void AssembleElementMatrix(const FiniteElement &el,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat) override {
    mfem_error("UltraweakDPGIntegrator::AssembleElementMatrix should not be "
               "called directly");
  }
};

// Domain integrator for first derivative terms in ultraweak form
class UltraweakDerivativeIntegrator : public UltraweakDPGIntegrator {
private:
  int di;

public:
  UltraweakDerivativeIntegrator(int dir, Coefficient *q = nullptr);

  virtual void AssembleElementMatrix(const FiniteElement &trial_fe,
                                     const FiniteElement &test_fe,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat);

  virtual void AssembleElementMatrix(const FiniteElement &el,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat) override {
    mfem_error("UltraweakDerivativeIntegrator: Use "
               "AssembleElementMatrix(trial, test, ...) version");
  }
};

// Integrator for second derivative terms in ultraweak form
class UltraweakSecondDerivativeIntegrator : public UltraweakDPGIntegrator {
private:
  int di, dj;

public:
  UltraweakSecondDerivativeIntegrator(int dir1, int dir2,
                                      Coefficient *q = nullptr);

  virtual void AssembleElementMatrix(const FiniteElement &trial_fe,
                                     const FiniteElement &test_fe,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat);

  virtual void AssembleElementMatrix(const FiniteElement &el,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat) override {
    mfem_error("UltraweakSecondDerivativeIntegrator: Use "
               "AssembleElementMatrix(trial, test, ...) version");
  }
};

// Face integrator for the ultraweak DPG formulation to handle trace terms
class UltraweakTraceIntegrator : public BilinearFormIntegrator {
private:
  int di;
  double alpha;

public:
  UltraweakTraceIntegrator(int dir, double alpha_ = -1.0);

  virtual void AssembleFaceMatrix(const FiniteElement &trial_fe1,
                                  const FiniteElement &trial_fe2,
                                  const FiniteElement &test_fe1,
                                  const FiniteElement &test_fe2,
                                  FaceElementTransformations &Trans,
                                  DenseMatrix &elmat);

  virtual void AssembleFaceMatrix(const FiniteElement &el1,
                                  const FiniteElement &el2,
                                  FaceElementTransformations &Trans,
                                  DenseMatrix &elmat) override {
    mfem_error("UltraweakTraceIntegrator::AssembleFaceMatrix(el1,el2) not "
               "implemented");
  }

  virtual void AssembleFaceMatrix(const FiniteElement &trial_fe,
                                  const FiniteElement &test_fe1,
                                  const FiniteElement &test_fe2,
                                  FaceElementTransformations &Trans,
                                  DenseMatrix &elmat) override {
    mfem_error("UltraweakTraceIntegrator::AssembleFaceMatrix(trial,test1,test2)"
               " not implemented");
  }
};

// Boundary integrator for the ultraweak DPG formulation
class UltraweakBoundaryIntegrator : public BilinearFormIntegrator {
private:
  int di;
  Coefficient *bdr_coef;

public:
  UltraweakBoundaryIntegrator(int dir, Coefficient *bdr_coef_ = nullptr);

  virtual void AssembleFaceMatrix(const FiniteElement &trial_fe,
                                  const FiniteElement &test_fe,
                                  FaceElementTransformations &Trans,
                                  DenseMatrix &elmat);
};

// Boundary linear form integrator for the ultraweak DPG formulation
class UltraweakBoundaryLFIntegrator : public LinearFormIntegrator {
private:
  int di;
  FunctionCoefficient *bdr_func;

public:
  UltraweakBoundaryLFIntegrator(int dir,
                                FunctionCoefficient *bdr_func_ = nullptr);

  virtual void AssembleFaceVector(const FiniteElement &el,
                                  FaceElementTransformations &Trans,
                                  Vector &elvect);
};

// Integrator for jump terms across element faces
class UltraweakJumpIntegrator : public BilinearFormIntegrator {
private:
  int direction;

public:
  UltraweakJumpIntegrator(int dir);

  virtual void AssembleFaceMatrix(const FiniteElement &trial_fe1,
                                  const FiniteElement &trial_fe2,
                                  const FiniteElement &test_fe,
                                  FaceElementTransformations &Trans,
                                  DenseMatrix &elmat);

  virtual void AssembleFaceMatrix(const FiniteElement &el1,
                                  const FiniteElement &el2,
                                  FaceElementTransformations &Trans,
                                  DenseMatrix &elmat) override {
    mfem_error("UltraweakJumpIntegrator::AssembleFaceMatrix(el1,el2) not "
               "implemented. Use version with explicit test_fe.");
  }

  virtual void AssembleFaceMatrix(const FiniteElement &trial_fe1,
                                  const FiniteElement &trial_fe2,
                                  const FiniteElement &test_fe1,
                                  const FiniteElement &test_fe2,
                                  FaceElementTransformations &Trans,
                                  DenseMatrix &elmat);
};

// Integrator for average terms across element faces
class UltraweakAverageIntegrator : public BilinearFormIntegrator {
private:
  int direction;
  double scaling;

public:
  UltraweakAverageIntegrator(int dir, double scale = 1.0);

  virtual void AssembleFaceMatrix(const FiniteElement &trial_fe1,
                                  const FiniteElement &trial_fe2,
                                  const FiniteElement &test_fe1,
                                  const FiniteElement &test_fe2,
                                  FaceElementTransformations &Trans,
                                  DenseMatrix &elmat);

  virtual void AssembleFaceMatrix(const FiniteElement &el1,
                                  const FiniteElement &el2,
                                  FaceElementTransformations &Trans,
                                  DenseMatrix &elmat) override {
    mfem_error("UltraweakAverageIntegrator::AssembleFaceMatrix(el1,el2) not "
               "implemented");
  }

  virtual void AssembleFaceMatrix(const FiniteElement &trial_fe,
                                  const FiniteElement &test_fe1,
                                  const FiniteElement &test_fe2,
                                  FaceElementTransformations &Trans,
                                  DenseMatrix &elmat) override {
    mfem_error("UltraweakAverageIntegrator::AssembleFaceMatrix(trial,test1,"
               "test2) not implemented");
  }
};

// Integrator for mixed derivative terms in Heston-SLV
class UltraweakMixedDerivativeIntegrator : public UltraweakDPGIntegrator {
public:
  UltraweakMixedDerivativeIntegrator(Coefficient *q = nullptr);

  virtual void AssembleElementMatrix(const FiniteElement &trial_fe,
                                     const FiniteElement &test_fe,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat);

  virtual void AssembleElementMatrix(const FiniteElement &el,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat) override {
    mfem_error("UltraweakMixedDerivativeIntegrator: Use "
               "AssembleElementMatrix(trial, test, ...) version");
  }
};

/**
 * @brief Class for computing optimal test functions in the DPG method
 *
 * This class computes the optimal test functions for the DPG method by solving
 * local problems.
 */
class OptimalTestFunctionComputer {
private:
  /// Trial finite element space
  FiniteElementSpace *trial_fes;

  /// Test finite element space
  FiniteElementSpace *test_fes;

  /// DPG bilinear form
  MixedBilinearForm *dpg_form;

  /// Test space inner product
  BilinearForm *test_ip;

public:
  /**
   * @brief Constructor for OptimalTestFunctionComputer
   *
   * @param tr_fes Trial finite element space
   * @param te_fes Test finite element space
   * @param dpg_form DPG bilinear form
   * @param ip Test space inner product
   * @throws std::invalid_argument If any required pointer is null
   */
  OptimalTestFunctionComputer(FiniteElementSpace *tr_fes,
                              FiniteElementSpace *te_fes,
                              MixedBilinearForm *dpg_form, BilinearForm *ip);

  /**
   * @brief Compute the optimal test function for a trial basis function
   *
   * @param elem_idx Element index
   * @param trial_dof_idx Local trial DOF index
   * @param optimal_test Vector to store the computed optimal test function
   * @throws std::out_of_range If indices are out of range
   */
  void ComputeOptimalTestFunction(int elem_idx, int trial_dof_idx,
                                  Vector &optimal_test);

  /**
   * @brief Compute all optimal test functions for an element
   *
   * @param elem_idx Element index
   * @param optimal_test_functions Matrix to store the computed optimal test
   * functions
   * @throws std::out_of_range If elem_idx is out of range
   */
  void ComputeAllOptimalTestFunctions(int elem_idx,
                                      DenseMatrix &optimal_test_functions);
};

// Drift integrator for handling S and V drift terms in Heston
class UltraweakDriftIntegrator : public UltraweakDPGIntegrator {
public:
  UltraweakDriftIntegrator(Coefficient *q = nullptr);

  virtual void AssembleElementMatrix(const FiniteElement &trial_fe,
                                     const FiniteElement &test_fe,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat);

  virtual void AssembleElementMatrix(const FiniteElement &el,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat) override {
    mfem_error("UltraweakDriftIntegrator: Use AssembleElementMatrix(trial, "
               "test, ...) version");
  }
};

// Reaction integrator for handling reaction terms in Heston (e.g., ru)
class UltraweakReactionIntegrator : public UltraweakDPGIntegrator {
public:
  UltraweakReactionIntegrator(Coefficient *q = nullptr);

  virtual void AssembleElementMatrix(const FiniteElement &trial_fe,
                                     const FiniteElement &test_fe,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat);

  virtual void AssembleElementMatrix(const FiniteElement &el,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat) override {
    mfem_error("UltraweakReactionIntegrator: Use AssembleElementMatrix(trial, "
               "test, ...) version");
  }
};

// Test inner product integrator for DPG method
class TestInnerProductIntegrator : public BilinearFormIntegrator {
private:
  double alpha_diff;
  double alpha_mass;

public:
  TestInnerProductIntegrator(double diff_coef = 1.0, double mass_coef = 1.0);

  virtual void AssembleElementMatrix(const FiniteElement &el,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat);
};

} // namespace mfem

#endif // DPG_INTEGRATORS_HPP