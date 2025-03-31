#ifndef HESTON_SLV_DPG_HPP
#define HESTON_SLV_DPG_HPP

#include "mfem.hpp"
#include <map>
#include <memory>
#include <string>
#include <vector>


namespace mfem {

class HestonSLVParameters {
public:
  double r;     // Risk-free interest rate
  double q;     // Dividend yield
  double kappa; // Mean reversion rate
  double theta; // Long-term variance level
  double xi;    // Volatility of volatility
  double rho;   // Correlation between asset and variance
  double L;     // Local volatility scaling factor
  double T;     // Time to maturity
  double K;     // Strike price

  HestonSLVParameters();
};

class HestonSLVCoefficients {
private:
  HestonSLVParameters &params;
  std::map<std::string, std::unique_ptr<Coefficient>> coeff_cache;

public:
  HestonSLVCoefficients(HestonSLVParameters &p);
  ~HestonSLVCoefficients();

  class SDriftCoefficient : public Coefficient {
  private:
    HestonSLVParameters &params;

  public:
    SDriftCoefficient(HestonSLVParameters &p);
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  };

  class VDriftCoefficient : public Coefficient {
  private:
    HestonSLVParameters &params;

  public:
    VDriftCoefficient(HestonSLVParameters &p);
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  };

  class SDiffusionCoefficient : public Coefficient {
  private:
    HestonSLVParameters &params;

  public:
    SDiffusionCoefficient(HestonSLVParameters &p);
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  };

  class VDiffusionCoefficient : public Coefficient {
  private:
    HestonSLVParameters &params;

  public:
    VDiffusionCoefficient(HestonSLVParameters &p);
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  };

  class MixedDerivativeCoefficient : public Coefficient {
  private:
    HestonSLVParameters &params;

  public:
    MixedDerivativeCoefficient(HestonSLVParameters &p);
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  };

  class InitialCondition : public Coefficient {
  private:
    HestonSLVParameters &params;

  public:
    InitialCondition(HestonSLVParameters &p);
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  };

  class FarFieldCondition : public Coefficient {
  private:
    HestonSLVParameters &params;
    double time;

  public:
    FarFieldCondition(HestonSLVParameters &p);
    void SetTime(double t);
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  };

  class ReactionCoefficient : public Coefficient {
  private:
    HestonSLVParameters &params;

  public:
    ReactionCoefficient(HestonSLVParameters &p);
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  };

  // Get or create coefficient instances with caching
  Coefficient *GetSDriftCoeff();
  Coefficient *GetVDriftCoeff();
  Coefficient *GetSDiffusionCoeff();
  Coefficient *GetVDiffusionCoeff();
  Coefficient *GetMixedDerivativeCoeff();
  Coefficient *GetInitialConditionCoeff();
  Coefficient *GetFarFieldConditionCoeff();
  Coefficient *GetReactionCoeff();

  // Set time for time-dependent coefficients
  void SetTime(double t);

  // Helper method to get all coefficients
  void GetCoefficients(std::unique_ptr<SDriftCoefficient> &s_drift,
                       std::unique_ptr<VDriftCoefficient> &v_drift,
                       std::unique_ptr<SDiffusionCoefficient> &s_diff,
                       std::unique_ptr<VDiffusionCoefficient> &v_diff,
                       std::unique_ptr<MixedDerivativeCoefficient> &mixed_diff,
                       std::unique_ptr<InitialCondition> &init_cond,
                       std::unique_ptr<FarFieldCondition> &far_field);
};

// Define a wrapper for IterativeSolver as Operator
class SolverOperator : public Operator {
private:
  IterativeSolver &solver_;
  const Operator *op_;

public:
  SolverOperator(IterativeSolver &solver)
      : Operator(solver.Height()), solver_(solver), op_(nullptr) {
    // We can't directly access the operator, but we know it exists
    // since the solver would fail without one
  }

  virtual void Mult(const Vector &x, Vector &y) const override {
    y = 0.0;
    solver_.Mult(x, y);
  }

  virtual void MultTranspose(const Vector &x, Vector &y) const override {
    mfem_error("SolverOperator::MultTranspose not implemented!");
  }
};

class UltraweakDPGHestonSLV {
public:
  enum Var {
    PRICE = 0,
    SIGMA_S = 1,
    SIGMA_V = 2,
    SIGMA_SS = 3,
    SIGMA_SV = 4,
    SIGMA_VV = 5,
    NUM_VARS = 6
  };

private:
  HestonSLVParameters params;

  double Smax, Vmax;
  int n_S, n_V;

  int trial_order;
  int test_order;
  int time_steps;
  double dt;

  // Parallel Mesh and FE Spaces
  ParMesh *mesh_ptr_;
  std::vector<FiniteElementCollection *> trial_fec_;
  std::vector<FiniteElementCollection *> test_fec_;
  std::vector<ParFiniteElementSpace *> trial_fes_;
  std::vector<ParFiniteElementSpace *> test_fes_;

  Array<int> trial_block_offsets_;
  Array<int> test_block_offsets_;

  // Solution Vectors
  BlockVector *U_;
  BlockVector *U_prev_;

  // Coefficient Management
  HestonSLVCoefficients *coeffs_;

  // Forms for system assembly
  std::vector<std::vector<std::unique_ptr<ParMixedBilinearForm>>> B_forms_;
  std::vector<std::unique_ptr<ParBilinearForm>> G_forms_;
  std::vector<std::unique_ptr<ParLinearForm>> F_forms_;

  // System operators and vectors
  std::unique_ptr<BlockOperator> B_;
  std::unique_ptr<BlockOperator> G_op_;
  std::unique_ptr<BlockVector> F_vec_;

  // DPG Solver components
  std::vector<std::unique_ptr<HypreBoomerAMG>> G_amg_precs_;
  std::unique_ptr<HyprePCG> G_solver_;
  std::unique_ptr<Operator> G_inv_op_;
  std::unique_ptr<Operator> B_T_op_;
  std::unique_ptr<Operator> A_;
  std::unique_ptr<HypreGMRES> A_solver_;

  // Adaptive refinement
  bool use_adaptive_refinement_;
  double error_fraction_;
  double threshold_factor_;
  int max_refinement_iterations_;

  // Visualization
  mfem::socketstream sout_;
  bool visualization_enabled_;

  // Solver settings
  double tolerance_;
  int max_iterations_;

  // Private methods
  void InitializeMesh();
  void InitializeSpaces();
  void InitializeSolution();

  void AssembleSystem(double current_time);
  void SetupSolvers();
  void SolveTimeStep(double current_time);
  void PerformAdaptiveRefinement();

  void SetupVisualization();
  void UpdateVisualization(double time);

public:
  UltraweakDPGHestonSLV(int trial_order_ = 1, int test_enrichment_ = 1);
  ~UltraweakDPGHestonSLV();

  // Configuration methods
  void SetParameters(const HestonSLVParameters &p);
  void SetDomainSize(double s_max, double v_max);
  void SetMeshResolution(int ns, int nv);
  void SetTimeSteps(int steps);
  void SetAdaptiveRefinement(bool adaptive, double frac = 0.7, double thr = 0.5,
                             int max_iter = 3);
  void EnableVisualization(bool viz = true);
  void SetLinearSolverOptions(double tol, int max_iter);

  // Main workflow methods
  void Initialize();
  void Solve();

  // Access methods
  const ParGridFunction *GetSolutionGridFunction(Var v = PRICE) const;
  const BlockVector *GetSolutionVector() const;
  const ParMesh *GetMesh() const;
  void SaveSolution(const std::string &filename_prefix);
  void SaveMesh(const std::string &filename);

  // Error estimation and adaptivity
  void ComputeErrorEstimates(Vector &errors);
  double ComputeTotalError();
};

} // namespace mfem

#endif // HESTON_SLV_DPG_HPP