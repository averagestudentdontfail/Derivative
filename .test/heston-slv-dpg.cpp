#include "heston-slv-dpg.hpp"
#include "dpg-error-estimator.hpp"
#include "dpg-integrators.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

using namespace std;
using namespace mfem;
using namespace std::chrono;

// HestonSLVParameters implementation
HestonSLVParameters::HestonSLVParameters()
    : r(0.05), q(0.0), kappa(2.0), theta(0.04), xi(0.3), rho(-0.7), L(1.0),
      T(1.0), K(100.0) {}

// HestonSLVCoefficients implementation
HestonSLVCoefficients::HestonSLVCoefficients(HestonSLVParameters &p)
    : params(p) {}

HestonSLVCoefficients::~HestonSLVCoefficients() { coeff_cache.clear(); }

// SDriftCoefficient implementation
HestonSLVCoefficients::SDriftCoefficient::SDriftCoefficient(
    HestonSLVParameters &p)
    : params(p) {}

double
HestonSLVCoefficients::SDriftCoefficient::Eval(ElementTransformation &T,
                                               const IntegrationPoint &ip) {
  Vector pos;
  T.Transform(ip, pos);
  double S = pos(0);
  return (params.r - params.q) * S;
}

// VDriftCoefficient implementation
HestonSLVCoefficients::VDriftCoefficient::VDriftCoefficient(
    HestonSLVParameters &p)
    : params(p) {}

double
HestonSLVCoefficients::VDriftCoefficient::Eval(ElementTransformation &T,
                                               const IntegrationPoint &ip) {
  Vector pos;
  T.Transform(ip, pos);
  double V = pos(1);
  return params.kappa * (params.theta - V);
}

// SDiffusionCoefficient implementation
HestonSLVCoefficients::SDiffusionCoefficient::SDiffusionCoefficient(
    HestonSLVParameters &p)
    : params(p) {}

double
HestonSLVCoefficients::SDiffusionCoefficient::Eval(ElementTransformation &T,
                                                   const IntegrationPoint &ip) {
  Vector pos;
  T.Transform(ip, pos);
  double S = pos(0);
  double V = pos(1);
  return 0.5 * params.L * params.L * std::max(V, 0.0) * S * S;
}

// VDiffusionCoefficient implementation
HestonSLVCoefficients::VDiffusionCoefficient::VDiffusionCoefficient(
    HestonSLVParameters &p)
    : params(p) {}

double
HestonSLVCoefficients::VDiffusionCoefficient::Eval(ElementTransformation &T,
                                                   const IntegrationPoint &ip) {
  Vector pos;
  T.Transform(ip, pos);
  double V = pos(1);
  return 0.5 * params.xi * params.xi * std::max(V, 0.0);
}

// MixedDerivativeCoefficient implementation
HestonSLVCoefficients::MixedDerivativeCoefficient::MixedDerivativeCoefficient(
    HestonSLVParameters &p)
    : params(p) {}

double HestonSLVCoefficients::MixedDerivativeCoefficient::Eval(
    ElementTransformation &T, const IntegrationPoint &ip) {
  Vector pos;
  T.Transform(ip, pos);
  double S = pos(0);
  double V = pos(1);
  return params.rho * params.xi * params.L * std::max(V, 0.0) * S;
}

// ReactionCoefficient implementation
HestonSLVCoefficients::ReactionCoefficient::ReactionCoefficient(
    HestonSLVParameters &p)
    : params(p) {}

double
HestonSLVCoefficients::ReactionCoefficient::Eval(ElementTransformation &T,
                                                 const IntegrationPoint &ip) {
  return -params.r; // Negative for the reaction term in the PDE
}

// InitialCondition implementation
HestonSLVCoefficients::InitialCondition::InitialCondition(
    HestonSLVParameters &p)
    : params(p) {}

double
HestonSLVCoefficients::InitialCondition::Eval(ElementTransformation &T,
                                              const IntegrationPoint &ip) {
  Vector pos;
  T.Transform(ip, pos);
  double S = pos(0);
  // European call payoff: max(S-K, 0)
  return std::max(S - params.K, 0.0);
}

// FarFieldCondition implementation
HestonSLVCoefficients::FarFieldCondition::FarFieldCondition(
    HestonSLVParameters &p)
    : params(p), time(0.0) {}

void HestonSLVCoefficients::FarFieldCondition::SetTime(double t) { time = t; }

double
HestonSLVCoefficients::FarFieldCondition::Eval(ElementTransformation &T,
                                               const IntegrationPoint &ip) {
  Vector pos;
  T.Transform(ip, pos);
  double S = pos(0);
  double tau = params.T - time;
  return S - params.K * exp(-params.r * std::max(tau, 0.0));
}

// Coefficient cache access methods
Coefficient *HestonSLVCoefficients::GetSDriftCoeff() {
  auto it = coeff_cache.find("s_drift");
  if (it == coeff_cache.end()) {
    coeff_cache["s_drift"] = std::make_unique<SDriftCoefficient>(params);
  }
  return coeff_cache["s_drift"].get();
}

Coefficient *HestonSLVCoefficients::GetVDriftCoeff() {
  auto it = coeff_cache.find("v_drift");
  if (it == coeff_cache.end()) {
    coeff_cache["v_drift"] = std::make_unique<VDriftCoefficient>(params);
  }
  return coeff_cache["v_drift"].get();
}

Coefficient *HestonSLVCoefficients::GetSDiffusionCoeff() {
  auto it = coeff_cache.find("s_diff");
  if (it == coeff_cache.end()) {
    coeff_cache["s_diff"] = std::make_unique<SDiffusionCoefficient>(params);
  }
  return coeff_cache["s_diff"].get();
}

Coefficient *HestonSLVCoefficients::GetVDiffusionCoeff() {
  auto it = coeff_cache.find("v_diff");
  if (it == coeff_cache.end()) {
    coeff_cache["v_diff"] = std::make_unique<VDiffusionCoefficient>(params);
  }
  return coeff_cache["v_diff"].get();
}

Coefficient *HestonSLVCoefficients::GetMixedDerivativeCoeff() {
  auto it = coeff_cache.find("mixed_diff");
  if (it == coeff_cache.end()) {
    coeff_cache["mixed_diff"] =
        std::make_unique<MixedDerivativeCoefficient>(params);
  }
  return coeff_cache["mixed_diff"].get();
}

Coefficient *HestonSLVCoefficients::GetInitialConditionCoeff() {
  auto it = coeff_cache.find("init_cond");
  if (it == coeff_cache.end()) {
    coeff_cache["init_cond"] = std::make_unique<InitialCondition>(params);
  }
  return coeff_cache["init_cond"].get();
}

Coefficient *HestonSLVCoefficients::GetFarFieldConditionCoeff() {
  auto it = coeff_cache.find("far_field");
  if (it == coeff_cache.end()) {
    coeff_cache["far_field"] = std::make_unique<FarFieldCondition>(params);
  }
  return coeff_cache["far_field"].get();
}

Coefficient *HestonSLVCoefficients::GetReactionCoeff() {
  auto it = coeff_cache.find("reaction");
  if (it == coeff_cache.end()) {
    coeff_cache["reaction"] = std::make_unique<ReactionCoefficient>(params);
  }
  return coeff_cache["reaction"].get();
}

void HestonSLVCoefficients::SetTime(double t) {
  auto it = coeff_cache.find("far_field");
  if (it != coeff_cache.end()) {
    dynamic_cast<FarFieldCondition *>(it->second.get())->SetTime(t);
  }
}

void HestonSLVCoefficients::GetCoefficients(
    std::unique_ptr<SDriftCoefficient> &s_drift,
    std::unique_ptr<VDriftCoefficient> &v_drift,
    std::unique_ptr<SDiffusionCoefficient> &s_diff,
    std::unique_ptr<VDiffusionCoefficient> &v_diff,
    std::unique_ptr<MixedDerivativeCoefficient> &mixed_diff,
    std::unique_ptr<InitialCondition> &init_cond,
    std::unique_ptr<FarFieldCondition> &far_field) {
  s_drift = std::make_unique<SDriftCoefficient>(params);
  v_drift = std::make_unique<VDriftCoefficient>(params);
  s_diff = std::make_unique<SDiffusionCoefficient>(params);
  v_diff = std::make_unique<VDiffusionCoefficient>(params);
  mixed_diff = std::make_unique<MixedDerivativeCoefficient>(params);
  init_cond = std::make_unique<InitialCondition>(params);
  far_field = std::make_unique<FarFieldCondition>(params);
}

// UltraweakDPGHestonSLV implementation
UltraweakDPGHestonSLV::UltraweakDPGHestonSLV(int trial_order_,
                                             int test_enrichment_)
    : Smax(300.0), Vmax(1.0), n_S(50), n_V(30), trial_order(trial_order_),
      test_order(trial_order_ + test_enrichment_), time_steps(50), dt(0.0),
      mesh_ptr_(nullptr), trial_fec_(Var::NUM_VARS, nullptr),
      test_fec_(Var::NUM_VARS, nullptr), trial_fes_(Var::NUM_VARS, nullptr),
      test_fes_(Var::NUM_VARS, nullptr), U_(nullptr), U_prev_(nullptr),
      coeffs_(nullptr), use_adaptive_refinement_(false), error_fraction_(0.7),
      threshold_factor_(0.5), max_refinement_iterations_(3),
      visualization_enabled_(false), tolerance_(1e-6), max_iterations_(100) {
  if (test_order <= trial_order) {
    test_order = trial_order + 1;
    cout << "Warning: Test order enrichment must be at least 1. Setting "
            "test_order to "
         << test_order << endl;
  }

  params = HestonSLVParameters();
  coeffs_ = new HestonSLVCoefficients(params);
  dt = params.T / time_steps;
}

UltraweakDPGHestonSLV::~UltraweakDPGHestonSLV() {
  delete coeffs_;
  delete U_;
  delete U_prev_;

  for (auto fes : trial_fes_) {
    delete fes;
  }
  for (auto fes : test_fes_) {
    delete fes;
  }

  if (trial_fec_[0])
    delete trial_fec_[0];
  if (test_fec_[0])
    delete test_fec_[0];

  delete mesh_ptr_;
}

void UltraweakDPGHestonSLV::SetParameters(const HestonSLVParameters &p) {
  params = p;
  delete coeffs_;
  coeffs_ = new HestonSLVCoefficients(params);
  dt = params.T / time_steps;
}

void UltraweakDPGHestonSLV::SetDomainSize(double s_max, double v_max) {
  Smax = s_max;
  Vmax = v_max;
}

void UltraweakDPGHestonSLV::SetMeshResolution(int ns, int nv) {
  n_S = ns;
  n_V = nv;
}

void UltraweakDPGHestonSLV::SetTimeSteps(int steps) {
  if (steps <= 0) {
    cout << "Warning: Number of time steps must be positive. Using default: "
         << time_steps << endl;
    return;
  }
  time_steps = steps;
  dt = params.T / time_steps;
}

void UltraweakDPGHestonSLV::SetAdaptiveRefinement(bool adaptive, double frac,
                                                  double thr, int max_iter) {
  use_adaptive_refinement_ = adaptive;
  error_fraction_ = frac;
  threshold_factor_ = thr;
  max_refinement_iterations_ = max_iter;
}

void UltraweakDPGHestonSLV::EnableVisualization(bool viz) {
  visualization_enabled_ = viz;
}

void UltraweakDPGHestonSLV::SetLinearSolverOptions(double tol, int max_iter) {
  tolerance_ = tol;
  max_iterations_ = max_iter;

  // Apply to existing solvers if they exist
  if (A_solver_) {
    A_solver_->SetTol(tolerance_);
    A_solver_->SetMaxIter(max_iterations_);
  }
  if (G_solver_) {
    G_solver_->SetTol(tolerance_ * 0.1);
    G_solver_->SetMaxIter(max_iterations_);
  }
}

void UltraweakDPGHestonSLV::Initialize() {
  cout << "Initializing Heston SLV DPG Solver..." << endl;

  // Initialize MPI if not already done
  int mpi_initialized;
  MPI_Initialized(&mpi_initialized);
  if (!mpi_initialized) {
    int argc = 0;
    char **argv = nullptr;
    MPI_Init(&argc, &argv);
  }

  InitializeMesh();
  InitializeSpaces();
  InitializeSolution();

  if (visualization_enabled_) {
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid == 0) {
      SetupVisualization();
    }
  }

  cout << "Initialization complete." << endl;
}

void UltraweakDPGHestonSLV::InitializeMesh() {
  cout << "Creating mesh with " << n_S << "x" << n_V << " elements..." << endl;

  delete mesh_ptr_;

  // First create a serial mesh
  Mesh *serial_mesh = new Mesh(Mesh::MakeCartesian2D(
      n_S, n_V, Element::QUADRILATERAL, true, Smax, Vmax));

  // Create a parallel mesh from the serial mesh
  mesh_ptr_ = new ParMesh(MPI_COMM_WORLD, *serial_mesh);
  delete serial_mesh;

  // Set up element attributes for boundaries
  // Boundary 1: S=0  (left)
  // Boundary 2: S=Smax (right)
  // Boundary 3: V=0 (bottom)
  // Boundary 4: V=Vmax (top)

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (myid == 0) {
    cout << "Parallel mesh created with global elements: " << mesh_ptr_->GetNE()
         << endl;
  }
}

void UltraweakDPGHestonSLV::InitializeSpaces() {
  cout << "Creating Finite Element Spaces..." << endl;

  MFEM_VERIFY(mesh_ptr_, "Mesh must be initialized before spaces");

  // Clean up existing spaces and collections
  for (auto fes : trial_fes_) {
    delete fes;
  }
  for (auto fes : test_fes_) {
    delete fes;
  }

  if (trial_fec_[0])
    delete trial_fec_[0];
  if (test_fec_[0])
    delete test_fec_[0];

  trial_fes_.assign(Var::NUM_VARS, nullptr);
  test_fes_.assign(Var::NUM_VARS, nullptr);
  trial_fec_.assign(Var::NUM_VARS, nullptr);
  test_fec_.assign(Var::NUM_VARS, nullptr);

  // Create FE Collections
  auto trial_coll = new DG_FECollection(trial_order, mesh_ptr_->Dimension());
  auto test_coll = new DG_FECollection(test_order, mesh_ptr_->Dimension());

  trial_fec_[0] = trial_coll;
  test_fec_[0] = test_coll;

  // Set up block offsets
  trial_block_offsets_.SetSize(Var::NUM_VARS + 1);
  test_block_offsets_.SetSize(Var::NUM_VARS + 1);
  trial_block_offsets_[0] = 0;
  test_block_offsets_[0] = 0;

  // Create FE spaces for each variable
  for (int i = 0; i < Var::NUM_VARS; ++i) {
    // Share FE collections for efficiency
    if (i > 0) {
      trial_fec_[i] = trial_fec_[0];
      test_fec_[i] = test_fec_[0];
    }

    trial_fes_[i] = new ParFiniteElementSpace(mesh_ptr_, trial_fec_[i]);
    test_fes_[i] = new ParFiniteElementSpace(mesh_ptr_, test_fec_[i]);

    // Store local size before converting to global offsets
    trial_block_offsets_[i + 1] = trial_fes_[i]->GetTrueVSize();
    test_block_offsets_[i + 1] = test_fes_[i]->GetTrueVSize();
  }

  // Convert sizes to offsets
  trial_block_offsets_.PartialSum();
  test_block_offsets_.PartialSum();

  long long trial_vdofs = trial_block_offsets_.Last();
  long long test_vdofs = test_block_offsets_.Last();

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (myid == 0) {
    cout << "  Trial Order: " << trial_order << ", Test Order: " << test_order
         << endl;
    cout << "  Total Global True Trial DOFs: " << trial_vdofs << endl;
    cout << "  Total Global True Test DOFs: " << test_vdofs << endl;
  }

  // Create BlockVectors for solution
  delete U_;
  delete U_prev_;
  U_ = new BlockVector(trial_block_offsets_);
  U_prev_ = new BlockVector(trial_block_offsets_);
  *U_ = 0.0;
  *U_prev_ = 0.0;
}

void UltraweakDPGHestonSLV::InitializeSolution() {
  cout << "Initializing Solution (Terminal Condition)..." << endl;

  MFEM_VERIFY(U_ && trial_fes_[Var::PRICE],
              "Spaces and solution vectors must be initialized first");

  // Project the terminal payoff condition
  ParGridFunction u_gf(trial_fes_[Var::PRICE]);

  // Get the initial condition coefficient
  Coefficient *init_coeff = coeffs_->GetInitialConditionCoeff();
  MFEM_VERIFY(init_coeff, "Failed to get InitialCondition coefficient");

  // Project the initial condition
  u_gf.ProjectCoefficient(*init_coeff);

  // Copy the projected data into the solution vector
  u_gf.GetTrueDofs(U_->GetBlock(Var::PRICE));

  // Initialize derivative variables to zero
  for (int i = Var::SIGMA_S; i <= Var::SIGMA_VV; ++i) {
    U_->GetBlock(i) = 0.0;
  }

  // Set previous solution
  *U_prev_ = *U_;

  cout << "Solution initialized with option payoff." << endl;
}

void UltraweakDPGHestonSLV::SetupVisualization() {
  if (!visualization_enabled_ || !mesh_ptr_)
    return;

  cout << "Setting up GLVis visualization..." << endl;
  char vishost[] = "localhost";
  int visport = 19916;

  sout_.open(vishost, visport);
  if (!sout_) {
    cout << "Warning: Unable to connect to GLVis server at " << vishost << ':'
         << visport << endl;
    cout << "Disabling visualization." << endl;
    visualization_enabled_ = false;
    return;
  }

  sout_.precision(8);

  // Send initial mesh and solution
  UpdateVisualization(params.T);
}

void UltraweakDPGHestonSLV::UpdateVisualization(double time) {
  if (!visualization_enabled_ || !mesh_ptr_ || !U_ || !trial_fes_[Var::PRICE])
    return;

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (myid != 0)
    return;

  if (!sout_)
    return;

  // Create a view of the price solution
  ParGridFunction u_gf_view(trial_fes_[Var::PRICE]);
  u_gf_view.Distribute(U_->GetBlock(Var::PRICE));

  // Send data to GLVis
  sout_ << "parallel " << mesh_ptr_->GetNRanks() << " "
        << mesh_ptr_->GetMyRank() << "\n";
  sout_ << "solution\n"
        << *mesh_ptr_ << u_gf_view
        << "window_title 'Option Price at t = " << time << "'" << std::flush;

  // Optionally add visualizations for derivatives (Delta, Gamma, etc.)
}

void UltraweakDPGHestonSLV::AssembleSystem(double current_time) {
  cout << "Assembling DPG System for time " << current_time << "..." << endl;

  MFEM_VERIFY(trial_fes_.size() == Var::NUM_VARS &&
                  test_fes_.size() == Var::NUM_VARS,
              "FE spaces not fully initialized");

  // Clear previous forms
  B_forms_.assign(
      Var::NUM_VARS,
      std::vector<std::unique_ptr<ParMixedBilinearForm>>(Var::NUM_VARS));
  G_forms_.assign(Var::NUM_VARS, nullptr);
  F_forms_.assign(Var::NUM_VARS, nullptr);

  // Define coefficients
  ConstantCoefficient one(1.0);
  ConstantCoefficient neg_one(-1.0);
  ConstantCoefficient zero(0.0);
  ConstantCoefficient dt_inv_coeff(1.0 / dt);

  // Get model coefficients
  coeffs_->SetTime(current_time);

  Coefficient *s_drift = coeffs_->GetSDriftCoeff();
  Coefficient *v_drift = coeffs_->GetVDriftCoeff();
  Coefficient *s_diff = coeffs_->GetSDiffusionCoeff();
  Coefficient *v_diff = coeffs_->GetVDiffusionCoeff();
  Coefficient *mixed_diff = coeffs_->GetMixedDerivativeCoeff();
  Coefficient *reaction = coeffs_->GetReactionCoeff();
  Coefficient *far_field = coeffs_->GetFarFieldConditionCoeff();

  // Time-dependent coefficient for the previous time step solution
  ParGridFunction u_prev_gf(trial_fes_[Var::PRICE]);
  u_prev_gf.Distribute(U_prev_->GetBlock(Var::PRICE));
  GridFunctionCoefficient u_prev_coeff(&u_prev_gf);
  ProductCoefficient time_term_coeff(dt_inv_coeff, u_prev_coeff);

  // Define boundary attributes
  int max_bdr_attr = mesh_ptr_->bdr_attributes.Max();
  Array<int> s0_bdr_attr, smax_bdr_attr, v0_bdr_attr, vmax_bdr_attr;

  if (max_bdr_attr >= 1)
    s0_bdr_attr.Append(1);
  if (max_bdr_attr >= 2)
    smax_bdr_attr.Append(2);
  if (max_bdr_attr >= 3)
    v0_bdr_attr.Append(3);
  if (max_bdr_attr >= 4)
    vmax_bdr_attr.Append(4);

  // Add integrators for B forms (trial × test)
  cout << "  Adding integrators for B forms..." << endl;

  // --- 1. Equation: sigma_S - du/dS = 0 ---
  B_forms_[Var::SIGMA_S][Var::SIGMA_S] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_S], test_fes_[Var::SIGMA_S]);
  B_forms_[Var::SIGMA_S][Var::SIGMA_S]->AddDomainIntegrator(
      new MassIntegrator(one));

  B_forms_[Var::SIGMA_S][Var::PRICE] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::PRICE], test_fes_[Var::SIGMA_S]);
  B_forms_[Var::SIGMA_S][Var::PRICE]->AddDomainIntegrator(
      new UltraweakDerivativeIntegrator(0));
  B_forms_[Var::SIGMA_S][Var::PRICE]->AddTraceFaceIntegrator(
      new UltraweakJumpIntegrator(0));

  // --- 2. Equation: sigma_V - du/dV = 0 ---
  B_forms_[Var::SIGMA_V][Var::SIGMA_V] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_V], test_fes_[Var::SIGMA_V]);
  B_forms_[Var::SIGMA_V][Var::SIGMA_V]->AddDomainIntegrator(
      new MassIntegrator(one));

  B_forms_[Var::SIGMA_V][Var::PRICE] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::PRICE], test_fes_[Var::SIGMA_V]);
  B_forms_[Var::SIGMA_V][Var::PRICE]->AddDomainIntegrator(
      new UltraweakDerivativeIntegrator(1));
  B_forms_[Var::SIGMA_V][Var::PRICE]->AddTraceFaceIntegrator(
      new UltraweakJumpIntegrator(1));

  // Add boundary integrators for V boundaries
  if (vmax_bdr_attr.Size() > 0) {
    B_forms_[Var::SIGMA_V][Var::PRICE]->AddBoundaryIntegrator(
        new BoundaryMassIntegrator(neg_one), vmax_bdr_attr);
  }
  if (v0_bdr_attr.Size() > 0) {
    B_forms_[Var::SIGMA_V][Var::PRICE]->AddBoundaryIntegrator(
        new BoundaryMassIntegrator(one), v0_bdr_attr);
  }

  // --- 3. Equation: sigma_SS - d(sigma_S)/dS = 0 ---
  B_forms_[Var::SIGMA_SS][Var::SIGMA_SS] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_SS], test_fes_[Var::SIGMA_SS]);
  B_forms_[Var::SIGMA_SS][Var::SIGMA_SS]->AddDomainIntegrator(
      new MassIntegrator(one));

  B_forms_[Var::SIGMA_SS][Var::SIGMA_S] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_S], test_fes_[Var::SIGMA_SS]);
  B_forms_[Var::SIGMA_SS][Var::SIGMA_S]->AddDomainIntegrator(
      new UltraweakDerivativeIntegrator(0));
  B_forms_[Var::SIGMA_SS][Var::SIGMA_S]->AddTraceFaceIntegrator(
      new UltraweakJumpIntegrator(0));

  // --- 4. Equation: sigma_SV - d(sigma_S)/dV = 0 ---
  B_forms_[Var::SIGMA_SV][Var::SIGMA_SV] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_SV], test_fes_[Var::SIGMA_SV]);
  B_forms_[Var::SIGMA_SV][Var::SIGMA_SV]->AddDomainIntegrator(
      new MassIntegrator(one));

  B_forms_[Var::SIGMA_SV][Var::SIGMA_S] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_S], test_fes_[Var::SIGMA_SV]);
  B_forms_[Var::SIGMA_SV][Var::SIGMA_S]->AddDomainIntegrator(
      new UltraweakDerivativeIntegrator(1));
  B_forms_[Var::SIGMA_SV][Var::SIGMA_S]->AddTraceFaceIntegrator(
      new UltraweakJumpIntegrator(1));

  // Add boundary integrators for V boundaries
  if (vmax_bdr_attr.Size() > 0) {
    B_forms_[Var::SIGMA_SV][Var::SIGMA_S]->AddBoundaryIntegrator(
        new BoundaryMassIntegrator(neg_one), vmax_bdr_attr);
  }
  if (v0_bdr_attr.Size() > 0) {
    B_forms_[Var::SIGMA_SV][Var::SIGMA_S]->AddBoundaryIntegrator(
        new BoundaryMassIntegrator(one), v0_bdr_attr);
  }

  // --- 5. Equation: sigma_VV - d(sigma_V)/dV = 0 ---
  B_forms_[Var::SIGMA_VV][Var::SIGMA_VV] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_VV], test_fes_[Var::SIGMA_VV]);
  B_forms_[Var::SIGMA_VV][Var::SIGMA_VV]->AddDomainIntegrator(
      new MassIntegrator(one));

  B_forms_[Var::SIGMA_VV][Var::SIGMA_V] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_V], test_fes_[Var::SIGMA_VV]);
  B_forms_[Var::SIGMA_VV][Var::SIGMA_V]->AddDomainIntegrator(
      new UltraweakDerivativeIntegrator(1));
  B_forms_[Var::SIGMA_VV][Var::SIGMA_V]->AddTraceFaceIntegrator(
      new UltraweakJumpIntegrator(1));

  // Add boundary integrators for V boundaries
  if (v0_bdr_attr.Size() > 0) {
    B_forms_[Var::SIGMA_VV][Var::SIGMA_V]->AddBoundaryIntegrator(
        new BoundaryMassIntegrator(one), v0_bdr_attr);
  }

  // --- 6. Equation: Primary PDE (Backward Euler) ---

  // Time derivative and reaction terms
  B_forms_[Var::PRICE][Var::PRICE] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::PRICE], test_fes_[Var::PRICE]);
  B_forms_[Var::PRICE][Var::PRICE]->AddDomainIntegrator(
      new MassIntegrator(dt_inv_coeff));
  B_forms_[Var::PRICE][Var::PRICE]->AddDomainIntegrator(
      new UltraweakReactionIntegrator(reaction));

  // Drift terms
  B_forms_[Var::PRICE][Var::SIGMA_S] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_S], test_fes_[Var::PRICE]);
  B_forms_[Var::PRICE][Var::SIGMA_S]->AddDomainIntegrator(
      new UltraweakDriftIntegrator(s_drift));

  B_forms_[Var::PRICE][Var::SIGMA_V] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_V], test_fes_[Var::PRICE]);
  B_forms_[Var::PRICE][Var::SIGMA_V]->AddDomainIntegrator(
      new UltraweakDriftIntegrator(v_drift));

  // Diffusion terms
  B_forms_[Var::PRICE][Var::SIGMA_SS] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_SS], test_fes_[Var::PRICE]);
  B_forms_[Var::PRICE][Var::SIGMA_SS]->AddDomainIntegrator(
      new UltraweakDriftIntegrator(s_diff));

  B_forms_[Var::PRICE][Var::SIGMA_VV] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_VV], test_fes_[Var::PRICE]);
  B_forms_[Var::PRICE][Var::SIGMA_VV]->AddDomainIntegrator(
      new UltraweakDriftIntegrator(v_diff));

  // Mixed derivative term
  B_forms_[Var::PRICE][Var::SIGMA_SV] = make_unique<ParMixedBilinearForm>(
      trial_fes_[Var::SIGMA_SV], test_fes_[Var::PRICE]);
  B_forms_[Var::PRICE][Var::SIGMA_SV]->AddDomainIntegrator(
      new UltraweakDriftIntegrator(mixed_diff));

  // Add integrators for G forms (test × test)
  cout << "  Adding integrators for G forms..." << endl;

  for (int i = 0; i < Var::NUM_VARS; ++i) {
    G_forms_[i] = make_unique<ParBilinearForm>(test_fes_[i]);
    G_forms_[i]->AddDomainIntegrator(new TestInnerProductIntegrator(1.0, 1.0));
  }

  // Add integrators for F forms (RHS)
  cout << "  Adding integrators for F forms..." << endl;

  for (int i = 0; i < Var::NUM_VARS; ++i) {
    F_forms_[i] = make_unique<ParLinearForm>(test_fes_[i]);
    *F_forms_[i] = 0.0;
  }

  // RHS for boundary conditions
  if (smax_bdr_attr.Size() > 0) {
    // Far-field boundary condition at S=Smax
    ProductCoefficient neg_far_field(-1.0, *far_field);
    F_forms_[Var::SIGMA_S]->AddBoundaryIntegrator(
        new BoundaryLFIntegrator(neg_far_field), smax_bdr_attr);

    // Delta (sigma_S) at S=Smax is 1
    ConstantCoefficient sigma_S_smax(1.0);
    ProductCoefficient neg_sigma_S_smax(-1.0, sigma_S_smax);
    F_forms_[Var::SIGMA_SS]->AddBoundaryIntegrator(
        new BoundaryLFIntegrator(neg_sigma_S_smax), smax_bdr_attr);
  }

  if (vmax_bdr_attr.Size() > 0) {
    // Sigma_VV at V=Vmax is 0 (Neumann condition)
    F_forms_[Var::SIGMA_VV]->AddBoundaryIntegrator(
        new BoundaryLFIntegrator(zero), vmax_bdr_attr);
  }

  // Time term for the option price equation
  F_forms_[Var::PRICE]->AddDomainIntegrator(
      new DomainLFIntegrator(time_term_coeff));

  // Assemble all forms
  cout << "  Assembling individual forms..." << endl;

  for (int i = 0; i < Var::NUM_VARS; ++i) {
    for (int j = 0; j < Var::NUM_VARS; ++j) {
      if (B_forms_[i][j]) {
        B_forms_[i][j]->Assemble();
        B_forms_[i][j]->Finalize();
      }
    }

    if (G_forms_[i]) {
      G_forms_[i]->Assemble();
      G_forms_[i]->Finalize();
    }

    if (F_forms_[i]) {
      F_forms_[i]->Assemble();
    }
  }

  // Create block operators and vectors
  cout << "  Creating block operators and RHS vector..." << endl;

  B_.reset(new BlockOperator(test_block_offsets_, trial_block_offsets_));
  G_op_.reset(new BlockOperator(test_block_offsets_, test_block_offsets_));
  F_vec_.reset(new BlockVector(test_block_offsets_));

  B_->owns_blocks = true;
  G_op_->owns_blocks = true;
  // No owns_blocks for BlockVector

  // Extract and set block matrices
  for (int i = 0; i < Var::NUM_VARS; ++i) {
    for (int j = 0; j < Var::NUM_VARS; ++j) {
      if (B_forms_[i][j]) {
        Array<int> trial_tdofs, test_tdofs;
        trial_fes_[j]->GetEssentialTrueDofs(Array<int>(), trial_tdofs);
        test_fes_[i]->GetEssentialTrueDofs(Array<int>(), test_tdofs);

        HypreParMatrix *mat = B_forms_[i][j]->ParallelAssemble();
        B_->SetBlock(i, j, mat);
      }
    }

    if (G_forms_[i]) {
      HypreParMatrix *mat = G_forms_[i]->ParallelAssemble();
      G_op_->SetBlock(i, i, mat);
    }

    if (F_forms_[i]) {
      HypreParVector *vec = F_forms_[i]->ParallelAssemble();
      // Correctly assign the vector to the block
      Vector &block = F_vec_->GetBlock(i);
      block = *vec;
      delete vec; // Free the newly allocated vector
    }
  }

  cout << "System assembly complete." << endl;
}

void UltraweakDPGHestonSLV::SetupSolvers() {
  cout << "Setting up DPG solvers..." << endl;

  // Check if operators are assembled
  if (!B_ || !G_op_) {
    throw std::runtime_error(
        "Operators B and G_op must be assembled before setting up solvers");
  }

  // Clear previous preconditioners
  G_amg_precs_.clear();

  // Create G Solver
  cout << "  Setting up G Solver (HyprePCG)..." << endl;
  G_solver_.reset(new HyprePCG(mesh_ptr_->GetComm()));
  G_solver_->SetTol(tolerance_ * 0.1); // Higher precision for inner solver
  G_solver_->SetMaxIter(max_iterations_);
  G_solver_->SetPrintLevel(0);
  G_solver_->SetOperator(*G_op_);

  // Set up AMG preconditioners for each diagonal block
  for (int i = 0; i < Var::NUM_VARS; ++i) {
    const HypreParMatrix *G_block =
        dynamic_cast<const HypreParMatrix *>(&G_op_->GetBlock(i, i));

    if (!G_block) {
      cout << "  Warning: Block " << i
           << " of G_op is NULL or not a HypreParMatrix" << endl;
      continue;
    }

    // Create AMG preconditioner for this diagonal block
    G_amg_precs_.push_back(std::make_unique<HypreBoomerAMG>(*G_block));
    auto &amg = G_amg_precs_.back();
    amg->SetPrintLevel(0);
    amg->SetMaxIter(1);
    amg->SetTol(0.0);

    // Use standard AMG parameters
    amg->SetStrengthThresh(0.25);
    amg->SetInterpType(0);
    amg->SetCoarsenType(6);
    amg->SetRelaxType(6); // Symmetric Gauss-Seidel
    amg->SetNumSweeps(1);
  }

  // Create a custom block-diagonal preconditioner that implements HypreSolver
  // interface
  class BlockDiagonalHypreSolver : public HypreSolver {
  private:
    const BlockOperator &A_;
    std::vector<HypreSolver *> block_solvers_;
    mutable BlockVector x_block_, y_block_;

  public:
    BlockDiagonalHypreSolver(const BlockOperator &A,
                             const std::vector<HypreSolver *> &block_solvers)
        : A_(A), block_solvers_(block_solvers), x_block_(A.RowOffsets()),
          y_block_(A.RowOffsets()) {
      MFEM_ASSERT(A.NumRowBlocks() == block_solvers.size(),
                  "Number of blocks doesn't match number of solvers");
      MFEM_ASSERT(A.NumRowBlocks() == A.NumColBlocks(),
                  "Only square block operators supported");
    }

    virtual void Mult(const HypreParVector &b,
                      HypreParVector &x) const override {
      // Convert to BlockVector
      x_block_.GetBlockView(0, x.Size(), x);
      y_block_.GetBlockView(0, b.Size(), const_cast<HypreParVector &>(b));

      // Apply block-diagonal preconditioner
      for (int i = 0; i < A_.NumRowBlocks(); i++) {
        const Operator *diag = &A_.GetBlock(i, i);
        if (diag && block_solvers_[i]) {
          // Get the block views
          int offset = A_.RowOffsets()[i];
          int size = A_.RowOffsets()[i + 1] - offset;

          // Make HypreParVector views for this block
          Vector b_block_view, x_block_view;
          y_block_.GetBlockView(i, b_block_view);
          x_block_.GetBlockView(i, x_block_view);

          // Cast to HypreParVector (assuming blocks are HypreParMatrix)
          HypreParVector *b_par = dynamic_cast<HypreParVector *>(&b_block_view);
          HypreParVector *x_par = dynamic_cast<HypreParVector *>(&x_block_view);

          if (b_par && x_par) {
            // Apply this block's solver
            block_solvers_[i]->Mult(*b_par, *x_par);
          }
        }
      }

      // Copy result back to x
      x = x_block_.GetBlock(0);
    }

    virtual void SetupOperator(const HypreParMatrix &op) override {
      // Nothing to do
    }
  };

  // Create vector of HypreSolver pointers for the custom preconditioner
  std::vector<HypreSolver *> block_solvers(Var::NUM_VARS, nullptr);
  for (int i = 0; i < Var::NUM_VARS; ++i) {
    if (i < G_amg_precs_.size()) {
      block_solvers[i] = G_amg_precs_[i].get();
    }
  }

  // Create the custom block-diagonal preconditioner
  auto block_prec = new BlockDiagonalHypreSolver(*G_op_, block_solvers);

  // Set the preconditioner for the G solver
  G_solver_->SetPreconditioner(*block_prec);

  // Create G^{-1} Operator wrapper
  cout << "  Creating G^{-1} Operator wrapper..." << endl;
  G_inv_op_.reset(new SolverOperator(*G_solver_));

  // Create B^T Operator
  cout << "  Creating B^T Operator..." << endl;
  B_T_op_.reset(new TransposeOperator(B_.get()));

  // Create A = B^T G^{-1} B Operator in two steps
  cout << "  Creating A = B^T G^{-1} B Operator..." << endl;

  // First create G^{-1} B
  Operator *GB_op =
      new ProductOperator(G_inv_op_.get(), B_.get(), false, false);

  // Then create B^T (G^{-1} B)
  A_.reset(new ProductOperator(B_T_op_.get(), GB_op, false,
                               true)); // true to own GB_op

  // Create A Solver
  cout << "  Setting up A Solver (HypreGMRES)..." << endl;
  A_solver_.reset(new HypreGMRES(mesh_ptr_->GetComm()));
  A_solver_->SetTol(tolerance_);
  A_solver_->SetMaxIter(max_iterations_);
  A_solver_->SetPrintLevel(1);
  A_solver_->SetOperator(*A_);

  cout << "DPG solvers set up successfully." << endl;
}

void UltraweakDPGHestonSLV::SolveTimeStep(double current_time) {
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == 0) {
    cout << "Solving time step for t = " << current_time << endl;
  }

  // Assemble the system for the current time step
  AssembleSystem(current_time);

  // Setup solvers based on the assembled system
  SetupSolvers();

  // Check if solvers and operators are ready
  if (!A_solver_ || !B_T_op_ || !G_inv_op_ || !F_vec_) {
    throw std::runtime_error(
        "Solvers and operators not ready for time step solve");
  }

  // Calculate the effective RHS: RHS_eff = B^T G^{-1} F
  BlockVector G_inv_F(test_block_offsets_);
  G_inv_F = 0.0;

  BlockVector RHS_eff(trial_block_offsets_);
  RHS_eff = 0.0;

  if (myid == 0) {
    cout << "  Computing RHS = B^T G^{-1} F..." << endl;
  }

  // First apply G^{-1} to F
  G_inv_op_->Mult(*F_vec_, G_inv_F);

  // Then apply B^T to G^{-1}F
  B_T_op_->Mult(G_inv_F, RHS_eff);

  // Solve the final system: A * U = RHS_eff
  if (myid == 0) {
    cout << "  Solving A * U = RHS..." << endl;
  }

  // Save the previous solution
  *U_prev_ = *U_;

  // Solve for the new solution
  A_solver_->Mult(RHS_eff, *U_);

  if (myid == 0) {
    cout << "  Solver iterations: " << A_solver_->GetNumIterations() << endl;
    cout << "  Final relative residual: " << A_solver_->GetFinalResidual()
         << endl;
    cout << "Time step solved." << endl;
  }

  // Update visualization if enabled
  if (visualization_enabled_) {
    if (myid == 0) {
      UpdateVisualization(current_time);
    }
  }
}

void UltraweakDPGHestonSLV::Solve() {
  cout << "\nStarting Heston-SLV option pricing (Ultraweak DPG)..." << endl;

  if (!U_ || !U_prev_) {
    throw std::runtime_error(
        "Solver must be initialized before calling Solve()");
  }

  double current_time = params.T;

  // Timing
  high_resolution_clock::time_point t_start = high_resolution_clock::now();

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Time stepping loop (backward in time)
  for (int step = 0; step < time_steps; step++) {
    current_time -= dt;

    if (myid == 0) {
      cout << "\n--- Time step " << step + 1 << " / " << time_steps
           << ", time = " << fixed << setprecision(4) << current_time << " ---"
           << endl;
    }

    // Perform one time step solve
    try {
      SolveTimeStep(current_time);

      // Perform adaptive refinement if enabled
      if (use_adaptive_refinement_ &&
          (step + 1) % max_refinement_iterations_ == 0) {
        PerformAdaptiveRefinement();
      }
    } catch (const std::exception &e) {
      cerr << "Error during time step " << step + 1 << ": " << e.what() << endl;
      throw;
    }
  }

  auto total_time =
      duration_cast<milliseconds>(high_resolution_clock::now() - t_start)
          .count();

  if (myid == 0) {
    cout << "\n-------------------------------------------------" << endl;
    cout << "Option pricing completed in " << total_time / 1000.0 << " seconds."
         << endl;
    cout << "-------------------------------------------------" << endl;
  }
}

void UltraweakDPGHestonSLV::PerformAdaptiveRefinement() {
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == 0) {
    cout << "Performing adaptive mesh refinement..." << endl;
  }

  // Compute error estimates
  Vector errors(mesh_ptr_->GetNE());
  ComputeErrorEstimates(errors);

  // Mark elements for refinement
  Array<int> marked_elements;
  double max_error = errors.Max();
  double threshold = threshold_factor_ * max_error;
  double total_error_sq = 0.0;

  for (int i = 0; i < errors.Size(); i++) {
    total_error_sq += errors(i) * errors(i);
  }

  double target_error_sq = error_fraction_ * total_error_sq;
  double accumulated_error_sq = 0.0;

  // Sort elements by error
  std::vector<std::pair<double, int>> error_indices(errors.Size());
  for (int i = 0; i < errors.Size(); i++) {
    error_indices[i] = std::make_pair(errors(i), i);
  }

  std::sort(error_indices.begin(), error_indices.end(),
            [](const std::pair<double, int> &a,
               const std::pair<double, int> &b) { return a.first > b.first; });

  // Accumulate error until we reach the target
  for (int i = 0; i < errors.Size(); i++) {
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

  if (myid == 0) {
    cout << "  Marked " << marked_elements.Size() << " elements for refinement."
         << endl;
  }

  // Refine the marked elements
  if (marked_elements.Size() > 0) {
    // Ensure mesh supports non-conforming refinement
    if (!mesh_ptr_->Nonconforming()) {
      mesh_ptr_->EnsureNCMesh(true);
    }

    // Apply refinement
    mesh_ptr_->GeneralRefinement(marked_elements);

    // Reinitialize spaces and solution
    InitializeSpaces();

    // Project solution from old to new spaces
    Vector old_sol(U_->Size());
    old_sol = *U_;

    // Project to new spaces
    for (int i = 0; i < Var::NUM_VARS; i++) {
      Vector old_block;
      U_prev_->GetBlock(i, old_block);

      int old_size = old_block.Size();
      int new_size = trial_fes_[i]->GetTrueVSize();

      if (old_size == new_size) {
        // Same size, direct copy
        U_->GetBlock(i) = old_block;
      } else {
        // Different size, need to project
        // For simplicity, we'll just initialize with zeros
        // A proper implementation would use projection operators
        U_->GetBlock(i) = 0.0;
      }
    }

    *U_prev_ = *U_;

    if (myid == 0) {
      cout << "  Mesh refined. New mesh size: " << mesh_ptr_->GetNE()
           << " elements." << endl;
    }
  }
}

void UltraweakDPGHestonSLV::ComputeErrorEstimates(Vector &errors) {
  // For each element, compute a local error estimate
  int ne = mesh_ptr_->GetNE();
  errors.SetSize(ne);
  errors = 0.0;

  // Create a simple error estimate based on the gradient of the solution
  ParGridFunction u_gf(trial_fes_[Var::PRICE]);
  u_gf.Distribute(U_->GetBlock(Var::PRICE));

  ParGridFunction sigma_s_gf(trial_fes_[Var::SIGMA_S]);
  sigma_s_gf.Distribute(U_->GetBlock(Var::SIGMA_S));

  ParGridFunction sigma_ss_gf(trial_fes_[Var::SIGMA_SS]);
  sigma_ss_gf.Distribute(U_->GetBlock(Var::SIGMA_SS));

  // Compute for local elements only
  for (int i = 0; i < mesh_ptr_->GetNE(); i++) {
    ElementTransformation *T = mesh_ptr_->GetElementTransformation(i);
    const FiniteElement *fe = trial_fes_[Var::PRICE]->GetFE(i);

    const IntegrationRule *ir =
        &IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder());

    double local_error = 0.0;

    for (int j = 0; j < ir->GetNPoints(); j++) {
      const IntegrationPoint &ip = ir->IntPoint(j);
      T->SetIntPoint(&ip);

      Vector pos;
      T->Transform(ip, pos);

      double S = pos(0);
      double V = pos(1);

      // Consider price, delta, and gamma in the error estimate
      double price_val = u_gf.GetValue(*T, ip);
      double delta_val = sigma_s_gf.GetValue(*T, ip);
      double gamma_val = sigma_ss_gf.GetValue(*T, ip);

      // Weight errors more near the strike price
      double weight =
          1.0 + 5.0 * exp(-0.5 * pow((S - params.K) / (0.2 * params.K), 2));

      // Higher weight for small V (more important region)
      weight *= 1.0 + 0.5 / (V + 0.05);

      // Combine error measures (this is a heuristic approach)
      local_error += weight * (fabs(gamma_val) + fabs(delta_val)) * ip.weight;
    }

    errors(i) = local_error;
  }
}

double UltraweakDPGHestonSLV::ComputeTotalError() {
  Vector errors;
  ComputeErrorEstimates(errors);

  double total_error = 0.0;
  for (int i = 0; i < errors.Size(); i++) {
    total_error += errors(i) * errors(i);
  }

  return sqrt(total_error);
}

const ParGridFunction *
UltraweakDPGHestonSLV::GetSolutionGridFunction(Var v) const {
  if (v < 0 || v >= Var::NUM_VARS || !trial_fes_[v] || !U_) {
    return nullptr;
  }

  // Create a new ParGridFunction (caller is responsible for deleting it)
  ParGridFunction *gf = new ParGridFunction(trial_fes_[v]);
  gf->Distribute(U_->GetBlock(v));

  return gf;
}

const BlockVector *UltraweakDPGHestonSLV::GetSolutionVector() const {
  return U_;
}

const ParMesh *UltraweakDPGHestonSLV::GetMesh() const { return mesh_ptr_; }

void UltraweakDPGHestonSLV::SaveSolution(const std::string &filename_prefix) {
  if (!U_ || !mesh_ptr_) {
    cout << "Warning: Cannot save solution - not initialized." << endl;
    return;
  }

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == 0) {
    cout << "Saving solution to " << filename_prefix << "..." << endl;
  }

  // Use ParaView data collection for parallel output
  ParaViewDataCollection paraview_dc(filename_prefix.c_str(), mesh_ptr_);
  paraview_dc.SetLevelsOfDetail(trial_order);
  paraview_dc.SetCycle(0);
  paraview_dc.SetTime(0.0);

  // Save each variable
  const char *var_names[Var::NUM_VARS] = {
      "Price", "Delta", "Vega", "Gamma", "Cross_Derivative", "Volga"};

  for (int i = 0; i < Var::NUM_VARS; i++) {
    ParGridFunction gf(trial_fes_[i]);
    gf.Distribute(U_->GetBlock(i));
    paraview_dc.RegisterField(var_names[i], &gf);
  }

  paraview_dc.Save();

  if (myid == 0) {
    cout << "Solution saved successfully." << endl;
  }
}

void UltraweakDPGHestonSLV::SaveMesh(const std::string &filename) {
  if (!mesh_ptr_) {
    cout << "Warning: Cannot save mesh - not initialized." << endl;
    return;
  }

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == 0) {
    cout << "Saving mesh to " << filename << "..." << endl;
    std::ofstream mesh_ofs(filename.c_str());
    mesh_ptr_->PrintAsSerial(mesh_ofs);
    mesh_ofs.close();
    cout << "Mesh saved successfully." << endl;
  }
}