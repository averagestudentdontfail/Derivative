#ifndef DPG_ERROR_ESTIMATOR_HPP
#define DPG_ERROR_ESTIMATOR_HPP

#include "mfem.hpp"

namespace mfem
{

// Forward declarations
class BilinearForm;
class MixedBilinearForm;
class LinearForm;
class GridFunction;
class FiniteElementSpace;

class UltraweakDPGErrorEstimator
{
private:
    // Pointers to DPG components (will be set later)
    FiniteElementSpace *trial_fes_;
    FiniteElementSpace *test_fes_;
    MixedBilinearForm  *b_form_; // B
    BilinearForm       *g_form_; // G (test inner product)
    LinearForm         *f_form_; // F (RHS)
    GridFunction       *u_sol_;  // Solution U

    Vector local_errors_; // Cached errors
    bool errors_computed_;

public:
    UltraweakDPGErrorEstimator(FiniteElementSpace *trial_fes = nullptr,
                               FiniteElementSpace *test_fes = nullptr,
                               MixedBilinearForm *b = nullptr,
                               BilinearForm *g = nullptr,
                               LinearForm *f = nullptr,
                               GridFunction *u = nullptr);

    virtual ~UltraweakDPGErrorEstimator() = default;

    // Set components (alternative to constructor)
    void SetTrialFESpace(FiniteElementSpace *fes) { trial_fes_ = fes; Reset(); }
    void SetTestFESpace(FiniteElementSpace *fes)  { test_fes_ = fes; Reset(); }
    void SetBForm(MixedBilinearForm *b) { b_form_ = b; Reset(); }
    void SetGForm(BilinearForm *g)      { g_form_ = g; Reset(); }
    void SetFForm(LinearForm *f)        { f_form_ = f; Reset(); }
    void SetSolution(GridFunction *u)   { u_sol_ = u; Reset(); }


    // Main methods (stubs for Phase 1)
    virtual const Vector &GetLocalErrors();
    virtual void Reset();
    virtual double GetLocalError(int elem_idx);

private:
    // Internal computation method
    virtual void ComputeEstimator();
};

// Placeholder for Refiner/Derefiner classes for Phase 1
class DPGAdaptiveRefiner {
    // To be implemented in Phase 1
};
class DPGAdaptiveDerefiner {
    // To be implemented in Phase 1
};


} // namespace mfem

#endif // DPG_ERROR_ESTIMATOR_HPP