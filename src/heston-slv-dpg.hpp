#ifndef HESTON_SLV_DPG_HPP
#define HESTON_SLV_DPG_HPP

#include "mfem.hpp"
#include <vector>
#include <memory>

namespace mfem
{

// Forward declarations if needed later
class HestonSLVParameters;
class HestonSLVCoefficients;

class UltraweakDPGHestonSLV
{
public:
    // Enum for variable indexing (will be expanded)
    enum Var { PRICE = 0, NUM_VARS = 1 }; // Start minimal

private:
    // Basic members needed later
    mfem::ParMesh *mesh_ptr_;
    int trial_order_;
    int test_order_;

    // FE Spaces (will be vectors later)
    mfem::ParFiniteElementSpace *trial_fes_;
    mfem::ParFiniteElementSpace *test_fes_;

    // Solution vector (will be BlockVector later)
    mfem::ParGridFunction *U_;

public:
    UltraweakDPGHestonSLV(int trial_order = 1, int test_enrichment = 1);
    virtual ~UltraweakDPGHestonSLV();

    // Basic methods (stubs for now)
    void Initialize();
    void Solve();

    // Add other configuration methods later...
};

} // namespace mfem

#endif // HESTON_SLV_DPG_HPP