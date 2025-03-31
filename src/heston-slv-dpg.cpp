#include "heston-slv-dpg.hpp"
#include <iostream>

namespace mfem
{

UltraweakDPGHestonSLV::UltraweakDPGHestonSLV(int trial_order, int test_enrichment)
    : mesh_ptr_(nullptr), trial_order_(trial_order),
      test_order_(trial_order + test_enrichment),
      trial_fes_(nullptr), test_fes_(nullptr), U_(nullptr)
{
    // Ensure test order is higher
    if (test_order_ <= trial_order_)
    {
        test_order_ = trial_order_ + 1;
        std::cout << "Warning: Test order enrichment must be at least 1. Setting test_order to "
                  << test_order_ << std::endl;
    }
    std::cout << "UltraweakDPGHestonSLV created (minimal)." << std::endl;
}

UltraweakDPGHestonSLV::~UltraweakDPGHestonSLV()
{
    delete U_;
    delete test_fes_;
    delete trial_fes_;
    delete mesh_ptr_; // We will own the mesh
    std::cout << "UltraweakDPGHestonSLV destroyed (minimal)." << std::endl;
}

void UltraweakDPGHestonSLV::Initialize()
{
    std::cout << "UltraweakDPGHestonSLV::Initialize() called (stub)." << std::endl;
    // Mesh, spaces, solution vector initialization will go here
}

void UltraweakDPGHestonSLV::Solve()
{
    std::cout << "UltraweakDPGHestonSLV::Solve() called (stub)." << std::endl;
    // Time stepping loop will go here
}

} // namespace mfem