#ifndef DPG_INTEGRATORS_HPP
#define DPG_INTEGRATORS_HPP

#include "mfem.hpp"

namespace mfem
{

// Base class (optional but good practice)
class UltraweakDPGIntegrator : public BilinearFormIntegrator
{
public:
    UltraweakDPGIntegrator() = default;
    virtual ~UltraweakDPGIntegrator() = default;
};

// Integrator for (u, div w) - trial u (scalar), test w (vector)
// Used in sigma - grad u = 0 -> (sigma, w) + (u, div w) - <{u}, [w.n]> = 0
class MixedScalarWeakGradientIntegrator : public BilinearFormIntegrator
{
public:
    MixedScalarWeakGradientIntegrator() = default;

    virtual void AssembleElementMatrix(const FiniteElement &trial_fe, // scalar u
                                     const FiniteElement &test_fe,  // vector w
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat) override;
};


// Integrator for <{u}, [w.n]> on interior faces and <u, w.n> on boundary
// Trial u (scalar), Test w (vector)
class UltraweakJumpIntegrator : public BilinearFormIntegrator
{
public:
    UltraweakJumpIntegrator() = default;

    // Interior faces
    virtual void AssembleFaceMatrix(const FiniteElement &trial_fe1, // u on K1
                                  const FiniteElement &trial_fe2, // u on K2
                                  const FiniteElement &test_fe1,  // w on K1
                                  const FiniteElement &test_fe2,  // w on K2
                                  FaceElementTransformations &Trans,
                                  DenseMatrix &elmat) override;

    // Boundary faces
    virtual void AssembleBoundaryFaceMatrix(const FiniteElement &trial_fe, // u on K
                                          const FiniteElement &test_fe,  // w on K
                                          FaceElementTransformations &Trans,
                                          DenseMatrix &elmat) override;
};


// Integrator for test space inner product G = (v,v') + (w,w')
// For Poisson: v is scalar, w is vector. Combine L2 norms for simplicity now.
class TestInnerProductIntegrator : public BilinearFormIntegrator
{
private:
    double alpha_mass_ = 1.0; // Coefficient for mass term (L2)
    // Add alpha_diffusion_ later if H1 norm is needed

public:
    TestInnerProductIntegrator(double alpha_mass = 1.0);

    virtual void AssembleElementMatrix(const FiniteElement &el,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat) override;
};


} // namespace mfem

#endif // DPG_INTEGRATORS_HPP