#include "dpg-integrators.hpp"

namespace mfem
{

// --- MixedScalarWeakGradientIntegrator --- (u, div w)
void MixedScalarWeakGradientIntegrator::AssembleElementMatrix(
    const FiniteElement &trial_fe, const FiniteElement &test_fe,
    ElementTransformation &Trans, DenseMatrix &elmat)
{
    int trial_dof = trial_fe.GetDof(); // u
    int test_dof = test_fe.GetDof();   // w
    int dim = Trans.GetSpaceDim();

    // MFEM_VERIFY(test_fe.GetRangeType() == FiniteElement::VECTOR,
    //             "Test FE must be vector type for Weak Gradient integrator.");
    // MFEM_VERIFY(trial_fe.GetRangeType() == FiniteElement::SCALAR,
    //             "Trial FE must be scalar type for Weak Gradient integrator.");

    Vector trial_shape(trial_dof);
    DenseMatrix test_dshape(test_dof, dim); // Dof x Dim
    Vector div_w(test_dof);

    elmat.SetSize(test_dof, trial_dof); // test_dof = rows, trial_dof = cols
    elmat = 0.0;

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = trial_fe.GetOrder() + test_fe.GetOrder(); // Sufficient order
        ir = &IntRules.Get(trial_fe.GetGeomType(), order);
    }

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);
        double w = ip.weight * Trans.Weight();

        trial_fe.CalcShape(ip, trial_shape); // u basis functions
        test_fe.CalcDivShape(ip, div_w);     // div(w) basis functions

        // elmat(j, k) += w * div_w(j) * trial_shape(k)
        AddMult_a_VWt(w, div_w, trial_shape, elmat);
    }
}

// --- UltraweakJumpIntegrator --- <{u}, [w.n]> and <u, w.n>
void UltraweakJumpIntegrator::AssembleFaceMatrix(
    const FiniteElement &trial_fe1, const FiniteElement &trial_fe2,
    const FiniteElement &test_fe1, const FiniteElement &test_fe2,
    FaceElementTransformations &Trans, DenseMatrix &elmat)
{
    int trial_dof1 = trial_fe1.GetDof(); // u on K1
    int trial_dof2 = trial_fe2.GetDof(); // u on K2
    int test_dof1 = test_fe1.GetDof();   // w on K1
    int test_dof2 = test_fe2.GetDof();   // w on K2
    int dim = Trans.GetSpaceDim();

    // MFEM_VERIFY(test_fe1.GetRangeType() == FiniteElement::VECTOR &&
    //             test_fe2.GetRangeType() == FiniteElement::VECTOR,
    //             "Test FEs must be vector type for Jump integrator.");
    // MFEM_VERIFY(trial_fe1.GetRangeType() == FiniteElement::SCALAR &&
    //             trial_fe2.GetRangeType() == FiniteElement::SCALAR,
    //             "Trial FEs must be scalar type for Jump integrator.");

    Vector trial_shape1(trial_dof1);
    Vector trial_shape2(trial_dof2);
    DenseMatrix test_shape1(test_dof1, dim); // Dof x Dim for vector FE
    DenseMatrix test_shape2(test_dof2, dim);

    // Result matrix: (test_dof1 + test_dof2) x (trial_dof1 + trial_dof2)
    elmat.SetSize(test_dof1 + test_dof2, trial_dof1 + trial_dof2);
    elmat = 0.0;

    // Pointers to submatrices of elmat
    DenseMatrix elmat11(elmat.Data(), test_dof1, trial_dof1);
    DenseMatrix elmat12(elmat.Data() + test_dof1 * trial_dof1, test_dof1, trial_dof2);
    DenseMatrix elmat21(elmat.Data() + test_dof1 * (trial_dof1 + trial_dof2), test_dof2, trial_dof1);
    DenseMatrix elmat22(elmat.Data() + test_dof1 * (trial_dof1 + trial_dof2) + test_dof2 * trial_dof1, test_dof2, trial_dof2);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = trial_fe1.GetOrder() + test_fe1.GetOrder();
        order = std::max(order, trial_fe2.GetOrder() + test_fe2.GetOrder());
        ir = &IntRules.Get(Trans.GetGeometryType(), order);
    }

    Vector normal(dim);
    Vector w_n(std::max(test_dof1, test_dof2)); // w . n

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.Face->SetIntPoint(&ip);
        CalcOrtho(Trans.Face->Jacobian(), normal); // Normal vector
        double w = ip.weight * Trans.Face->Weight();

        // Get basis functions evaluated on the face
        IntegrationPoint eip1, eip2;
        Trans.Loc1.Transform(ip, eip1);
        Trans.Loc2.Transform(ip, eip2);

        trial_fe1.CalcShape(eip1, trial_shape1); // u basis on K1 side
        trial_fe2.CalcShape(eip2, trial_shape2); // u basis on K2 side
        test_fe1.CalcVShape(Trans.Loc1, eip1, test_shape1); // w basis on K1 side
        test_fe2.CalcVShape(Trans.Loc2, eip2, test_shape2); // w basis on K2 side


        // Term <{u}, [w.n]> = <0.5*(u1+u2), w1.n1 + w2.n2>
        // Note: n1 = normal, n2 = -normal

        // Contribution to elmat11: <0.5*u1, w1.n>
        w_n.SetSize(test_dof1);
        test_shape1.MultTranspose(normal, w_n); // w1.n
        AddMult_a_VWt(0.5 * w, w_n, trial_shape1, elmat11);

        // Contribution to elmat12: <0.5*u2, w1.n>
        AddMult_a_VWt(0.5 * w, w_n, trial_shape2, elmat12);

        // Contribution to elmat21: <0.5*u1, w2.(-n)> = <-0.5*u1, w2.n>
        w_n.SetSize(test_dof2);
        test_shape2.MultTranspose(normal, w_n); // w2.n
        AddMult_a_VWt(-0.5 * w, w_n, trial_shape1, elmat21);

        // Contribution to elmat22: <0.5*u2, w2.(-n)> = <-0.5*u2, w2.n>
        AddMult_a_VWt(-0.5 * w, w_n, trial_shape2, elmat22);

    }
}

void UltraweakJumpIntegrator::AssembleBoundaryFaceMatrix(
    const FiniteElement &trial_fe, const FiniteElement &test_fe,
    FaceElementTransformations &Trans, DenseMatrix &elmat)
{
    int trial_dof = trial_fe.GetDof(); // u
    int test_dof = test_fe.GetDof();   // w
    int dim = Trans.GetSpaceDim();

    // MFEM_VERIFY(test_fe.GetRangeType() == FiniteElement::VECTOR,
    //             "Test FE must be vector type for Boundary Jump integrator.");
    // MFEM_VERIFY(trial_fe.GetRangeType() == FiniteElement::SCALAR,
    //             "Trial FE must be scalar type for Boundary Jump integrator.");

    Vector trial_shape(trial_dof);
    DenseMatrix test_shape(test_dof, dim); // Dof x Dim for vector FE

    elmat.SetSize(test_dof, trial_dof);
    elmat = 0.0;

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = trial_fe.GetOrder() + test_fe.GetOrder();
        ir = &IntRules.Get(Trans.GetGeometryType(), order);
    }

    Vector normal(dim);
    Vector w_n(test_dof); // w . n

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.Face->SetIntPoint(&ip);
        CalcOrtho(Trans.Face->Jacobian(), normal); // Normal vector
        double w = ip.weight * Trans.Face->Weight();

        // Get basis functions evaluated on the face
        IntegrationPoint eip;
        Trans.Loc1.Transform(ip, eip);

        trial_fe.CalcShape(eip, trial_shape);       // u basis
        test_fe.CalcVShape(Trans.Loc1, eip, test_shape); // w basis

        // Term <u, w.n>
        w_n.SetSize(test_dof);
        test_shape.MultTranspose(normal, w_n); // w.n
        AddMult_a_VWt(w, w_n, trial_shape, elmat);
    }
}


// --- TestInnerProductIntegrator --- (v,v') + (w,w') L2 norms
TestInnerProductIntegrator::TestInnerProductIntegrator(double alpha_mass)
    : alpha_mass_(alpha_mass) {}

void TestInnerProductIntegrator::AssembleElementMatrix(
    const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
{
    int ndofs = el.GetDof();
    int dim = el.GetDim();
    bool is_scalar = (el.GetRangeType() == FiniteElement::SCALAR);
    int vdim = is_scalar ? 1 : dim; // Vector dimension

    Vector shape_scalar;
    DenseMatrix shape_vector; // Dof x vdim

    if (is_scalar) {
        shape_scalar.SetSize(ndofs);
    } else {
        shape_vector.SetSize(ndofs, vdim);
    }

    elmat.SetSize(ndofs);
    elmat = 0.0;

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr) {
        int order = 2 * el.GetOrder();
        ir = &IntRules.Get(el.GetGeomType(), order);
    }

    DenseMatrix M_loc(ndofs); // Local mass matrix

    for (int i = 0; i < ir->GetNPoints(); i++) {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);
        double w = ip.weight * Trans.Weight() * alpha_mass_;

        if (is_scalar) {
            el.CalcShape(ip, shape_scalar);
            M_loc = 0.0;
            AddMult_a_VVt(1.0, shape_scalar, M_loc); // M_loc = shape * shape^T
        } else {
            el.CalcVShape(Trans, ip, shape_vector); // Dof x vdim
            // We want sum_d (shape_d * shape_d^T)
            M_loc = 0.0;
            Vector shape_d(ndofs);
            for(int d=0; d<vdim; ++d) {
                shape_vector.GetColumn(d, shape_d);
                AddMult_a_VVt(1.0, shape_d, M_loc);
            }
        }
        elmat.Add(w, M_loc);
    }
    // Add H1 term later if needed: elmat.Add(alpha_diffusion_, K_loc);
}

} // namespace mfem