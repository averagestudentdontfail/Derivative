#include "dpg-integrators.hpp"
#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// UltraweakDPGIntegrator
UltraweakDPGIntegrator::UltraweakDPGIntegrator(Coefficient *q) : Q(q) {}

UltraweakDPGIntegrator::~UltraweakDPGIntegrator() {}

// UltraweakDerivativeIntegrator
UltraweakDerivativeIntegrator::UltraweakDerivativeIntegrator(int dir,
                                                             Coefficient *q)
    : UltraweakDPGIntegrator(q), di(dir) {}

void UltraweakDerivativeIntegrator::AssembleElementMatrix(
    const FiniteElement &trial_fe, const FiniteElement &test_fe,
    ElementTransformation &Trans, DenseMatrix &elmat) {
  int trial_ndofs = trial_fe.GetDof();
  int test_ndofs = test_fe.GetDof();
  int dim = trial_fe.GetDim();

  MFEM_VERIFY(di >= 0 && di < dim, "Invalid derivative direction");

  Vector shape(trial_ndofs);
  DenseMatrix dshape(test_ndofs, dim);
  Vector dshape_i(test_ndofs);

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = trial_fe.GetOrder() + test_fe.GetOrder();
    ir = &IntRules.Get(trial_fe.GetGeomType(), order);
  }

  elmat.SetSize(test_ndofs, trial_ndofs);
  elmat = 0.0;

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    Trans.SetIntPoint(&ip);

    trial_fe.CalcShape(ip, shape);
    test_fe.CalcDShape(ip, dshape);

    for (int j = 0; j < test_ndofs; j++) {
      dshape_i(j) = dshape(j, di);
    }

    double coef = 1.0;
    if (Q) {
      coef = Q->Eval(Trans, ip);
    }

    double w = ip.weight * coef * Trans.Weight();

    for (int j = 0; j < test_ndofs; j++) {
      for (int k = 0; k < trial_ndofs; k++) {
        elmat(j, k) += dshape_i(j) * shape(k) * w;
      }
    }
  }
}

// UltraweakSecondDerivativeIntegrator implementation
UltraweakSecondDerivativeIntegrator::UltraweakSecondDerivativeIntegrator(
    int dir1, int dir2, Coefficient *q)
    : UltraweakDPGIntegrator(q), di(dir1), dj(dir2) {}

void UltraweakSecondDerivativeIntegrator::AssembleElementMatrix(
    const FiniteElement &trial_fe, const FiniteElement &test_fe,
    ElementTransformation &Trans, DenseMatrix &elmat) {
  int trial_ndofs = trial_fe.GetDof();
  int test_ndofs = test_fe.GetDof();
  int dim = trial_fe.GetDim();

  MFEM_VERIFY(di >= 0 && di < dim && dj >= 0 && dj < dim,
              "Invalid derivative directions");

  Vector shape_trial(trial_ndofs);
  DenseMatrix dshape(test_ndofs, dim);
  DenseMatrix d2shape;

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = trial_fe.GetOrder() + test_fe.GetOrder();
    ir = &IntRules.Get(trial_fe.GetGeomType(), order);
  }

  elmat.SetSize(test_ndofs, trial_ndofs);
  elmat = 0.0;

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    Trans.SetIntPoint(&ip);

    trial_fe.CalcShape(ip, shape_trial);

    if (di == dj) {
      // Second derivative in same direction
      d2shape.SetSize(test_ndofs, dim * dim);
      test_fe.CalcHessian(ip, d2shape);

      // Extract the specific second derivative
      Vector d2shape_ij(test_ndofs);
      for (int j = 0; j < test_ndofs; j++) {
        d2shape_ij(j) = d2shape(j, di * dim + dj);
      }

      double coef = 1.0;
      if (Q) {
        coef = Q->Eval(Trans, ip);
      }

      double w = ip.weight * coef * Trans.Weight();

      for (int j = 0; j < test_ndofs; j++) {
        for (int k = 0; k < trial_ndofs; k++) {
          elmat(j, k) += d2shape_ij(j) * shape_trial(k) * w;
        }
      }
    } else {
      // Mixed derivative
      d2shape.SetSize(test_ndofs, dim * dim);
      test_fe.CalcHessian(ip, d2shape);

      // Extract the mixed derivative
      Vector d2shape_ij(test_ndofs);
      for (int j = 0; j < test_ndofs; j++) {
        d2shape_ij(j) = d2shape(j, di * dim + dj);
      }

      double coef = 1.0;
      if (Q) {
        coef = Q->Eval(Trans, ip);
      }

      double w = ip.weight * coef * Trans.Weight();

      for (int j = 0; j < test_ndofs; j++) {
        for (int k = 0; k < trial_ndofs; k++) {
          elmat(j, k) += d2shape_ij(j) * shape_trial(k) * w;
        }
      }
    }
  }
}

// UltraweakTraceIntegrator implementation
UltraweakTraceIntegrator::UltraweakTraceIntegrator(int dir, double alpha_)
    : di(dir), alpha(alpha_) {}

void UltraweakTraceIntegrator::AssembleFaceMatrix(
    const FiniteElement &trial_fe1, const FiniteElement &trial_fe2,
    const FiniteElement &test_fe1, const FiniteElement &test_fe2,
    FaceElementTransformations &Trans, DenseMatrix &elmat) {
  int trial1_ndofs = trial_fe1.GetDof();
  int trial2_ndofs = trial_fe2.GetDof();
  int test1_ndofs = test_fe1.GetDof();
  int test2_ndofs = test_fe2.GetDof();
  int dim = trial_fe1.GetDim();

  MFEM_VERIFY(di >= 0 && di < dim, "Invalid derivative direction");

  Vector shape1(trial1_ndofs);
  Vector shape2(trial2_ndofs);
  Vector test_shape1(test1_ndofs);
  Vector test_shape2(test2_ndofs);

  elmat.SetSize((test1_ndofs + test2_ndofs), (trial1_ndofs + trial2_ndofs));
  elmat = 0.0;

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = max(trial_fe1.GetOrder() + test_fe1.GetOrder(),
                    trial_fe2.GetOrder() + test_fe2.GetOrder());
    ir = &IntRules.Get(Trans.GetGeometryType(), order);
  }

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    IntegrationPoint eip1, eip2;

    Trans.Loc1.Transform(ip, eip1);
    Trans.Loc2.Transform(ip, eip2);

    trial_fe1.CalcShape(eip1, shape1);
    trial_fe2.CalcShape(eip2, shape2);
    test_fe1.CalcShape(eip1, test_shape1);
    test_fe2.CalcShape(eip2, test_shape2);

    Trans.Face->SetIntPoint(&ip);

    Vector nor(dim);
    CalcOrtho(Trans.Face->Jacobian(), nor);
    nor /= nor.Norml2(); // Unit normal

    double ni = nor(di); // Normal component in the derivative direction
    double w = ip.weight * Trans.Face->Weight();

    // Assemble the four sub-blocks of the element matrix

    // Block (test1, trial1): Average flux term
    for (int j = 0; j < test1_ndofs; j++) {
      for (int k = 0; k < trial1_ndofs; k++) {
        elmat(j, k) += 0.5 * w * ni * test_shape1(j) * shape1(k);
      }
    }

    // Block (test1, trial2): Average flux term
    for (int j = 0; j < test1_ndofs; j++) {
      for (int k = 0; k < trial2_ndofs; k++) {
        elmat(j, trial1_ndofs + k) -= 0.5 * w * ni * test_shape1(j) * shape2(k);
      }
    }

    // Block (test2, trial1): Average flux term
    for (int j = 0; j < test2_ndofs; j++) {
      for (int k = 0; k < trial1_ndofs; k++) {
        elmat(test1_ndofs + j, k) -= 0.5 * w * ni * test_shape2(j) * shape1(k);
      }
    }

    // Block (test2, trial2): Average flux term
    for (int j = 0; j < test2_ndofs; j++) {
      for (int k = 0; k < trial2_ndofs; k++) {
        elmat(test1_ndofs + j, trial1_ndofs + k) +=
            0.5 * w * ni * test_shape2(j) * shape2(k);
      }
    }

    // Add upwinding term if needed
    if (alpha != 0.0) {
      double upwind_term = 0.5 * alpha * fabs(ni);

      for (int j = 0; j < test1_ndofs; j++) {
        for (int k = 0; k < trial1_ndofs; k++) {
          elmat(j, k) += upwind_term * w * test_shape1(j) * shape1(k);
        }
        for (int k = 0; k < trial2_ndofs; k++) {
          elmat(j, trial1_ndofs + k) -=
              upwind_term * w * test_shape1(j) * shape2(k);
        }
      }

      for (int j = 0; j < test2_ndofs; j++) {
        for (int k = 0; k < trial1_ndofs; k++) {
          elmat(test1_ndofs + j, k) -=
              upwind_term * w * test_shape2(j) * shape1(k);
        }
        for (int k = 0; k < trial2_ndofs; k++) {
          elmat(test1_ndofs + j, trial1_ndofs + k) +=
              upwind_term * w * test_shape2(j) * shape2(k);
        }
      }
    }
  }
}

// UltraweakBoundaryIntegrator implementation
UltraweakBoundaryIntegrator::UltraweakBoundaryIntegrator(int dir,
                                                         Coefficient *bdr_coef_)
    : di(dir), bdr_coef(bdr_coef_) {}

void UltraweakBoundaryIntegrator::AssembleFaceMatrix(
    const FiniteElement &trial_fe, const FiniteElement &test_fe,
    FaceElementTransformations &Trans, DenseMatrix &elmat) {
  int trial_ndofs = trial_fe.GetDof();
  int test_ndofs = test_fe.GetDof();
  int dim = trial_fe.GetDim();

  MFEM_VERIFY(di >= 0 && di < dim, "Invalid derivative direction");

  Vector shape(trial_ndofs);
  Vector test_shape(test_ndofs);

  elmat.SetSize(test_ndofs, trial_ndofs);
  elmat = 0.0;

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = trial_fe.GetOrder() + test_fe.GetOrder();
    ir = &IntRules.Get(Trans.GetGeometryType(), order);
  }

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    IntegrationPoint eip;

    Trans.Loc1.Transform(ip, eip);
    Trans.Face->SetIntPoint(&ip);

    trial_fe.CalcShape(eip, shape);
    test_fe.CalcShape(eip, test_shape);

    Vector nor(dim);
    CalcOrtho(Trans.Face->Jacobian(), nor);
    nor /= nor.Norml2(); // Unit normal

    double ni = nor(di); // Normal component in the derivative direction

    double coef = 1.0;
    if (bdr_coef) {
      coef = bdr_coef->Eval(*Trans.Face, ip);
    }

    double w = ip.weight * coef * ni * Trans.Face->Weight();

    for (int j = 0; j < test_ndofs; j++) {
      for (int k = 0; k < trial_ndofs; k++) {
        elmat(j, k) += test_shape(j) * shape(k) * w;
      }
    }
  }
}

// UltraweakBoundaryLFIntegrator implementation
UltraweakBoundaryLFIntegrator::UltraweakBoundaryLFIntegrator(
    int dir, FunctionCoefficient *bdr_func_)
    : di(dir), bdr_func(bdr_func_) {}

void UltraweakBoundaryLFIntegrator::AssembleFaceVector(
    const FiniteElement &el, FaceElementTransformations &Trans,
    Vector &elvect) {
  int dofs = el.GetDof();
  int dim = el.GetDim();

  MFEM_VERIFY(di >= 0 && di < dim, "Invalid derivative direction");

  Vector shape(dofs);

  elvect.SetSize(dofs);
  elvect = 0.0;

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = 2 * el.GetOrder();
    ir = &IntRules.Get(Trans.GetGeometryType(), order);
  }

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    IntegrationPoint eip;

    Trans.Loc1.Transform(ip, eip);
    el.CalcShape(eip, shape);

    Trans.Face->SetIntPoint(&ip);

    Vector nor(dim);
    CalcOrtho(Trans.Face->Jacobian(), nor);
    nor /= nor.Norml2(); // Unit normal

    double ni = nor(di); // Normal component in the derivative direction

    double val = 0.0;
    if (bdr_func) {
      val = bdr_func->Eval(*Trans.Face, ip);
    }

    double w = ip.weight * val * ni * Trans.Face->Weight();

    shape *= w;
    elvect += shape;
  }
}

// UltraweakJumpIntegrator implementation
UltraweakJumpIntegrator::UltraweakJumpIntegrator(int dir) : direction(dir) {}

void UltraweakJumpIntegrator::AssembleFaceMatrix(
    const FiniteElement &trial_fe1, const FiniteElement &trial_fe2,
    const FiniteElement &test_fe1, const FiniteElement &test_fe2,
    FaceElementTransformations &Trans, DenseMatrix &elmat) {
  // Assuming the test space is conforming, test_fe1 and test_fe2 are the same
  // FE evaluated using the respective element transformations.
  MFEM_ASSERT(&test_fe1 == &test_fe2,
              "UltraweakJumpIntegrator assumes conforming test space (test_fe1 "
              "== test_fe2)");

  int trial1_ndofs = trial_fe1.GetDof();
  int trial2_ndofs = trial_fe2.GetDof();
  int test_ndofs = test_fe1.GetDof();
  int dim = trial_fe1.GetDim();

  MFEM_VERIFY(direction >= 0 && direction < dim, "Invalid jump direction");

  Vector shape1(trial1_ndofs);
  Vector shape2(trial2_ndofs);
  Vector test_shape1(test_ndofs);
  Vector test_shape2(test_ndofs);

  elmat.SetSize(2 * test_ndofs, trial1_ndofs + trial2_ndofs);
  elmat = 0.0;

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = max(trial_fe1.GetOrder() + test_fe1.GetOrder(),
                    trial_fe2.GetOrder() + test_fe2.GetOrder());
    ir = &IntRules.Get(Trans.GetGeometryType(), order);
  }

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    IntegrationPoint eip1, eip2;

    Trans.Loc1.Transform(ip, eip1);
    Trans.Loc2.Transform(ip, eip2);

    trial_fe1.CalcShape(eip1, shape1);
    trial_fe2.CalcShape(eip2, shape2);
    test_fe1.CalcShape(eip1, test_shape1);
    test_fe2.CalcShape(eip2, test_shape2);

    Trans.Face->SetIntPoint(&ip);

    Vector nor(dim);
    CalcOrtho(Trans.Face->Jacobian(), nor);
    nor /= nor.Norml2(); // Unit normal

    double n_dir = nor(direction); // Normal component in jump direction
    double w = ip.weight * Trans.Face->Weight();

    // Jump term [[v]] = v^+ * n^+ + v^- * n^-
    // For element 1 (v^+)
    for (int j = 0; j < test_ndofs; j++) {
      for (int k = 0; k < trial1_ndofs; k++) {
        elmat(j, k) += w * n_dir * test_shape1(j) * shape1(k);
      }
    }

    // For element 2 (v^-)
    for (int j = 0; j < test_ndofs; j++) {
      for (int k = 0; k < trial2_ndofs; k++) {
        elmat(test_ndofs + j, trial1_ndofs + k) -=
            w * n_dir * test_shape2(j) * shape2(k);
      }
    }
  }
}

void UltraweakJumpIntegrator::AssembleFaceMatrix(
    const FiniteElement &trial_fe1, const FiniteElement &trial_fe2,
    const FiniteElement &test_fe, FaceElementTransformations &Trans,
    DenseMatrix &elmat) {
  // This is assuming test_fe is shared between the two elements
  int trial1_ndofs = trial_fe1.GetDof();
  int trial2_ndofs = trial_fe2.GetDof();
  int test_ndofs = test_fe.GetDof();
  int dim = trial_fe1.GetDim();

  MFEM_VERIFY(direction >= 0 && direction < dim, "Invalid jump direction");

  Vector shape1(trial1_ndofs);
  Vector shape2(trial2_ndofs);
  Vector test_shape(test_ndofs);

  elmat.SetSize(test_ndofs, trial1_ndofs + trial2_ndofs);
  elmat = 0.0;

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = max(trial_fe1.GetOrder() + test_fe.GetOrder(),
                    trial_fe2.GetOrder() + test_fe.GetOrder());
    ir = &IntRules.Get(Trans.GetGeometryType(), order);
  }

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    IntegrationPoint eip1, eip2;

    Trans.Loc1.Transform(ip, eip1);
    Trans.Loc2.Transform(ip, eip2);

    trial_fe1.CalcShape(eip1, shape1);
    trial_fe2.CalcShape(eip2, shape2);
    test_fe.CalcShape(eip1, test_shape);

    Trans.Face->SetIntPoint(&ip);

    Vector nor(dim);
    CalcOrtho(Trans.Face->Jacobian(), nor);
    nor /= nor.Norml2(); // Unit normal

    double n_dir = nor(direction); // Normal component in jump direction
    double w = ip.weight * Trans.Face->Weight();

    // Jump term [[v]] = v^+ * n^+ + v^- * n^-
    // For element 1 (v^+)
    for (int j = 0; j < test_ndofs; j++) {
      for (int k = 0; k < trial1_ndofs; k++) {
        elmat(j, k) += w * n_dir * test_shape(j) * shape1(k);
      }
    }

    // For element 2 (v^-)
    for (int j = 0; j < test_ndofs; j++) {
      for (int k = 0; k < trial2_ndofs; k++) {
        elmat(j, trial1_ndofs + k) -= w * n_dir * test_shape(j) * shape2(k);
      }
    }
  }
}

// UltraweakAverageIntegrator implementation
UltraweakAverageIntegrator::UltraweakAverageIntegrator(int dir, double scale)
    : direction(dir), scaling(scale) {}

void UltraweakAverageIntegrator::AssembleFaceMatrix(
    const FiniteElement &trial_fe1, const FiniteElement &trial_fe2,
    const FiniteElement &test_fe1, const FiniteElement &test_fe2,
    FaceElementTransformations &Trans, DenseMatrix &elmat) {
  int trial1_ndofs = trial_fe1.GetDof();
  int trial2_ndofs = trial_fe2.GetDof();
  int test1_ndofs = test_fe1.GetDof();
  int test2_ndofs = test_fe2.GetDof();
  int dim = Trans.GetSpaceDim();

  Vector shape1(trial1_ndofs);
  Vector shape2(trial2_ndofs);
  Vector test_shape1(test1_ndofs);
  Vector test_shape2(test2_ndofs);

  elmat.SetSize(test1_ndofs + test2_ndofs, trial1_ndofs + trial2_ndofs);
  elmat = 0.0;

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = max(trial_fe1.GetOrder() + test_fe1.GetOrder(),
                    trial_fe2.GetOrder() + test_fe2.GetOrder());
    ir = &IntRules.Get(Trans.GetGeometryType(), order);
  }

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    IntegrationPoint eip1, eip2;

    Trans.Loc1.Transform(ip, eip1);
    Trans.Loc2.Transform(ip, eip2);

    trial_fe1.CalcShape(eip1, shape1);
    trial_fe2.CalcShape(eip2, shape2);
    test_fe1.CalcShape(eip1, test_shape1);
    test_fe2.CalcShape(eip2, test_shape2);

    Trans.Face->SetIntPoint(&ip);

    Vector nor(dim);
    CalcOrtho(Trans.Face->Jacobian(), nor);
    nor /= nor.Norml2(); // Unit normal

    double n_dir = nor(direction);
    double w = ip.weight * Trans.Face->Weight() * scaling;

    // Implement {{u}}·[[v]] term where {{u}} = (u1+u2)/2

    // Block (test1, trial1)
    for (int j = 0; j < test1_ndofs; j++) {
      for (int k = 0; k < trial1_ndofs; k++) {
        elmat(j, k) += 0.5 * w * n_dir * test_shape1(j) * shape1(k);
      }
    }

    // Block (test1, trial2)
    for (int j = 0; j < test1_ndofs; j++) {
      for (int k = 0; k < trial2_ndofs; k++) {
        elmat(j, trial1_ndofs + k) +=
            0.5 * w * n_dir * test_shape1(j) * shape2(k);
      }
    }

    // Block (test2, trial1)
    for (int j = 0; j < test2_ndofs; j++) {
      for (int k = 0; k < trial1_ndofs; k++) {
        elmat(test1_ndofs + j, k) -=
            0.5 * w * n_dir * test_shape2(j) * shape1(k);
      }
    }

    // Block (test2, trial2)
    for (int j = 0; j < test2_ndofs; j++) {
      for (int k = 0; k < trial2_ndofs; k++) {
        elmat(test1_ndofs + j, trial1_ndofs + k) -=
            0.5 * w * n_dir * test_shape2(j) * shape2(k);
      }
    }
  }
}

// UltraweakMixedDerivativeIntegrator implementation
UltraweakMixedDerivativeIntegrator::UltraweakMixedDerivativeIntegrator(
    Coefficient *q)
    : UltraweakDPGIntegrator(q) {}

void UltraweakMixedDerivativeIntegrator::AssembleElementMatrix(
    const FiniteElement &trial_fe, const FiniteElement &test_fe,
    ElementTransformation &Trans, DenseMatrix &elmat) {
  int trial_ndofs = trial_fe.GetDof();
  int test_ndofs = test_fe.GetDof();
  int dim = trial_fe.GetDim();

  MFEM_VERIFY(
      dim >= 2,
      "UltraweakMixedDerivativeIntegrator requires at least 2D elements");

  Vector shape(trial_ndofs);
  DenseMatrix dshape(test_ndofs, dim);
  DenseMatrix d2shape;
  Vector d2shape_mixed(test_ndofs);

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = trial_fe.GetOrder() + test_fe.GetOrder();
    ir = &IntRules.Get(trial_fe.GetGeomType(), order);
  }

  elmat.SetSize(test_ndofs, trial_ndofs);
  elmat = 0.0;

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    Trans.SetIntPoint(&ip);

    trial_fe.CalcShape(ip, shape);

    // Calculate mixed derivative d^2/dSdV (direction 0 and 1)
    d2shape.SetSize(test_ndofs, dim * dim);
    test_fe.CalcHessian(ip, d2shape);

    // Extract the mixed derivative (d^2/dSdV)
    for (int j = 0; j < test_ndofs; j++) {
      d2shape_mixed(j) = d2shape(j, 0 * dim + 1); // S=0, V=1
    }

    double coef = 1.0;
    if (Q) {
      coef = Q->Eval(Trans, ip);
    }

    double w = ip.weight * coef * Trans.Weight();

    for (int j = 0; j < test_ndofs; j++) {
      for (int k = 0; k < trial_ndofs; k++) {
        elmat(j, k) += d2shape_mixed(j) * shape(k) * w;
      }
    }
  }
}

// OptimalTestFunctionComputer implementation
OptimalTestFunctionComputer::OptimalTestFunctionComputer(
    FiniteElementSpace *tr_fes, FiniteElementSpace *te_fes,
    MixedBilinearForm *dpg, BilinearForm *ip)
    : trial_fes(tr_fes), test_fes(te_fes), dpg_form(dpg), test_ip(ip) {
  MFEM_VERIFY(trial_fes != nullptr, "Trial FiniteElementSpace is null");
  MFEM_VERIFY(test_fes != nullptr, "Test FiniteElementSpace is null");
  MFEM_VERIFY(dpg_form != nullptr, "MixedBilinearForm is null");
  MFEM_VERIFY(test_ip != nullptr, "BilinearForm is null");

  MFEM_VERIFY(trial_fes->GetMesh() == test_fes->GetMesh(),
              "Trial and test spaces must use the same mesh");
}

void OptimalTestFunctionComputer::ComputeOptimalTestFunction(
    int elem_idx, int trial_dof_idx, Vector &optimal_test) {
  // Check if the element index is valid
  Mesh *mesh = trial_fes->GetMesh();
  MFEM_VERIFY(mesh != nullptr, "Mesh is null");
  MFEM_VERIFY(elem_idx >= 0 && elem_idx < mesh->GetNE(),
              "Element index out of range");

  // Make sure forms are assembled
  dpg_form->Assemble();
  test_ip->Assemble();

  // Get local DOFs
  Array<int> trial_dofs, test_dofs;
  trial_fes->GetElementDofs(elem_idx, trial_dofs);
  test_fes->GetElementDofs(elem_idx, test_dofs);

  // Check trial DOF index
  MFEM_VERIFY(trial_dof_idx >= 0 && trial_dof_idx < trial_dofs.Size(),
              "Trial DOF index out of range");

  // Extract the local B matrix
  DenseMatrix B(test_dofs.Size(), trial_dofs.Size());
  B = 0.0;

  SparseMatrix &B_mat = dpg_form->SpMat();
  for (int i = 0; i < test_dofs.Size(); i++) {
    int row = test_dofs[i];
    if (row >= B_mat.Height())
      continue;

    for (int j = 0; j < trial_dofs.Size(); j++) {
      int col = trial_dofs[j];
      if (col >= B_mat.Width())
        continue;

      double val = B_mat.Elem(row, col);
      if (val != 0.0) {
        B(i, j) = val;
      }
    }
  }

  // Extract the local G matrix (test inner product)
  DenseMatrix G(test_dofs.Size(), test_dofs.Size());
  G = 0.0;

  SparseMatrix &G_mat = test_ip->SpMat();
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
        G(i, j) = val;
      }
    }
  }

  // Add regularization to ensure invertibility
  double min_diag = std::numeric_limits<double>::max();
  for (int i = 0; i < G.Height(); i++) {
    if (G(i, i) < min_diag && G(i, i) > 0) {
      min_diag = G(i, i);
    }

    if (G(i, i) <= 1e-14) {
      G(i, i) = 1e-10;
    }
  }

  // Create a unit vector for the trial basis function
  Vector e_j(trial_dofs.Size());
  e_j = 0.0;
  e_j(trial_dof_idx) = 1.0;

  // Compute B*e_j
  Vector B_ej(test_dofs.Size());
  B.Mult(e_j, B_ej);

  // Solve G*v_j = B*e_j for the optimal test function
  optimal_test.SetSize(test_dofs.Size());

  // Use direct solver for the local problem
  DenseMatrixInverse G_inv(G);
  try {
    G_inv.Mult(B_ej, optimal_test);
  } catch (const std::exception &e) {
    std::cerr << "Warning: Exception when computing optimal test function: "
              << e.what() << std::endl;

    // Use diagonal approximation as fallback
    optimal_test = 0.0;
    for (int i = 0; i < test_dofs.Size(); i++) {
      if (G(i, i) > 1e-14) {
        optimal_test(i) = B_ej(i) / G(i, i);
      }
    }
  }
}

void OptimalTestFunctionComputer::ComputeAllOptimalTestFunctions(
    int elem_idx, DenseMatrix &optimal_test_functions) {
  // Get local DOFs
  Array<int> trial_dofs, test_dofs;
  trial_fes->GetElementDofs(elem_idx, trial_dofs);
  test_fes->GetElementDofs(elem_idx, test_dofs);

  optimal_test_functions.SetSize(test_dofs.Size(), trial_dofs.Size());

  // Compute each column of the optimal test function matrix
  Vector optimal_test;
  for (int j = 0; j < trial_dofs.Size(); j++) {
    ComputeOptimalTestFunction(elem_idx, j, optimal_test);

    // Set column j of the matrix
    for (int i = 0; i < test_dofs.Size(); i++) {
      optimal_test_functions(i, j) = optimal_test(i);
    }
  }
}

// UltraweakDriftIntegrator implementation
UltraweakDriftIntegrator::UltraweakDriftIntegrator(Coefficient *q)
    : UltraweakDPGIntegrator(q) {}

void UltraweakDriftIntegrator::AssembleElementMatrix(
    const FiniteElement &trial_fe, const FiniteElement &test_fe,
    ElementTransformation &Trans, DenseMatrix &elmat) {
  int trial_ndofs = trial_fe.GetDof();
  int test_ndofs = test_fe.GetDof();

  Vector shape_trial(trial_ndofs);
  Vector shape_test(test_ndofs);

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = trial_fe.GetOrder() + test_fe.GetOrder();
    ir = &IntRules.Get(trial_fe.GetGeomType(), order);
  }

  elmat.SetSize(test_ndofs, trial_ndofs);
  elmat = 0.0;

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    Trans.SetIntPoint(&ip);

    trial_fe.CalcShape(ip, shape_trial);
    test_fe.CalcShape(ip, shape_test);

    double coef = 1.0;
    if (Q) {
      coef = Q->Eval(Trans, ip);
    }

    double w = ip.weight * coef * Trans.Weight();

    for (int j = 0; j < test_ndofs; j++) {
      for (int k = 0; k < trial_ndofs; k++) {
        elmat(j, k) += shape_test(j) * shape_trial(k) * w;
      }
    }
  }
}

// UltraweakReactionIntegrator implementation
UltraweakReactionIntegrator::UltraweakReactionIntegrator(Coefficient *q)
    : UltraweakDPGIntegrator(q) {}

void UltraweakReactionIntegrator::AssembleElementMatrix(
    const FiniteElement &trial_fe, const FiniteElement &test_fe,
    ElementTransformation &Trans, DenseMatrix &elmat) {
  int trial_ndofs = trial_fe.GetDof();
  int test_ndofs = test_fe.GetDof();

  Vector shape_trial(trial_ndofs);
  Vector shape_test(test_ndofs);

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = trial_fe.GetOrder() + test_fe.GetOrder();
    ir = &IntRules.Get(trial_fe.GetGeomType(), order);
  }

  elmat.SetSize(test_ndofs, trial_ndofs);
  elmat = 0.0;

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    Trans.SetIntPoint(&ip);

    trial_fe.CalcShape(ip, shape_trial);
    test_fe.CalcShape(ip, shape_test);

    double coef = 1.0;
    if (Q) {
      coef = Q->Eval(Trans, ip);
    }

    double w = ip.weight * coef * Trans.Weight();

    for (int j = 0; j < test_ndofs; j++) {
      for (int k = 0; k < trial_ndofs; k++) {
        elmat(j, k) += shape_test(j) * shape_trial(k) * w;
      }
    }
  }
}

// TestInnerProductIntegrator implementation
TestInnerProductIntegrator::TestInnerProductIntegrator(double diff_coef,
                                                       double mass_coef)
    : alpha_diff(diff_coef), alpha_mass(mass_coef) {}

void TestInnerProductIntegrator::AssembleElementMatrix(
    const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat) {
  int ndofs = el.GetDof();
  int dim = el.GetDim();

  Vector shape(ndofs);
  DenseMatrix dshape(ndofs, dim);

  const IntegrationRule *ir = IntRule;
  if (ir == nullptr) {
    int order = 2 * el.GetOrder();
    ir = &IntRules.Get(el.GetGeomType(), order);
  }

  elmat.SetSize(ndofs);
  elmat = 0.0;

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    Trans.SetIntPoint(&ip);

    el.CalcShape(ip, shape);
    el.CalcDShape(ip, dshape);

    double w = ip.weight * Trans.Weight();

    // Mass term: alpha_mass * (u,v)
    if (alpha_mass != 0.0) {
      for (int j = 0; j < ndofs; j++) {
        for (int k = 0; k < ndofs; k++) {
          elmat(j, k) += alpha_mass * shape(j) * shape(k) * w;
        }
      }
    }

    // Diffusion term: alpha_diff * (grad u, grad v)
    if (alpha_diff != 0.0) {
      for (int j = 0; j < ndofs; j++) {
        for (int k = 0; k < ndofs; k++) {
          double grad_dot = 0.0;
          for (int d = 0; d < dim; d++) {
            grad_dot += dshape(j, d) * dshape(k, d);
          }
          elmat(j, k) += alpha_diff * grad_dot * w;
        }
      }
    }
  }
} // namespace mfem