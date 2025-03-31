#ifndef HESTON_SLV_MULTIGRID_HPP
#define HESTON_SLV_MULTIGRID_HPP

#include "mfem.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace mfem {

// Mixed derivative integrator for the Heston PDE
class MixedCrossDerivativeIntegrator : public BilinearFormIntegrator {
private:
  Coefficient &Q;

public:
  MixedCrossDerivativeIntegrator(Coefficient &q);

  virtual void AssembleElementMatrix(const FiniteElement &el,
                                     ElementTransformation &Trans,
                                     DenseMatrix &elmat);
};

// Multigrid solver for the Heston SLV PDE
class HestonSLVMultigrid : public GeometricMultigrid {
private:
  // Heston model parameters
  double r, q, kappa, theta, xi, rho;

  // Storage for bilinear forms at each level
  Array<BilinearForm *> bfs;

public:
  HestonSLVMultigrid(FiniteElementSpaceHierarchy &fespaces, Array<int> &ess_bdr,
                     double r_, double q_, double kappa_, double theta_,
                     double xi_, double rho_);

  virtual ~HestonSLVMultigrid();

  void ConstructBilinearForm(FiniteElementSpace &fespace);
  void ConstructCoarseOperatorAndSolver(FiniteElementSpace &coarse_fespace);
  void ConstructOperatorAndSmoother(FiniteElementSpace &fespace, int level);
};

// Block multigrid preconditioner for the DPG system
class DPGSystemMultigrid : public Solver {
private:
  BlockOperator *A;
  Array<int> offsets;
  int num_vars;
  mutable BlockVector r, z;
  std::vector<std::unique_ptr<HestonSLVMultigrid>> mg_solvers;

public:
  DPGSystemMultigrid(BlockOperator *block_operator, Array<int> &block_offsets,
                     Array<Mesh *> &meshes, Array<double> &heston_params,
                     int max_levels = 3);

  virtual ~DPGSystemMultigrid();

  virtual void Mult(const Vector &x, Vector &y) const;
  void SetPrintLevel(int print_lvl);
};

// Helper function to create FE space hierarchy
FiniteElementSpaceHierarchy *
CreateHestonHPHierarchy(Mesh *mesh, int max_order, int h_levels, int p_levels);

} // namespace mfem

#endif // HESTON_SLV_MULTIGRID_HPP