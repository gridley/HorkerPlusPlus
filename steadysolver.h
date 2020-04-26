#pragma once
#include "geometry.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>

using SpMat = Eigen::SparseMatrix<double>;
using Trip = Eigen::Triplet<double>;

// Ng = number of neutron groups
template<typename RefElement, typename Quadrature, unsigned Ng>
class SteadySolver
{
  static constexpr unsigned LocalDim = Ng * RefElement::Npts;

  // Shape function ordering offsets. These are the indices in the
  // basis array where the first corner, edge, and internal basis functions
  // occur. These are defined for convenience.
  static constexpr unsigned OffsetCorner = 0;
  static constexpr unsigned OffsetEdge = 4;
  static constexpr unsigned OffsetInternal = 4*(RefElement::NEdgeNodes+1);

  using LocalMat = Eigen::Matrix<double, LocalDim, LocalDim>;

  ReactorGeometry<RefElement>& geom;
  unsigned matdim; // number of rows and columns of whole-problem matrix

  SpMat lhs_matrix;
  SpMat rhs_matrix;

  // It's useful to pre-compute an LU of the LHS matrix.
  Eigen::SparseLU<SpMat> lhs_lu;


public:

  Eigen::VectorXd soln;
  Eigen::VectorXd fiss_source;

  // constructor takes a geometry specification object
  SteadySolver(ReactorGeometry<RefElement>& geom_) : geom(geom_),
    matdim(geom.points.size() * Ng),
    lhs_matrix(matdim, matdim),
    rhs_matrix(matdim, matdim),
    soln(matdim),
    fiss_source(matdim)
  {
    // Construct sparse matrices
    std::vector<Trip> lhs_matrix_entries;
    std::vector<Trip> rhs_matrix_entries;
    lhs_matrix_entries.reserve(Ng*geom.points.size());
    rhs_matrix_entries.reserve(Ng*geom.points.size());

    // Loop over elements, compute local mass/stiffness matrices,
    // and put them back to the lhs and rhs entries using the
    // local-to-global mapping.
    for (auto& q: geom.quads) {
      // Local and RHS matrix to calculate
      LocalMat local_lhs = LocalMat::Zero();
      LocalMat local_rhs = LocalMat::Zero();

      // TODO: for general deformed quadrilaterals, this approach doesn't work,
      // and the true jacobian must be used to scale the derivative terms
      double dx2 = geom.dx_nom * geom.dx_nom;

      // Loop over quadrature points:
      for (int w=0; w<Quadrature::points.size(); ++w) {
        // Loop over test functions:
        double x = Quadrature::points[w].x;
        double y = Quadrature::points[w].y;
        std::array<double, RefElement::Npts> basis = RefElement::evaluate_basis(x, y);
        std::array<double, RefElement::Npts> der_x = RefElement::evaluate_deriv_x(x, y);
        std::array<double, RefElement::Npts> der_y = RefElement::evaluate_deriv_y(x, y);

        for (int test=0; test<RefElement::Npts; ++test) {
          for (int trial=0; trial<RefElement::Npts; ++trial) {
            for (int group=0; group<Ng; ++group) {
              double entry = q.xs.diffcoeffs[group] * (der_x[test]*der_x[trial]+
                                                       der_y[test]*der_y[trial])/dx2;
              double mass = basis[test] * basis[trial];
              entry += q.xs.sigma_r[group] * mass;
              local_lhs(Ng*test+group, Ng*trial+group) += Quadrature::weights[w] * entry;

              // Now loop over in-scattering to this group:
              for (int gprime=0; gprime<Ng; ++gprime) {
                local_lhs(Ng*test+group, Ng*test+gprime) -= Quadrature::weights[w] * mass * q.xs.sig_s[Ng*group+gprime];
                local_rhs(Ng*test+group, Ng*test+gprime) += Quadrature::weights[w] * mass * q.xs.chi[group] * q.xs.nu_sig_f[gprime];
              }
            }
          }
        }
      }

      // Now multiply by the Jacobian of the transformation
      // matrix (i.e., the area of the quadrilateral)
      double A = q.area();
      local_lhs *= A;
      local_rhs *= A;

      // Edit in vacuum boundary conditions. Neumann boundary conditions should
      // come naturally, if not, test functions should be modified in some way?
      // Moreover, this has to be done for each group.
      for (int n=0; n<q.neighbors.size(); ++n) {
        auto& bdry = q.neighbors[n];
        if (bdry.type == ElementBoundary<RefElement>::BoundaryType::Vacuum) {
          std::array<int, 2> zero_corners = {n, (n+1)%4};
          for (int c: zero_corners) {
            for (int g=0; g<Ng; ++g) {
              local_lhs.row(Ng*c+g).setZero();
              local_rhs.row(Ng*c+g).setZero();
              local_lhs(Ng*c+g, Ng*c+g) = 1.0;
            }
          }

          // Now all corresponding edge indices must have the same operation done
          for (int e=OffsetEdge+n*RefElement::NEdgeNodes; e<OffsetEdge+(n+1)*RefElement::NEdgeNodes; ++e) {
            for (int g=0; g<Ng; ++g) {
              local_lhs.row(Ng*e+g).setZero();
              local_rhs.row(Ng*e+g).setZero();
              local_lhs(Ng*e+g, Ng*e+g) = 1.0;
            }
          }
        }
      }

      // Put local matrix to the global sparse matrix specification
      for (int test=0; test<RefElement::Npts; ++test) {
        for (int trial=0; trial<RefElement::Npts; ++trial) {
          for (int g=0; g<Ng; ++g) {
            for (int gprime=0; gprime<Ng; ++gprime) {
              lhs_matrix_entries.emplace_back(Ng*q.loc2glb(test)+g, Ng*q.loc2glb(trial)+gprime, local_lhs(Ng*test+g, Ng*trial+gprime));
              rhs_matrix_entries.emplace_back(Ng*q.loc2glb(test)+g, Ng*q.loc2glb(trial)+gprime, local_rhs(Ng*test+g, Ng*trial+gprime));
            }
          }
        }
      }
    }
    lhs_matrix.setFromTriplets(lhs_matrix_entries.begin(), lhs_matrix_entries.end());
    rhs_matrix.setFromTriplets(rhs_matrix_entries.begin(), rhs_matrix_entries.end());
  }

  void showLhsMatrix() {
    std::cout << Eigen::MatrixXd(lhs_matrix) << std::endl;
  }
  void showRhsMatrix() {
    std::cout << Eigen::MatrixXd(rhs_matrix) << std::endl;
  }

  void inversePowerIteration() {
    // Note, Eigen-lib sadly has no routines for doing eigenproblems with sparse matrices
    std::cout << "Computing LHS LU..." << std::endl;
    lhs_lu.analyzePattern(lhs_matrix);
    soln = Eigen::MatrixXd::Constant(matdim, 1, 1.0);

    lhs_lu.factorize(lhs_matrix);
    int error = lhs_lu.info();
    if (error != 0) {
      std::cerr << "Unable to LU factorize matrix. Error value: " << error << std::endl;
      exit(1);
    }

    // Do a couple power iterations
    double k=1.0;
    int max_outer = 10000;
    Eigen::VectorXd residual(soln.size());
    for (int outer=0; outer<max_outer; ++outer) {
      std::cout << k << std::endl;
      soln /= soln.norm();
      fiss_source = rhs_matrix * soln / k;
      soln = lhs_lu.solve(fiss_source);
      double norm2 = soln.norm();
      k *= norm2;

      // Every ten iterations, check the eigenproblem residual.
      if (outer%10==0 and outer>0) {
        residual = (lhs_matrix * soln - rhs_matrix * soln/k);
        double resid_norm = residual.norm();
        // std::cout << "eigenproblem residual norm: " << resid_norm << std::endl;
        if (resid_norm < geom.eig_tol) break;
      }
    }
    geom.k = k;
  }

  Eigen::VectorXd getGroupFlux(int g) {
    if (g >= Ng) { std::cerr << "asking for too high of a group flux!!!" << std::endl; exit(1); }
    Eigen::VectorXd result(matdim/Ng);
    int k=0;
    for (int i=g; i<soln.size(); i+=Ng) {
      result(k) = soln(i);
      k++;
    }
    return result;
  }
};
