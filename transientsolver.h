#pragma once
#include "geometry.h"
#include "steadysolver.h"
#include <string>
#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <experimental/filesystem>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using SpMat = Eigen::SparseMatrix<double>;
using Trip = Eigen::Triplet<double>;

template<typename RefElement, typename Quadrature, unsigned Ng>
class TransientSolver
{
  static constexpr unsigned LocalDim = Ng * RefElement::Npts;
  static constexpr unsigned OffsetCorner = 0;
  static constexpr unsigned OffsetEdge = 4;
  static constexpr unsigned OffsetInternal = 4*(RefElement::NEdgeNodes+1);

  // Element mass matrix
  using LocalMat = Eigen::Matrix<double, LocalDim, LocalDim>;

  ReactorGeometry<RefElement>& geom;

  unsigned matdim; // number of rows and columns of whole-problem matrix
  unsigned n_precursor_groups;

  double dt;
  double tfinal;

  std::vector<int> zero_rhs_indices;

  SpMat lhs_matrix;
  Eigen::VectorXd rhs;
  Eigen::VectorXd soln;
  Eigen::VectorXd precursors;
  Eigen::VectorXd fiss_source;

  Eigen::BiCGSTAB<SpMat> solver;
  std::string output_directory;

  double integrated_fission_source;
  double initial_fission_source;

public:

  TransientSolver(ReactorGeometry<RefElement>& geom_, SteadySolver<RefElement, Quadrature, Ng>& sted,
      double dt_, double tfinal_, std::string output_directory_) : geom(geom_),
    matdim(geom.npoints() * Ng),
    n_precursor_groups(geom.n_precursors),
    lhs_matrix(matdim, matdim),
    rhs(matdim),
    soln(matdim),
    dt(dt_),
    tfinal(tfinal_),
    precursors(n_precursor_groups * geom.npoints()),
    output_directory(output_directory_),
    fiss_source(geom.npoints())
  {

    // Copy eigenvector as initial condition
    soln = sted.soln;

    // Go ahead and allocate the structure of the sparse matrix.
    // The values will change over the course of the transient, but
    // the sparsity pattern will not.

    std::vector<Trip> lhs_matrix_entries;
    lhs_matrix_entries.reserve(Ng*geom.npoints());

    if (geom.n_precursors == 0) {
      std::cerr << "Either no transient data was supplied in the input, or it was supplied incorrectly." << std::endl;
      exit(1);
    }

    for (auto& q: geom.quads) {
      LocalMat local_lhs = LocalMat::Zero();

      // Factor to adjust nu-sigma-f by (transient-specific)
      double nsf_factor = 0.0;
      for (int i=0; i<n_precursor_groups; ++i) {
        nsf_factor += q.xs.betas[i] * q.xs.lambdas[i] / (1.0 + dt*q.xs.lambdas[i]);
      }
      nsf_factor *= dt;
      nsf_factor += 1.0-q.xs.beta;
      nsf_factor /= geom.k;


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
              entry += (q.xs.sigma_r[group]+dt/q.xs.velocities[group]) * mass;
              local_lhs(Ng*test+group, Ng*trial+group) += Quadrature::weights[w] * entry;

              // Now loop over in-scattering to this group:
              for (int gprime=0; gprime<Ng; ++gprime) {
                local_lhs(Ng*test+group, Ng*test+gprime) -= Quadrature::weights[w] * mass * q.xs.sig_s[Ng*group+gprime];
                local_lhs(Ng*test+group, Ng*test+gprime) -= Quadrature::weights[w] * mass * q.xs.chi[group] * q.xs.nu_sig_f[gprime] * nsf_factor;
              }
            }
          }
        }
      }

      // Now multiply by the Jacobian of the transformation
      // matrix (i.e., the area of the quadrilateral)
      double A = q.area();
      local_lhs *= A;

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
              zero_rhs_indices.push_back(Ng*q.loc2glb(c)+g);
              local_lhs(Ng*c+g, Ng*c+g) = 1.0;
            }
          }

          // Now all corresponding edge indices must have the same operation done
          for (int e=OffsetEdge+n*RefElement::NEdgeNodes; e<OffsetEdge+(n+1)*RefElement::NEdgeNodes; ++e) {
            for (int g=0; g<Ng; ++g) {
              local_lhs.row(Ng*e+g).setZero();
              zero_rhs_indices.push_back(Ng*q.loc2glb(e)+g);
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
            }
          }
        }
      }
    }
    lhs_matrix.setFromTriplets(lhs_matrix_entries.begin(), lhs_matrix_entries.end());

    std::cout << "initializing precursor fields..." << std::endl;
    // The fission source to each precursor has to be integrated on its
    // own because beta may be an arbitrary function in space.
    Eigen::VectorXd this_prec(geom.npoints());
    for (int prec=0; prec<n_precursor_groups; ++prec) {
      SpMat lhs_prec_mat(geom.npoints(), geom.npoints());
      std::vector<Trip> lhs_entries;
      fiss_source = Eigen::VectorXd::Zero(geom.npoints()); 
      // Integrate product of fission source and beta
      for (auto& q: geom.quads) {
        double A = q.area();
        double this_beta = q.xs.betas[prec];
        for (int w=0; w<Quadrature::points.size(); ++w) {
          double x = Quadrature::points[w].x;
          double y = Quadrature::points[w].y;
          std::array<double, RefElement::Npts> basis = RefElement::evaluate_basis(x, y);
          for (int test=0; test<RefElement::Npts; ++test) {
            int test_glb_indx = q.loc2glb(test);
            for (int trial=0; trial<RefElement::Npts; ++trial) {
              int tri_glb_indx = q.loc2glb(trial);
              double mass = basis[test] * basis[trial];
              lhs_entries.emplace_back(test_glb_indx, tri_glb_indx, mass * q.xs.lambdas[prec] * A * Quadrature::weights[w]);
              for (int g=0; g<Ng; ++g) {
                fiss_source(test_glb_indx) += this_beta * q.xs.nu_sig_f[g] * soln(tri_glb_indx*Ng+g) * mass * A * Quadrature::weights[w] / geom.k;
              }
            }
          }
        }
      }

      // OK, now solve a simple linear system to get the precursor concentrations.
      // This corresponds to the L2 projection for spatially variable lambdas and betas.
      lhs_prec_mat.setFromTriplets(lhs_entries.begin(), lhs_entries.end());
      solver.compute(lhs_prec_mat);
      this_prec = solver.solve(fiss_source);
      std::cout << "Solved for precursor group " << prec << "with approximate error of " << solver.error() << std::endl;

      // Put temporary precursors to main array
      for (int i=0; i<geom.npoints(); ++i) {
        precursors(i*n_precursor_groups+prec) = this_prec(i);
      }
    }

    // Check if the output directory exists. If not, create it.
    std::experimental::filesystem::create_directory(output_directory);
  }

  // TODO
  void updateMatrix() {
  }

  void formRhs() {

    // This method also calculates the beginning-of-step fission source.
    integrated_fission_source = 0;

    // Have to integrate over whole domain, each time
    rhs = Eigen::VectorXd::Zero(rhs.size());
    for (auto& q: geom.quads) {
      double A = q.area();
      for (int w=0; w<Quadrature::points.size(); ++w) {
        double x = Quadrature::points[w].x;
        double y = Quadrature::points[w].y;
        std::array<double, RefElement::Npts> basis = RefElement::evaluate_basis(x, y);
        for (int test=0; test<RefElement::Npts; ++test) {
          int test_glb_indx = q.loc2glb(test);

          // Keep integrating fission source over the whole domain
          for (int g=0; g<Ng; ++g) integrated_fission_source += Quadrature::weights[w] * A * soln(test_glb_indx*Ng+g) * basis[test] * q.xs.nu_sig_f[g];

          for (int trial=0; trial<RefElement::Npts; ++trial) {
            int tri_glb_indx = q.loc2glb(trial);
            double mass = basis[test] * basis[trial];
            double prec_source = 0;
            for (int prec=0; prec<n_precursor_groups; ++prec) {
              prec_source += q.xs.lambdas[prec] * precursors[n_precursor_groups*tri_glb_indx+prec] / (1+q.xs.lambdas[prec]*dt);
            }
            for (int g=0; g<Ng; ++g) {
              rhs(Ng*test_glb_indx+g) += soln(tri_glb_indx*Ng+g) * mass * A / (q.xs.velocities[g]*dt) * Quadrature::weights[w];
              rhs(Ng*test_glb_indx+g) += q.xs.chi_d[g] * prec_source * mass * A * Quadrature::weights[w];
            }
          }
        }
      }
    }

    // Now zero out the dirichlet BC spots:
    for (int i: zero_rhs_indices) {
      rhs(i) = 0;
    }
  }

  void solveNextFlux() {
    solver.compute(lhs_matrix);
    soln = solver.solve(rhs);
    std::cout << "converged in " << solver.iterations() << std::endl;
    std::cout << "  estimated error: " << solver.error() << std::endl;
  }

  void updatePrecursors() {
  }

  void writeOutput(double t) {
    std::string filename = output_directory + "/" + std::to_string(t) + ".vtu";
    std::ofstream f;
    f.open(filename, std::fstream::out);
    geom.printUpToPointData(f);
    // Get min/max for each group:
    std::array<double, Ng> group_mins = {1e100};
    std::array<double, Ng> group_maxs = {0.0};
    std::vector<double> prec_mins(n_precursor_groups, 1e100);
    std::vector<double> prec_maxs(n_precursor_groups, 0.0);
    for (int pt=0; pt<geom.npoints(); ++pt) {
      for (int g=0; g<Ng; ++g) {
        if (soln(Ng*pt+g) > group_maxs[g]) group_maxs[g] = soln(Ng*pt+g);
        if (soln(Ng*pt+g) < group_mins[g]) group_mins[g] = soln(Ng*pt+g);
      }
      for (int prec=0; prec<n_precursor_groups; ++prec) {
        if (precursors(n_precursor_groups*pt+prec) > prec_maxs[prec]) prec_maxs[prec] = precursors(n_precursor_groups*pt+prec);
        if (precursors(n_precursor_groups*pt+prec) < prec_mins[prec]) prec_mins[prec] = precursors(n_precursor_groups*pt+prec);
      }
    }
    for (int g=0; g<Ng; ++g) {
      f << "    <DataArray type=\"Float32\" Name=\"" << "group" + std::to_string(g) << "\" format=\"ascii\" RangeMin=\"" << group_mins[g] << "\" RangeMax=\"" << group_maxs[g] << "\">" << std::endl;
      for (int pt=0; pt<geom.npoints(); ++pt) {
        f << soln(Ng*pt+g) << " ";
        if (pt % 10 == 0 and pt > 0) f << std::endl;
      }
      f << "    </DataArray>" << std::endl;
    }
    for (int prec=0; prec<n_precursor_groups; ++prec) {
      f << "    <DataArray type=\"Float32\" Name=\"" << "prec" + std::to_string(prec) << "\" format=\"ascii\" RangeMin=\"" << prec_mins[prec] << "\" RangeMax=\"" << prec_maxs[prec] << "\">" << std::endl;
      for (int pt=0; pt<geom.npoints(); ++pt) {
        f << precursors(n_precursor_groups*pt+prec) << " ";
        if (pt % 10 == 0 and pt > 0) f << std::endl;
      }
      f << "    </DataArray>" << std::endl;
    }
    geom.printAfterPointData(f);
  }

  // Executes the transient
  void run() {
    double time = 0;

    while (time < tfinal) {
      std::cout << "t=" << time << "s" << std::endl;
      writeOutput(time);
      updateMatrix();

      formRhs();
      if (time == 0) initial_fission_source = integrated_fission_source;
      std::cout << "P = " << integrated_fission_source / initial_fission_source << std::endl;

      solveNextFlux();
      updatePrecursors();
      time += dt;
    }
  }
};
