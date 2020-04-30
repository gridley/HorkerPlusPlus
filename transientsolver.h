#pragma once
#include "geometry.h"
#include "steadysolver.h"
#include "input_parsing.h" // for the material perturbation list
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
  double sqrtT0;

  std::vector<int> zero_rhs_indices;

  SpMat lhs_matrix;
  SpMat precursor_update_matrix;
  SpMat temperature_matrix;
  SpMat mass_matrix;
  Eigen::VectorXd rhs;
  Eigen::VectorXd precursor_rhs;
  Eigen::VectorXd soln;
  Eigen::VectorXd precursors;
  Eigen::VectorXd fiss_source;
  Eigen::VectorXd temperature;
  Eigen::VectorXd temperature_rhs;
  Eigen::VectorXd power;

  Eigen::BiCGSTAB<SpMat> solver;
  std::string output_directory;

  double integrated_fission_source;
  double initial_fission_source;
  double power_scale_factor {1};

  // Stores cross section perturbation data that initiate transients
  std::vector<MaterialPerturbationTable>& perts;
public:

  TransientSolver(ReactorGeometry<RefElement>& geom_, SteadySolver<RefElement, Quadrature, Ng>& sted,
      double dt_, double tfinal_, std::string output_directory_, std::vector<MaterialPerturbationTable>& perts_a, double tol=0) : geom(geom_),
    matdim(geom.npoints() * Ng),
    n_precursor_groups(geom.n_precursors),
    lhs_matrix(matdim, matdim),
    temperature_matrix(geom.npoints(), geom.npoints()),
    rhs(matdim),
    soln(matdim),
    dt(dt_),
    tfinal(tfinal_),
    precursors(n_precursor_groups * geom.npoints()),
    output_directory(output_directory_),
    fiss_source(geom.npoints()),
    precursor_update_matrix(n_precursor_groups*geom.npoints(), n_precursor_groups*geom.npoints()),
    precursor_rhs(n_precursor_groups*geom.npoints()),
    perts(perts_a),
    temperature(geom.npoints()),
    temperature_rhs(geom.npoints()),
    mass_matrix(geom.npoints(), geom.npoints()),
    power(geom.npoints())
  {
    // set initial core temperature
    temperature = Eigen::VectorXd::Ones(geom.npoints()) * geom.input_.T0;

    // set reference temperature for doppler effect
    sqrtT0 = std::sqrt(geom.input_.T0);

    // Possibly set linear solver tolerance
    if (tol != 0) {
      solver.setTolerance(tol);
    }

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
      double dx2 = geom.dx* geom.dx;

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
            int trial_glb = q.loc2glb(trial);
            for (int group=0; group<Ng; ++group) {
              double entry = q.xs.diffcoeffs[group] * (der_x[test]*der_x[trial]+
                                                       der_y[test]*der_y[trial])/dx2;
              double mass = basis[test] * basis[trial];
              double T = temperature[trial_glb]>geom.input_.T0 ? temperature[trial_glb] : geom.input_.T0;
              entry += (q.xs.sigma_a[group]*(1+q.xs.gam[group]*(std::sqrt(T)-sqrtT0))
                  +q.xs.tot_scattering[group]+1/(dt*q.xs.velocities[group])) * mass;
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

    if (not geom.input_.quiet)
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
      if (not geom.input_.quiet)
        std::cout << "Solved for precursor group " << prec << "with approximate error of " << solver.error() << std::endl;

      // Put temporary precursors to main array
      for (int i=0; i<geom.npoints(); ++i) {
        precursors(i*n_precursor_groups+prec) = this_prec(i);
      }
    }

    // Check if the output directory exists. If not, create it.
    std::experimental::filesystem::create_directory(output_directory);

    // Now create the precursor update matrix, which won't be changing in
    // time as long as the lambda values aren't also changing in time.
    std::vector<Trip> prec_mat_values;
    for (auto& q: geom.quads) {
      for (int prec=0; prec<n_precursor_groups; ++prec) {
        double A = q.area();
        double this_lambda = q.xs.lambdas[prec];
        for (int w=0; w<Quadrature::points.size(); ++w) {
          double x = Quadrature::points[w].x;
          double y = Quadrature::points[w].y;
          std::array<double, RefElement::Npts> basis = RefElement::evaluate_basis(x, y);
          for (int test=0; test<RefElement::Npts; ++test) {
            int test_glb_indx = q.loc2glb(test);
            for (int trial=0; trial<RefElement::Npts; ++trial) {
              int tri_glb_indx = q.loc2glb(trial);
              double mass = basis[test] * basis[trial];
              prec_mat_values.emplace_back(n_precursor_groups*test_glb_indx+prec, n_precursor_groups*tri_glb_indx+prec, mass * (q.xs.lambdas[prec]+1/dt) * A * Quadrature::weights[w]);
            }
          }
        }
      }
    }
    precursor_update_matrix.setFromTriplets(prec_mat_values.begin(), prec_mat_values.end());

    // Now the temperature update mass matrix can be calculated
    std::vector<Trip> temperature_matrix_entries;
    std::vector<Trip> mass_matrix_entries;
    temperature_matrix_entries.reserve(geom.npoints() * 4);
    mass_matrix_entries.reserve(geom.npoints() * 4);
    for (auto& q: geom.quads) {
      double A = q.area();
      for (int w=0; w<Quadrature::points.size(); ++w) {
        double x = Quadrature::points[w].x;
        double y = Quadrature::points[w].y;
        std::array<double, RefElement::Npts> basis = RefElement::evaluate_basis(x, y);
        for (int test=0; test<RefElement::Npts; ++test) {
          int test_glb_indx = q.loc2glb(test);
          for (int trial=0; trial<RefElement::Npts; ++trial) {
            int tri_glb_indx = q.loc2glb(trial);
            double mass = basis[test] * basis[trial];
            // Note: could add spatially dependent rho cp in the future
            temperature_matrix_entries.emplace_back(test_glb_indx, tri_glb_indx, mass * geom.input_.rhocp * A * Quadrature::weights[w]);
            mass_matrix_entries.emplace_back(test_glb_indx, tri_glb_indx, mass * A * Quadrature::weights[w]);
          }
        }
      }
    }
    temperature_matrix.setFromTriplets(temperature_matrix_entries.begin(), temperature_matrix_entries.end());
    mass_matrix.setFromTriplets(mass_matrix_entries.begin(), mass_matrix_entries.end());
  }

  void updateMatrix() {
    // First zero-out the old matrix entries
    for (int k=0; k<lhs_matrix.outerSize(); ++k)
      for (SpMat::InnerIterator it(lhs_matrix, k); it; ++it)
        it.valueRef() = 0;

    // Then re-create the LHS matrix, just as in problem initialization
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
      double dx2 = geom.dx* geom.dx;

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
            int trial_glb = q.loc2glb(trial);
            for (int group=0; group<Ng; ++group) {
              double entry = q.xs.diffcoeffs[group] * (der_x[test]*der_x[trial]+
                                                       der_y[test]*der_y[trial])/dx2;
              double mass = basis[test] * basis[trial];
              double T = temperature[trial_glb]>geom.input_.T0 ? temperature[trial_glb] : geom.input_.T0;
              entry += (q.xs.sigma_a[group]*(1+q.xs.gam[group]*(std::sqrt(T)-sqrtT0))
                  +q.xs.tot_scattering[group]+1/(dt*q.xs.velocities[group])) * mass;
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
              local_lhs(Ng*c+g, Ng*c+g) = 1.0;
            }
          }

          // Now all corresponding edge indices must have the same operation done
          for (int e=OffsetEdge+n*RefElement::NEdgeNodes; e<OffsetEdge+(n+1)*RefElement::NEdgeNodes; ++e) {
            for (int g=0; g<Ng; ++g) {
              local_lhs.row(Ng*e+g).setZero();
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
              lhs_matrix.coeffRef(Ng*q.loc2glb(test)+g, Ng*q.loc2glb(trial)+gprime) += local_lhs(Ng*test+g, Ng*trial+gprime);
            }
          }
        }
      }
    }
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
    soln = solver.solveWithGuess(rhs, soln);
    if (not geom.input_.quiet) {
      std::cout << "  converged in " << solver.iterations() << std::endl;
      std::cout << "  estimated error: " << solver.error() << std::endl;
    }
  }

  void updatePrecursors() {
    // Compute RHS for precursor update equation consisting of fission source and previous timestep values
    precursor_rhs = Eigen::VectorXd::Zero(precursor_rhs.size());
    for (auto& q: geom.quads) {
      for (int prec=0; prec<n_precursor_groups; ++prec) {
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
              double mass = basis[test] * basis[trial] * A * Quadrature::weights[w];
              for (int g=0; g<Ng; ++g) {
                precursor_rhs(test_glb_indx*n_precursor_groups+prec) += this_beta * mass * q.xs.nu_sig_f[g]/geom.k * soln(tri_glb_indx*Ng+g);
              }
              precursor_rhs(test_glb_indx*n_precursor_groups+prec) += mass * precursors(n_precursor_groups*tri_glb_indx+prec) / dt;
            }
          }
        }
      }
    }
    solver.compute(precursor_update_matrix);
    precursors = solver.solveWithGuess(precursor_rhs, precursors);
    // std::cout << "    precursors solved in " << solver.iterations() << " with error " << solver.error() << std::endl;;
  }

  void writeOutput(int ti) {
    std::string filename = output_directory + "/"+ geom.input_.title+ std::to_string(ti) + ".vtu";
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

    // print out temperature
    f << "    <DataArray type=\"Float32\" Name=\"temperature\" format=\"ascii\" RangeMin=\"" << temperature.minCoeff() << "\" RangeMax=\"" << temperature.maxCoeff() << "\">" << std::endl;
    for (int pt=0; pt<geom.npoints(); ++pt) {
      f << temperature(pt) << " ";
      if (pt % 10 == 0 and pt > 0) f << std::endl;
    }
    f << "    </DataArray>" << std::endl;

    // print out power
    f << "    <DataArray type=\"Float32\" Name=\"power\" format=\"ascii\" RangeMin=\"" << power.minCoeff() << "\" RangeMax=\"" << power.maxCoeff() << "\">" << std::endl;
    for (int pt=0; pt<geom.npoints(); ++pt) {
      f << power(pt) << " ";
      if (pt % 10 == 0 and pt > 0) f << std::endl;
    }
    f << "    </DataArray>" << std::endl;

    geom.printAfterPointData(f);
  }

  void computeFissSource() {
    // Computes the fission rate at each node
    fiss_source = Eigen::VectorXd::Zero(fiss_source.size());
    for (auto& q: geom.quads) {
      double A = q.area();
      for (int w=0; w<Quadrature::points.size(); ++w) {
        double x = Quadrature::points[w].x;
        double y = Quadrature::points[w].y;
        std::array<double, RefElement::Npts> basis = RefElement::evaluate_basis(x, y);
        for (int test=0; test<RefElement::Npts; ++test) {
          int test_glb_indx = q.loc2glb(test);
          for (int trial=0; trial<RefElement::Npts; ++trial) {
            int tri_glb_indx = q.loc2glb(trial);
            double mass = basis[test] * basis[trial];
            for (int g=0; g<Ng; ++g) {
              fiss_source(test_glb_indx) += q.xs.nu_sig_f[g] * soln(tri_glb_indx*Ng+g) * mass * A * Quadrature::weights[w] / geom.k * power_scale_factor;
            }
          }
        }
      }
    }

    // Solve for power
    solver.compute(mass_matrix);
    power = solver.solve(fiss_source);
  }

  void computeNewTemperature() {
    temperature_rhs = Eigen::VectorXd::Zero(temperature_rhs.size());
    for (auto& q: geom.quads) {
      double A = q.area();
      for (int w=0; w<Quadrature::points.size(); ++w) {
        double x = Quadrature::points[w].x;
        double y = Quadrature::points[w].y;
        std::array<double, RefElement::Npts> basis = RefElement::evaluate_basis(x, y);
        for (int test=0; test<RefElement::Npts; ++test) {
          int test_glb_indx = q.loc2glb(test);
          for (int trial=0; trial<RefElement::Npts; ++trial) {
            int tri_glb_indx = q.loc2glb(trial);
            double mass = basis[test] * basis[trial];
            temperature_rhs(test_glb_indx) += temperature(tri_glb_indx) * mass * A * Quadrature::weights[w] * geom.input_.rhocp;
          }
        }
      }
    }
    temperature_rhs += dt * fiss_source; // scaled power
    solver.compute(temperature_matrix);
    temperature = solver.solve(temperature_rhs);
  }

  void calculatePowerScaleFactor() {
    // Calculates a power scaling factor so that the average initial core power matches that
    // given in the input file

    // First, the core volume (not reflector) must be calculated:
    double vol_core = 0;
    for (auto& q: geom.quads) {
      if (q.xs.nu_sig_f[0] > 0) vol_core += q.area();
    }
    power_scale_factor = geom.input_.P0 * vol_core / integrated_fission_source;
  }

  // Executes the transient
  void run() {
    std::ofstream power_file;
    power_file.open(output_directory+"/power_history");
    double time = 0;
    int ti = 0;
    while (time < tfinal) {

      if (not geom.input_.quiet) std::cout << "t=" << time << "s ";

      // Lookup material changes
      for (auto& pert: perts) pert.setPert(time);

      updateMatrix();

      formRhs();

      if (ti == 0) calculatePowerScaleFactor();
      computeFissSource();
      computeNewTemperature();

      if (time == 0) initial_fission_source = integrated_fission_source;
      if (not geom.input_.quiet) std::cout << "P = " << integrated_fission_source / initial_fission_source << std::endl;
      power_file << time << " " << integrated_fission_source/initial_fission_source << std::endl;

      if (ti%10==0) {
        writeOutput(ti);
        if (not geom.input_.quiet) std::cout << "    writing VTU..." << std::endl;
      }

      solveNextFlux();
      updatePrecursors();

      time += dt;
      ti++;
    }
  }
};
