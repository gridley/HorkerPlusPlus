// HORKER++
//
// High
// Order
// Reactor
// Kinetics and
// Evaluation of
// Reactivity
//
// Written by Gavin Ridley for an MIT 22.213 final project, solving the
// LRA transient benchmark with a new method. LRA, as detailed in:
//
//   Benchmark Problem Book, supplement 2, benchmark problem 14, Tech.
//   Report ANL-7416, Argonne National Laboratory, 1977.
//
// The code here is capable of generating code for arbitrarily high order bases,
// on quadrilaterals, and using them for reactor transient analysis. 
// The limiting factor is the quadrilateral cubature, which goes up
// to 50th order in the script generate_quadrature.py. This means that
// approximately an approximately 25th order polynomial basis can be
// employed, at most.
//
// A feature not commonly found in nodal codes is that the assembly mesh
// can be moved dynamically to simulate effects like core flowering.

#include "input_parsing.h"
#include "geometry.h"
#include "openmc/position.h"
#include "reference_elements.h"
#include "steadysolver.h"
#include "transientsolver.h"
#include "quadrature.h"
#include <string>
#include <iostream>

constexpr unsigned Degree = 4;
constexpr unsigned Ngroups = 2;
typedef Sommariva16 Quad;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Must provide _one_ command line argument, the textual input file."
              << std::endl;
    exit(1);
  }

  // Read input file
  std::string filename = argv[1];
  ParsedInput input(filename);
  std::cout << "--- " << input.title << " ---" << std::endl;

  // Initialize core geometry
  ReactorGeometry<ReferenceElement<Degree>> geom(input);

  // Create steady-state solver
  // This is always done for both steady and transient solves right
  // now. In the future, I could have output files for the steady
  // solve for transient solves to read in, but that's a bit much
  // effort for now.
  SteadySolver<ReferenceElement<Degree>, Quad, Ngroups> solver(geom);

  if (input.showmatrix) {
    std::cout << "LHS" << std::endl;
    solver.showLhsMatrix();
    std::cout << "RHS" << std::endl;
    solver.showRhsMatrix();
    exit(0);
  }

  solver.inversePowerIteration();

  // Prepare output data
  std::vector<decltype(solver.getGroupFlux(0))> group_vectors(Ngroups);
  std::vector<std::string> group_names;
  for (int g=0; g<Ngroups; ++g) {
    group_vectors[g] = solver.getGroupFlux(g);
    group_names.push_back("group" + std::to_string(g));
  }
  std::vector<decltype(group_vectors)::value_type*> group_pointers(group_vectors.size());
  std::vector<std::string*> group_name_ptr(group_vectors.size());
  for (int g=0; g<Ngroups; ++g) {
    group_pointers[g] = &group_vectors[g];
    group_name_ptr[g] = &group_names[g];
  }

  // Output steady-state results
  geom.printDataToVTK("steady_result.vtu", group_pointers, group_name_ptr);

  // Now prepare to run a transient, if that mode was set in input
  if (input.runmode == ParsedInput::RunMode::transient) {
    std::cout << "running in transient mode..." << std::endl;
    TransientSolver<ReferenceElement<Degree>, Quad, Ngroups> trans_solver(geom, solver, input.dt,
        input.tfinal, "transient_results", input.perturbations, input.linsolve_tol);
    trans_solver.run();
  }
}
