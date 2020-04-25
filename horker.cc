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
#include "quadrature.h"
#include <string>
#include <iostream>

constexpr unsigned Degree = 5;
constexpr unsigned Ngroups = 2;

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

  // Create solver
  SteadySolver<ReferenceElement<Degree>, RabinowitzRichter6, Ngroups> solver(geom);

  if (input.showmatrix) {
    std::cout << "LHS" << std::endl;
    solver.showLhsMatrix();
    std::cout << "RHS" << std::endl;
    solver.showRhsMatrix();
    exit(0);
  }

  solver.inversePowerIteration();
  auto group1 = solver.getGroupFlux(0);
  auto group2 = solver.getGroupFlux(1);
  std::string name = "group1";
  std::string name2 = "group2";
  geom.printDataToVTK("test.vtu", {&group1, &group2}, {&name, &name2});

  // geom.printGeometryToVTK("hope.vtu");
}
