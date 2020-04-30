#pragma once

// Parses input files
#include <exception>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

template <typename C>
bool all(C& array) {
  for (int i=0; i<array.size(); ++i) {
    if (not array[i]) return false;
  }
  return true;
}

class ParseError: public std::exception
{
public:
    virtual const char* what() const throw (){
       return "No keyword found after \"set\" keyword!";
    }
};
class InvalidKeyword : public std::exception
{
public:
    explicit InvalidKeyword(const std::string& message):
      word(message) {}
    virtual const char* what() const throw (){
       return word.c_str();
    }
protected:
    std::string word;
};

struct ParsedAssembly
{
  unsigned n_delayed_groups;
  double beta;

  std::string id;
  std::vector<double> diffcoeffs;
  std::vector<double> sigma_a;
  std::vector<double> nu_sig_f;
  std::vector<double> sig_s;
  std::vector<double> chi;

  std::vector<double> tot_scattering; // total scattering cross section, calculated after-the-fact

  std::vector<double> lambdas;
  std::vector<double> betas;
  std::vector<double> chi_d;
  std::vector<double> velocities;
  std::vector<double> gam; // doppler feedback coefficients

  // sums columns of scattering matrix
  void calculate_tot_scattering() {
    int ng = sigma_a.size();
    tot_scattering.resize(ng);

    for (int g=0; g<ng; ++g)
      tot_scattering[g] = 0;

    for (int row=0; row<ng; ++row) {
      for (int col=0; col<ng; ++col) {
        tot_scattering[col] += sig_s[row*ng+col];
      }
    }
  }

  unsigned checkTransientData() {
    // Checks that all transient data has been provided, if one piece was given
    // Returns the number of delayed neutron groups provided if all requisite data's there.
    unsigned size = lambdas.size();
    if (size == 0) return 0;
    if (betas.size() != size) return 0;
    if (chi_d.size() != diffcoeffs.size()) return 0;
    if (velocities.size() != diffcoeffs.size()) return 0;
    if (gam.size() != diffcoeffs.size()) return 0;

    // sum delayed neutron fractions
    beta = 0.0;
    for (auto b: betas) beta += b;

    return size;
  }

  void setDelayedGroups() {
    // Sets the delayed group count, only if all required data was provided
    n_delayed_groups = checkTransientData();
  }
};

struct MaterialPerturbationTable {
  std::string mat_id;
  ParsedAssembly* mat_ptr {nullptr}; // this gets set after all materials are read
  std::vector<double> times;
  std::vector<double> values;

  enum class XSType {
    DiffCoeff,
    SigA,
    NuSigF
  } xstype;

  int group;
  double orig_xs;

  // Save the original cross section value
  void setOriginal() {
    if (mat_ptr) {
      if (xstype == XSType::DiffCoeff) orig_xs = mat_ptr->diffcoeffs[group];
      else if (xstype == XSType::SigA) orig_xs = mat_ptr->sigma_a[group];
      else if (xstype == XSType::NuSigF) orig_xs = mat_ptr->nu_sig_f[group];
      else std::cerr << "WTF in pert table" << std::endl;
    } else
      std::cerr << "must set material pointer first" << std::endl;
  }

  // Sets the perturbed cross section at time t
  void setPert(double t) {
    // Lookup multiplication factor. Linear search cuz the table is usually very small
    double fac;
    int ti;
    for (ti=0; ti<times.size(); ++ti)
      if (times[ti] > t) { ti--; break; }
    if (ti == times.size()) fac = values[ti-1];
    else if (ti == -1) fac = values[0];
    else fac = (values[ti+1]-values[ti])/(times[ti+1]-times[ti])*(t-times[ti]) + values[ti];

    if (xstype == XSType::DiffCoeff) mat_ptr->diffcoeffs[group] = orig_xs * fac;
    else if (xstype == XSType::SigA) mat_ptr->sigma_a[group] = orig_xs * fac;
    else if (xstype == XSType::NuSigF) mat_ptr->nu_sig_f[group] = orig_xs * fac;
    else std::cerr << "WTF set pert" << std::endl;

    // std::cout << "fac = " << fac << std::endl;
  }

};

class ParsedInput
{
public:
  ParsedInput(std::string filename);

  enum class RunMode {
    eigenvalue,
    transient } runmode;

  int refine {1}; // mesh refinement level

  std::string title;

  std::vector<ParsedAssembly> assemblies;
  std::vector<std::vector<std::string>> core_geom;

  double axial_buckling;
  double assembly_size;

  double dt {1e-3};
  double tfinal {0.1};
  bool quiet {false};

  double eig_tol {1e-6}; // eigenproblem residual L2 norm to stop at
  double linsolve_tol {0}; // zero means use Eigen's default epsilon
  double P0 {1}; // initial core power in transient calculations
  double rhocp {1}; // initial core power in transient calculations
  double T0 {273}; // initial core temperature in transient calculations

  bool showmatrix {false};

  unsigned npoints;

  ParsedAssembly& lookupMaterial(std::string& id);

  // This checks if transient data is present for each assembly, and
  // also checks if the same number of delayed neutron group data for
  // each assembly is present. If any of those conditions are not met,
  // this returns false.
  bool hasTransientData();

  unsigned getPrecursorCount();

  std::vector<MaterialPerturbationTable> perturbations;
};
