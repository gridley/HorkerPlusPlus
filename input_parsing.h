#pragma once

// Parses input files
#include <exception>
#include <fstream>
#include <string>
#include <vector>

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

  std::vector<double> sigma_r; // removal cross section, calculated after-the-fact

  std::vector<double> lambdas;
  std::vector<double> betas;
  std::vector<double> chi_d;
  std::vector<double> velocities;

  // Converts the absorption XS to a removal XS
  void calculate_removal_xs() {
    int ng = sigma_a.size();
    sigma_r.resize(ng);
    std::vector<double> tot_scattering(ng, 0.0);
    for (int row=0; row<ng; ++row) {
      for (int col=0; col<ng; ++col) {
        tot_scattering[col] += sig_s[row*ng+col];
      }
    }
    for (int g=0; g<ng; ++g)
      sigma_r[g] = tot_scattering[g] + sigma_a[g];
  }

  unsigned checkTransientData() {
    // Checks that all transient data has been provided, if one piece was given
    // Returns the number of delayed neutron groups provided if all requisite data's there.
    unsigned size = lambdas.size();
    if (size == 0) return 0;
    if (betas.size() != size) return 0;
    if (chi_d.size() != diffcoeffs.size()) return 0;
    if (velocities.size() != diffcoeffs.size()) return 0;

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

class ParsedInput
{
public:
  ParsedInput(std::string filename);

  enum class RunMode {
    eigenvalue,
    transient } runmode;

  std::string title;

  std::vector<ParsedAssembly> assemblies;
  std::vector<std::vector<std::string>> core_geom;

  double axial_buckling;
  double assembly_size;

  double eig_tol {1e-6}; // eigenproblem residual L2 norm to stop at

  bool showmatrix {false};

  unsigned npoints;

  ParsedAssembly& lookupMaterial(std::string& id);

  // This checks if transient data is present for each assembly, and
  // also checks if the same number of delayed neutron group data for
  // each assembly is present. If any of those conditions are not met,
  // this returns false.
  bool hasTransientData();

  unsigned getPrecursorCount();
};
