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
  std::string id;
  std::vector<double> diffcoeffs;
  std::vector<double> sigma_a;
  std::vector<double> nu_sig_f;
  std::vector<double> sig_s;
  std::vector<double> chi;

  std::vector<double> sigma_r; // removal cross section, calculated after-the-fact

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

  bool showmatrix {false};

  unsigned npoints;

  ParsedAssembly& lookupMaterial(std::string& id);
};
