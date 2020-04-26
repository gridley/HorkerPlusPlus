#include "input_parsing.h"
#include <iostream>
#include <sstream>
#include <array>

// TODO add warning if Horker++ is not compiled to work with the given
// amount of groups.

ParsedAssembly& ParsedInput::lookupMaterial(std::string& id) {
  for (auto& assem: assemblies) if (assem.id == id) return assem;

  std::cerr << "Unable to find assembly with material ID " << id << std::endl;
  exit(1);
  return assemblies[0];
}

ParsedInput::ParsedInput(std::string filename) {
  std::ifstream infile;
  infile.open(filename, std::ios::in);
  if (not infile.good()) {
    std::cerr << "unable to open " << filename << std::endl;
    exit(1);
  }

  std::string thisword;
  std::string nextword;
  while (infile >> thisword) {
    if (thisword == "set") {
      if (not (infile >> thisword)) throw ParseError();
      if (thisword == "title") {
        // Get full title
        std::getline(std::getline(infile, nextword, '"'), thisword, '"');
        title = thisword;
      } else if (thisword == "mode") {
        if (not (infile >> thisword)) throw ParseError();
        if (thisword == "eigenvalue") {
          runmode = RunMode::eigenvalue;
        } else if (thisword == "transient") {
          runmode = RunMode::transient;
        } else {
          std::cerr << "unrecognized run mode!!" << std::endl;
          exit(1);
        }
      } else if (thisword == "axial_buckling") {
        if (not (infile >> thisword)) throw ParseError();
        axial_buckling = std::stod(thisword);
      } else if (thisword == "assembly_size") {
        if (not (infile >> thisword)) throw ParseError();
        assembly_size = std::stod(thisword);
      } else if (thisword == "showmatrix") {
        showmatrix = true;
      } else if (thisword == "eig_tol") {
        if (not (infile >> thisword)) throw ParseError();
        eig_tol = std::stod(thisword);
      } else {
        throw InvalidKeyword(thisword);
      }
    } else if (thisword.at(0) == '%') {
      // comment character
      std::getline(infile, thisword);
    } else if (thisword == "assem") {
      ParsedAssembly this_assem;
      if (not (infile >> thisword)) throw ParseError();
      this_assem.id = thisword;

      std::array<bool, 5> found_cross_sections = {false};
      std::array<bool, 4> transient_data = {false};

      while (not all(found_cross_sections) or not (all(transient_data) or runmode==RunMode::eigenvalue)) {
        if (not (infile >> thisword)) throw ParseError();
        std::getline(infile, nextword);
        std::stringstream substring(nextword);

        if (thisword == "d") {
          while(substring >> thisword) this_assem.diffcoeffs.push_back(std::stod(thisword));
          found_cross_sections[0] = true;
        } else if (thisword == "sig_a") {
          while(substring >> thisword) this_assem.sigma_a.push_back(std::stod(thisword));
          found_cross_sections[1] = true;
        } else if (thisword == "nu_sig_f") {
          while(substring >> thisword) this_assem.nu_sig_f.push_back(std::stod(thisword));
          found_cross_sections[2] = true;
        } else if (thisword == "sig_s") {
          while(substring >> thisword) this_assem.sig_s.push_back(std::stod(thisword));
          found_cross_sections[3] = true;
        } else if (thisword == "chi") {
          while(substring >> thisword) this_assem.chi.push_back(std::stod(thisword));
          found_cross_sections[4] = true;
        } else if (thisword == "lambdas") {
          while(substring >> thisword) this_assem.lambdas.push_back(std::stod(thisword));
          transient_data[0] = true;
        } else if (thisword == "betas") {
          while(substring >> thisword) this_assem.betas.push_back(std::stod(thisword));
          transient_data[1] = true;
        } else if (thisword == "chi_d") {
          while(substring >> thisword) this_assem.chi_d.push_back(std::stod(thisword));
          transient_data[2] = true;
        } else if (thisword == "velocities") {
          while(substring >> thisword) this_assem.velocities.push_back(std::stod(thisword));
          transient_data[3] = true;
        } else throw InvalidKeyword(thisword);
      }

      this_assem.calculate_removal_xs();
      this_assem.setDelayedGroups();
      assemblies.push_back(this_assem);

    } else if (thisword == "corelayout") {
      if (not (infile >> thisword)) throw ParseError();
      int nrows = std::stoi(thisword);
      if (not (infile >> thisword)) throw ParseError();
      int ncols = std::stoi(thisword);
      core_geom.resize(nrows);
      for (auto& row: core_geom) row.resize(ncols);
      int row = nrows;
      while (row --> 0) {
        for (int col=0; col<ncols; ++col) {
          if (not (infile >> thisword)) throw ParseError();
          core_geom[row][col] = thisword;
        }
      }
    }
  }
  infile.close();

  // Now check that all assemblies in the core layout have specified XS
  for (auto& row: core_geom) {
    for (std::string& id: row) {
      if (id == "r" or id == "v") continue; // skip reflective and vacuum BC blocks
      bool found = false;
      for (auto& assem: assemblies) {
        if (assem.id == id) {
          found = true;
          break;
        }
      }
      if (not found) {
        std::cerr << "assembly from core geometry called " << id << " has no XS." << std::endl;
        exit(1);
      }
    }
  }

  if (assemblies.size() == 0) {
    std::cerr << "no assemblies!!!" << std::endl;
    exit(1);
  }

  // Check that all assemblies have the same number of groups.
  unsigned ng = assemblies[0].diffcoeffs.size();
  for (auto& assem: assemblies) {
    if (ng != assem.diffcoeffs.size() or
        ng != assem.sigma_a.size() or
        ng != assem.nu_sig_f.size() or
        ng*ng != assem.sig_s.size()) {
      std::cerr << "all assemblies must have the same number of groups" << std::endl;
    }
  }

  // Check that whole outer boundary is of BC type, not assembly
  bool good = true;
  int ncol = core_geom[0].size();
  std::string name;
  for (int i=1; i<core_geom.size()-1; ++i) {
    name = core_geom[i][0];
    if (name != "r" and name != "v") good = false;
    name = core_geom[i][ncol-1];
    if (name != "r" and name != "v") good = false;
  }
  for (int i=0; i<ncol; ++i) {
    name = core_geom[0][i];
    if (name != "r" and name != "v") good = false;
    name = core_geom[core_geom.size()-1][i];
    if (name != "r" and name != "v") good = false;
  }
  if (not good) {
    std::cerr << "geometry must be bordered only by r's and v's to define the boundary condition"
      << std::endl;
    exit(1);
  }

  // Add axial buckling to all assembly total cross sections
  for (auto& assem : assemblies) {
    for (int i=0; i<assem.diffcoeffs.size(); ++i) {
      assem.sigma_r[i] += axial_buckling * assem.diffcoeffs[i];
    }
  }

  // Now the number of unique points in the geometry can be calculated.
  npoints = ncol * core_geom[0].size();
}

bool ParsedInput::hasTransientData() {
  bool result = true;
  for (auto& assem: assemblies) {
    result = result || assem.n_delayed_groups;
    if (assem.n_delayed_groups == 0) return false;
  }
  return result;
}

unsigned ParsedInput::getPrecursorCount() {
  if (hasTransientData()) {
    return assemblies[0].n_delayed_groups;
  } else return 0;
}
