#pragma once
// Responsible for storing geometry data, and how the finite elements nodes
// are globally enumerated. Currently, only quadrilaterals are handled.
#include <vector>
#include <array>
#include <exception>
#include <utility>
#include <string>
#include <fstream>
#include <iostream>
#include "openmc/position.h"
#include "input_parsing.h"
#include <Eigen/Dense>

template<typename RefElement>
class Quadrilateral;

template<typename RefElement>
class ElementBoundary
{
  // If it does border another element, store a pointer to the neighbor
  Quadrilateral<RefElement>* neighbor;

public:

  enum class BoundaryType {
    Vacuum,
    Reflective,
    OtherElement
  } type;

  ElementBoundary() : type(BoundaryType::Vacuum) {}

  ElementBoundary(BoundaryType type_) : type(type_) {
    if (type == BoundaryType::OtherElement) {
      std::cerr << "Constructing an ElementBoundary with a neighbor, but no neighbor given!" << std::endl;
      exit(1);
    }
  }

  ElementBoundary(Quadrilateral<RefElement>* q) : type(BoundaryType::OtherElement), neighbor(q) {
    if (q == nullptr) std::cerr << "warning: passed null pointer as a neighbor" << std::endl;
  }

  Quadrilateral<RefElement>* getNeighbor() {
    if (type == BoundaryType::OtherElement) {
        return neighbor;
    } else {
      return nullptr;
    }
  }

  bool isVac() { return type == BoundaryType::Vacuum; }
  bool isRefl() { return type == BoundaryType::Reflective; }
  bool isNeighb() { return type == BoundaryType::OtherElement; }
  std::string id() {
    if (type == BoundaryType::Vacuum) return "v";
    else if (type == BoundaryType::Reflective) return "r";
    else return neighbor->material_id;
  }

  // If a neighbor is present, indicate that the corner point is in the
  // points array with index pos_i
  void setCorner(int corner_i, int pos_i) {
    if (type == BoundaryType::OtherElement) {
      neighbor->corners[corner_i] = pos_i;
      neighbor->corners_set[corner_i] = true;
    }
  }
};

template<typename RefElement>
class Quadrilateral
{
public:
  static constexpr unsigned OffsetCorner = 0;
  static constexpr unsigned OffsetEdge = 4;
  static constexpr unsigned OffsetInternal = 4*(RefElement::NEdgeNodes+1);

  Quadrilateral(ParsedAssembly& xs_arg, std::string id,
      std::vector<openmc::Position>& points_, int ix, int iy) :
    material_id(id),
    xs(xs_arg),
    neighbors_set({false}),
    edge_nodes_created({false}),
    corners_set({false}),
    geom_spec_coords({ix, iy}),
    points(points_)
    {}

  std::vector<openmc::Position>& points;
  std::string material_id;
  std::array<int, 4> corners; // indices in points array
  std::array<bool, 4> corners_set;

  // x and y indices of the lattice in the geometry specification
  std::pair<int, int> geom_spec_coords;

  // May make a dedicated material class in the future for more
  // complicated XS calculations.
  ParsedAssembly& xs;

  // The constructor for ReactorGeometry fills out what the neighbors are.
  std::array<ElementBoundary<RefElement>, 4> neighbors;
  std::array<bool, 4> neighbors_set;

  // Each edge will share some nodes with a neighbor if an element greater
  // than order 1 is being used. Internal nodes are not shared, but it's
  // useful to keep them as pointers for eventually printing stuff out
  // to a VTK output file.
  std::array< std::array<int, RefElement::NEdgeNodes>, 4> edge_nodes;
  std::array<bool, 4> edge_nodes_created;
  std::array<int, RefElement::NInternalNodes> internal_nodes;

  // These methods make accessing the point list
  // more expressive.
  int bottomleft() { return corners[0]; }
  int bottomright() { return corners[1]; }
  int topright() { return corners[2]; }
  int topleft() { return corners[3]; }

  openmc::Position deltabottom() { return points[bottomright()]-points[bottomleft()]; }
  openmc::Position deltaright() { return points[topright()]-points[bottomright()]; }
  openmc::Position deltatop() { return points[topright()]-points[topleft()]; }
  openmc::Position deltaleft() { return points[topleft()]-points[bottomleft()]; }

  void setCorner(int corner_i, int pos_i) {
    corners[corner_i] = pos_i;
    corners_set[corner_i] = true;
  }

  // Access to ordinal direction neighbors. If the neighbor doesn't exist,
  // this returns a null pointer.
  Quadrilateral* getSouthWest() {
    Quadrilateral* left = neighbors[3].getNeighbor();
    Quadrilateral* down = neighbors[0].getNeighbor();
    if (left) return left->neighbors[0].getNeighbor();
    else if (down) return down->neighbors[3].getNeighbor();
    else return nullptr; // conclude that such a neighbor DNE
  }
  Quadrilateral* getSouthEast() {
    Quadrilateral* down = neighbors[0].getNeighbor();
    Quadrilateral* right = neighbors[1].getNeighbor();
    if (down) return down->neighbors[1].getNeighbor();
    else if (right) return right->neighbors[0].getNeighbor();
    else return nullptr;
  }
  Quadrilateral* getNorthEast() {
    Quadrilateral* right = neighbors[1].getNeighbor();
    Quadrilateral* up = neighbors[2].getNeighbor();
    if (right) return right->neighbors[2].getNeighbor();
    else if (up) return up->neighbors[1].getNeighbor();
    else return nullptr;
  }
  Quadrilateral* getNorthWest() {
    Quadrilateral* up = neighbors[2].getNeighbor();
    Quadrilateral* left = neighbors[3].getNeighbor();
    if (up) return up->neighbors[3].getNeighbor();
    else if (left) return left->neighbors[2].getNeighbor();
    else return nullptr;
  }

  void markCornerSet(int c) {
    Quadrilateral* ordinal;
    switch (c) {
      case 0:
        neighbors[0].setCorner(3, corners[c]);
        neighbors[3].setCorner(1, corners[c]);
        ordinal = getSouthWest();
        if (ordinal) ordinal->setCorner(2, corners[c]);
        break;
      case 1:
        neighbors[1].setCorner(0, corners[c]);
        neighbors[0].setCorner(2, corners[c]);
        ordinal = getSouthEast();
        if (ordinal) ordinal->setCorner(3, corners[c]);
        break;
      case 2:
        neighbors[1].setCorner(3, corners[c]);
        neighbors[2].setCorner(1, corners[c]);
        ordinal = getNorthEast();
        if (ordinal) ordinal->setCorner(0, corners[c]);
        break;
      case 3:
        neighbors[2].setCorner(0, corners[c]);
        neighbors[3].setCorner(2, corners[c]);
        ordinal = getNorthWest();
        if (ordinal) ordinal->setCorner(1, corners[c]);
        break;
      default:
        std::cerr << "WTF????" << std::endl;
    }
  }

  void insertEdgeNodes(int c) {
    Quadrilateral* neighbor;
    edge_nodes_created[c] = true;
    neighbor = neighbors[c].getNeighbor();
    openmc::Position dx;
    openmc::Position initial;
    switch (c) {
      case 0:
        dx = (points[corners[1]]-points[corners[0]])/RefElement::Degree;
        initial = points[corners[0]]+dx;
        if (neighbor) neighbor->edge_nodes_created[2]=true;
        break;
      case 1:
        dx = (points[corners[2]]-points[corners[1]])/RefElement::Degree;
        initial = points[corners[1]]+dx;
        if (neighbor) neighbor->edge_nodes_created[3]=true;
        break;
      case 2:
        dx = (points[corners[2]]-points[corners[3]])/RefElement::Degree;
        initial = points[corners[3]]+dx;
        if (neighbor) neighbor->edge_nodes_created[0]=true;
        break;
      case 3:
        dx = (points[corners[3]]-points[corners[0]])/RefElement::Degree;
        initial = points[corners[0]]+dx;
        if (neighbor) neighbor->edge_nodes_created[1]=true;
        break;
      default:
        std::cerr << "WTF?!?!?!" << std::endl;
    }
    for (int i=0; i<RefElement::NEdgeNodes; ++i) {
      edge_nodes[c][i] = points.size();
      points.push_back(initial);
      initial += dx;
    }
    int nc;
    switch (c) {
      case 0:
        nc = 2;
        break;
      case 1:
        nc = 3;
        break;
      case 2:
        nc = 0;
        break;
      case 3:
        nc = 1;
        break;
    }
    if (neighbor) neighbor->edge_nodes[nc] = edge_nodes[c];
  }

  void insertInternalNodes() {
    openmc::Position dx = (points[corners[1]]-points[corners[0]])/RefElement::Degree;
    openmc::Position dy = (points[corners[3]]-points[corners[0]])/RefElement::Degree;
    openmc::Position original = points[corners[0]] + dx + dy;
    for (int iy=0; iy<RefElement::NEdgeNodes; ++iy) {
      for (int ix=0; ix<RefElement::NEdgeNodes; ++ix) {
        internal_nodes[iy*RefElement::NEdgeNodes+ix] = points.size();
        points.push_back(original + ix * dx + iy * dy);
      }
    }
  }

  double area() {
    double result=0.0;
    unsigned last = corners.size()-1;
    for (int i=0; i<last; ++i) {
      result += points[corners[i]].x * points[corners[i+1]].y - points[corners[i+1]].x * points[corners[i]].y;
    }
    result += points[corners[last]].x * points[corners[0]].y - points[corners[0]].x * points[corners[last]].y;
    return 0.5 * result;
  }

  // Local-to-global enumeration of nodes
  int loc2glb(int in) {
    if (in < 0 or in >= RefElement::Npts) {
      std::cerr << "ERROR: invalid loc2glb call" << std::endl;
      exit(1);
    }
    if (in < OffsetEdge) {
      return corners[in];
    } else if (in < OffsetInternal) {
      // Determine edge it lies on
      int e = (in-OffsetEdge)/RefElement::NEdgeNodes;

      // determine index along edge
      int ei = (in-OffsetEdge)%RefElement::NEdgeNodes;

      return edge_nodes[e][ei];
    } else {
      return internal_nodes[in-OffsetInternal];
    }
  }

};


template<typename RefElement>
class ReactorGeometry
{

public:
  std::vector<openmc::Position> points;
  std::vector<Quadrilateral<RefElement>> quads;
  double dx_nom; // nominal assembly size
  double k; // eigenvalue for this geometry, saved for transients
  double eig_tol;
  unsigned n_precursors;

  // This is useful for debugging geometry. Just prints
  // each assembly geometry index, and what its neighbors are.
  void printNeighborMaterials() {
    for (auto& q:quads) {
      std::cout << q.geom_spec_coords.first << " " << q.geom_spec_coords.second << std::endl;
      std::cout << "    B:" << q.neighbors[0].id() << std::endl;
      std::cout << "    R:" << q.neighbors[1].id() << std::endl;
      std::cout << "    T:" << q.neighbors[2].id() << std::endl;
      std::cout << "    L:" << q.neighbors[3].id() << std::endl;
    }
  }

  using Bdry = ElementBoundary<RefElement>;

  ReactorGeometry(ParsedInput& input) : dx_nom(input.assembly_size),
  eig_tol(input.eig_tol)
  {
    points.reserve(RefElement::Npts * input.npoints);

    // Read in precursor count, and double-check that all materials are defining the
    // same number of precursors to be used.
    if (input.hasTransientData()) {
      n_precursors = input.getPrecursorCount();
    } else {
      n_precursors = 0;
    }

    // The geometry set routine goes like this:
    // 1) Create all the quads, assign materials, assign neighbor
    //    relations.
    // 2) Next, create all the spatial points in the geometry. Doing
    //    this second, a nicely clustered spatial ordering can be done
    //    in order to obtain a matrix with a diagonally-clustered sparsity
    //    pattern
    // 3) ??? 
    // 4) profit
    for (int row=1; row<input.core_geom.size()-1; ++row) {
      for (int col=1; col<input.core_geom[0].size()-1; ++col) {
        std::string& this_id = input.core_geom[row][col];
        if (this_id == "r" or this_id == "v") continue;
        ParsedAssembly& mat = input.lookupMaterial(this_id);
        quads.emplace_back(mat, this_id, points, col-1, row-1);
      }
    }

    // Now set the neighbors and boundary conditions of each quadrilateral (but still not any spatial points!)
    for (auto& q: quads) {
      // Need to check each neighbor. The order I'm doing is bottom, right, top, left.

      // bottom
      if (not q.neighbors_set[0]) {
        std::string& this_id = input.core_geom[q.geom_spec_coords.second][q.geom_spec_coords.first+1];
        if (this_id == "r") {
          // set reflective boundary condition
          q.neighbors[0] = {Bdry::BoundaryType::Reflective};
        } else if (this_id == "v") {
          // set vacuum boundary condition
          q.neighbors[0] = {Bdry::BoundaryType::Vacuum};
        } else {
          // Need find pointer to neighbor. Double-check it's the right one by matching up points on
          // the top and bottom. Also, set the neighbor's top neighbor to be "q".
          bool found = false;
          for (auto& qprime: quads) {
            if (qprime.geom_spec_coords.first == q.geom_spec_coords.first and q.geom_spec_coords.second-1 == qprime.geom_spec_coords.second) {
              q.neighbors[0] = {&qprime};
              qprime.neighbors[2] = {&q};
              qprime.neighbors_set[2] = true;
              found = true;
              break;
            }
          }
          if (not found) std::cerr << "can't find a bottom neighbor!!!" << std::endl;
        }
        q.neighbors_set[0] = true;
      }

      // right
      if (not q.neighbors_set[1]) {
        std::string& this_id = input.core_geom[q.geom_spec_coords.second+1][q.geom_spec_coords.first+2];
        if (this_id == "r") {
          q.neighbors[1] = {Bdry::BoundaryType::Reflective};
        } else if (this_id == "v") {
          q.neighbors[1] = {Bdry::BoundaryType::Vacuum};
        } else {
          bool found = false;
          for (auto& qprime: quads) {
            if (qprime.geom_spec_coords.second == q.geom_spec_coords.second and qprime.geom_spec_coords.first-1==q.geom_spec_coords.first) {
              q.neighbors[1] = {&qprime};
              qprime.neighbors[3] = {&q};
              qprime.neighbors_set[3] = true;
              found = true;
              break;
            }
          }
          if (not found) std::cerr << "can't find a right neighbor!!!" << std::endl;
        }
        q.neighbors_set[1] = true;
      }

      // top
      if (not q.neighbors_set[2]) {
        std::string& this_id = input.core_geom[q.geom_spec_coords.second+2][q.geom_spec_coords.first+1];
        if (this_id == "r") {
          q.neighbors[2] = {Bdry::BoundaryType::Reflective};
        } else if (this_id == "v") {
          q.neighbors[2] = {Bdry::BoundaryType::Vacuum};
        } else {
          bool found = false;
          for (auto& qprime: quads) {
            if (qprime.geom_spec_coords.second-1==q.geom_spec_coords.second and qprime.geom_spec_coords.first==qprime.geom_spec_coords.first) {
              q.neighbors[2] = {&qprime};
              qprime.neighbors[0] = {&q};
              qprime.neighbors_set[0] = true;
              found = true;
              break;
            }
          }
          if (not found) std::cerr << "can't find a top neighbor!" << std::endl;
        }
        q.neighbors_set[2] = true;
      }

      // Left
      if (not q.neighbors_set[3]) {
        std::string& this_id = input.core_geom[q.geom_spec_coords.second+1][q.geom_spec_coords.first];
        if (this_id == "r") {
          q.neighbors[3] = {Bdry::BoundaryType::Reflective};
        } else if (this_id == "v") {
          q.neighbors[3] = {Bdry::BoundaryType::Vacuum};
        } else {
          bool found = false;
          for (auto& qprime: quads) {
            if (qprime.geom_spec_coords.first+1==q.geom_spec_coords.first and qprime.geom_spec_coords.second == q.geom_spec_coords.second) {
              q.neighbors[3] = {&qprime};
              qprime.neighbors[1] = {&q};
              qprime.neighbors_set[1] = true;
              found = true;
              break;
            }
          }
          if (not found) std::cerr << "unable to find a left neighbor!" << std::endl;
        }
        q.neighbors_set[3] = true;
      }
    }

    // Now all of the points will be inserted.
    double dx = input.assembly_size;
    for (auto& q:quads) {

      // corners
      if (not q.corners_set[0]) {
        q.setCorner(0, points.size());
        q.markCornerSet(0); // this propagates that corner to all relevant neighbors who share it
        points.emplace_back(q.geom_spec_coords.first*dx, q.geom_spec_coords.second*dx, 0);
      } if (not q.corners_set[1]) {
        q.setCorner(1, points.size());
        q.markCornerSet(1);
        points.emplace_back((q.geom_spec_coords.first+1)*dx, q.geom_spec_coords.second*dx, 0);
      } if (not q.corners_set[2]) {
        q.setCorner(2, points.size());
        q.markCornerSet(2);
        points.emplace_back((q.geom_spec_coords.first+1)*dx, (q.geom_spec_coords.second+1)*dx, 0);
      } if (not q.corners_set[3]) {
        q.setCorner(3, points.size());
        q.markCornerSet(3);
        points.emplace_back(q.geom_spec_coords.first*dx, (q.geom_spec_coords.second+1)*dx, 0);
      }

      // Now add edge nodes
      if (not q.edge_nodes_created[0]) {
        q.insertEdgeNodes(0);
      } if (not q.edge_nodes_created[1]) {
        q.insertEdgeNodes(1);
      } if (not q.edge_nodes_created[2]) {
        q.insertEdgeNodes(2);
      } if (not q.edge_nodes_created[3]) {
        q.insertEdgeNodes(3);
      }

      // Finally, dealing with internal nodes is way easier since they're not shared with others.
      q.insertInternalNodes();
    }
    
    // a quick double-check for any redundant points. I made a coding
    // mistake if this loop trips on the else if
    for (int i=0; i<points.size(); ++i) {
      for (int j=0; j<points.size(); ++j) {
        if (i==j) continue;
        else if (points[i] == points[j]) std::cerr << points[i] << " " << points[j] << std::endl;
      }
    }

    // This is a simple check that no points are created which aren't referred to
    // by any quadrilateral in the mesh.
    // for (int pt=0; pt<points.size(); ++pt) {
    //   bool found = false;
    //   for (auto& q: quads) {
    //     for (auto& crn: q.corners)
    //       if (crn == pt) {found=true; break;}
    //     for (auto& neigh_array: q.edge_nodes)
    //       for (auto& pos: neigh_array)
    //         if (pos == pt) {found=true; break;}
    //     for (auto& n: q.internal_nodes)
    //       if (n == pt) {found=true;break;}
    //   }
    //   if (not found) {
    //     std::cerr << "Found an unused point!!!" << std::endl;
    //     std::cerr << pt << std::endl;
    //   }
    // }
  }

  double getMaxPointCoordinate() {
    double max = 0.0;
    for (auto& pt: points) {
      if (pt.norm() > max) max = pt.norm();
    } return max;
  }

  double getMinPointCoordinate() {
    double min = 1e100;
    for (auto& pt: points) {
      if (pt.norm() < min) min = pt.norm();
    } return min;
  }

  void printUpToPointData(std::ofstream& f) {
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
    f << "  <UnstructuredGrid>" << std::endl;
    f << "  <Piece NumberOfPoints=\""<<points.size()<<"\" NumberOfCells=\""<<quads.size()<<"\">"<<std::endl;
    f << "  <PointData>" << std::endl;
  }

  void printAfterPointData(std::ofstream& f) {
    f << "  </PointData>" << std::endl;
    f << "  <CellData>" << std::endl;
    f << "    <DataArray type=\"Int32\" Name=\"Materials\" NumberOfComponents=\"1\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"7\">" << std::endl;
    for (int i=0; i<quads.size(); ++i) {
      f << quads[i].material_id << " ";
      if (i%5 == 0 and i>0) f << std::endl;
    }
    f << "    </DataArray>" << std::endl;
    f << "  </CellData>" << std::endl;
    f << "<Points>" << std::endl;
    f << "<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\""<<getMinPointCoordinate()<<"\" RangeMax=\""<<getMaxPointCoordinate()<<"\">" << std::endl;
    for (int i=0; i<points.size(); ++i) {
      f << "    " << points[i].x << " " << points[i].y << " " << points[i].z << std::endl;
    }
    f << "</DataArray>" << std::endl;
    f << "</Points>" << std::endl;
    f << "<Cells>" << std::endl;
    f << "  <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"" << points.size() << "\">" << std::endl;
    for (auto& q:quads) {
      for (int i=0; i<q.corners.size(); ++i) f << q.corners[i] << " ";
      for (auto& arr: q.edge_nodes)
        for (int j: arr) f << j << " ";
      for (int j: q.internal_nodes) f << j << " " << std::endl;
    }
    f << "  </DataArray>" << std::endl;
    f << "  <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"" << RefElement::Npts << "\" RangeMax=\"" << RefElement::Npts*quads.size() << "\">" << std::endl;
    for (int i=1; i<quads.size()+1; ++i) f << RefElement::Npts*i << " ";
    f << std::endl;
    f << "  </DataArray>" << std::endl;
    f << "  <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"70\" RangeMax=\"70\">" << std::endl;
    for (int i=1; i<quads.size()+1; ++i) f << "70 ";
    f << std::endl;
    f << "  </DataArray>" << std::endl;
    f << "</Cells>" << std::endl;
    f << "</Piece>" << std::endl;
    f << "</UnstructuredGrid>" << std::endl;
    f << "</VTKFile>" << std::endl;
  }

  void printGeometryToVTK(std::string const& filename) {
    std::ofstream f;
    std::cout << "Writing geometry specification to " << filename << std::endl;
    f.open(filename, std::ios::out);
    if (not f.good()) std::cerr << "can't open VTK output " << filename << "!"<<std::endl;
    printUpToPointData(f);
    printAfterPointData(f);
  }

  void printDataToVTK(std::string const& filename, std::vector<Eigen::VectorXd*> data_arrays,
      std::vector<std::string*> data_names) {
    // Check that all data arrays are of the correct length
    for (auto arr: data_arrays) {
      if (arr->size() != points.size()) {
        std::cerr << "Passed an array of a size not corresponding to the DOF count in the mesh to VTK output!" << std::endl;
        std::cerr << arr->size() << " " << points.size() << std::endl;
        exit(1);
      }
    }
    if (data_names.size() != data_arrays.size()) { std::cerr << "must pass same data array name count as data arrays" << std::endl; exit(1); }
    std::ofstream f;
    std::cout << "Writing data specification to " << filename << std::endl;
    f.open(filename, std::ios::out);
    if (not f.good()) std::cerr << "can't open VTK output " << filename << "!"<<std::endl;
    printUpToPointData(f);
    for (int i=0; i<data_arrays.size(); ++i) {
      f << "    <DataArray type=\"Float32\" Name=\"" << *data_names[i] << "\" format=\"ascii\" RangeMin=\"" << data_arrays[i]->minCoeff() << "\" RangeMax=\"" << data_arrays[i]->maxCoeff() << "\">" << std::endl;
      f << "      ";
      for (int j=0; j<data_arrays[i]->size(); ++j) {
        f << data_arrays[i]->operator()(j) << " ";
        if (j%10==0 and j!=0) f << std::endl << "      ";
      }
      f << "    </DataArray>" << std::endl;
    }
    printAfterPointData(f);
  }
  unsigned npoints() { return points.size(); }
};
