#include "unstructured_toylib.h"

#include <sstream>
int Vertex::idGen_;
vmap_t Vertex::vmap_;

Triangle::Triangle(int id, int dom_size_horizontal) : id_(id) {
  int gridIdx = id / 2;
  int subGridIdx = id % 2;
  int gridI = gridIdx % dom_size_horizontal;
  int gridJ = gridIdx / dom_size_horizontal;

  v[0] = Vertex::get(gridI, gridJ);
  v[1] = subGridIdx == 0 ? Vertex::get(gridI + 1, gridJ) : Vertex::get(gridI, gridJ + 1);
  v[2] = Vertex::get(gridI + 1, gridJ + 1);
}

std::string Grid::toVtk() {
  std::string output;
  output += "# vtk DataFile Version 3.0\n2D scalar data\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  output += printVertices();
  output += printTriangles();

  return output;
}

std::string Grid::printVertices() {
  std::string output = "POINTS " + std::to_string(vertices_.size()) + " float\n";
  for(auto v : vertices_)
    output += std::to_string(v->i_) + " " + std::to_string(v->j_) + " " + "0" + "\n";
  return output;
}

std::string Grid::printTriangles() {
  std::ostringstream os; 
  os << "CELLS " << triangles_.size() << " " << triangles_.size() * 4 << std::endl;
  for(auto tri : triangles_) {
    auto v = tri->getVertices();
    os << "3 " << v[0]->id_ << " " << v[1]->id_ << " " << v[2]->id_ << std::endl;
  }
  os << "CELL_TYPES " << triangles_.size() << std::endl;
  for(auto tri : triangles_) {
    os << "5" << std::endl;
  }
  os << "CELL_DATA " << triangles_.size() << "\nSCALARS temperature  float 1\nLOOKUP_TABLE default\n";
  for(auto tri : triangles_) {
    os << tri->data_ << std::endl;
  }
  return os.str();
}