#include "unstructured_toylib.h"
int Vertex::idGen_;
vmap_t Vertex::vmap_;

Triangle::Triangle(int id) : id_(id) {
  int gridIdx = id / 2;
  int subGridIdx = id % 2;
  int gridI = gridIdx % DOM_SIZE_HORIZONTAL;
  int gridJ = gridIdx / DOM_SIZE_HORIZONTAL;

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

std::string Grid::printTriangles() { return "hello"; }