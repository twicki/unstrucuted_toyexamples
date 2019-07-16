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

std::ostream& operator<<(std::ostream& os, const std::list<Vertex*>& list) {
  return printVertices(os, list);
};
std::ostream& operator<<(std::ostream& os, const std::set<Vertex*, VertexCompare>& list) {
  return printVertices(os, list);
};
std::ostream& operator<<(std::ostream& os, const std::list<Triangle*>& list) {
  return printTriangles(os, list);
};
std::ostream& operator<<(std::ostream& os, const std::vector<Triangle*>& list) {
  return printTriangles(os, list);
};

void Grid::toVtk() {
  std::cout << "# vtk DataFile Version 3.0\n2D scalar data\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  // print vertices
  std::cout << vertices_ << std::endl;
  // print triangles
  std::cout << triangles_ << std::endl;
}