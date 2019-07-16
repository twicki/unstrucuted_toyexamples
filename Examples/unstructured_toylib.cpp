#include "unstructured_toylib.h"

int Vertex::idGen_;
vmap_t Vertex::vmap_;

Data::Data(Grid& grid) : grid_(grid) {
  for(const auto& vtx : grid.getVertices()) {
    vertexData_.emplace(vtx, 0);
  }
  for(const auto& cell : grid.getTriangles()) {
    triangleData_.emplace(cell, 0);
  }
}