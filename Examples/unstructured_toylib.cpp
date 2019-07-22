#include "unstructured_toylib.h"

Edge::Edge(int grid_i, int grid_j, int color, int idGen, Grid& grid)
    : grid_i_(grid_i), grid_j_(grid_j), color_(color), id_(idGen) {
  assert(color >= 0 && color < 3);
  v_[0] = color == 0
              ? grid.getVertex(grid_i + 1, grid_j + 1)
              : (color == 1 ? grid.getVertex(grid_i, grid_j + 1) : grid.getVertex(grid_i, grid_j));
  v_[1] = color == 0 ? grid.getVertex(grid_i, grid_j + 1)
                     : (color == 1 ? grid.getVertex(grid_i, grid_j)
                                   : grid.getVertex(grid_i + 1, grid_j + 1));
};

Triangle::Triangle(int id, Grid& grid) : id_(id) {
  int gridIdx = id / 2;
  int triangleColor = id % 2; // 0 for upper, 1 for lower
  int gridI = gridIdx % grid.size_horizontal_;
  int gridJ = gridIdx / grid.size_horizontal_;

  v_[0] = grid.getVertex(gridI, gridJ);
  v_[1] = triangleColor == 0 ? grid.getVertex(gridI + 1, gridJ) : grid.getVertex(gridI, gridJ + 1);
  v_[2] = grid.getVertex(gridI + 1, gridJ + 1);

  e_[0] = triangleColor == 0 ? grid.getEdge(gridI + 1, gridJ, 1) : grid.getEdge(gridI, gridJ, 0);
  e_[1] = triangleColor == 0 ? grid.getEdge(gridI, gridJ - 1, 0) : grid.getEdge(gridI, gridJ, 1);
  e_[2] = grid.getEdge(gridI, gridJ, 2);
};

std::list<Triangle*> Grid::cellNeighboursOfCell(Triangle* center) {
  std::list<Triangle*> out;
  int adjId;
  if(center->getId() % 2 == 0) { // even

    out.push_back(triangles_[center->getId() + 1]);

    adjId = center->getId() - (size_horizontal_ * 2 - 1);
    if(triangleIdxValid(adjId))
      out.push_back(triangles_[adjId]);
    adjId = center->getId() + 3;
    if(triangleIdxValid(adjId))
      out.push_back(triangles_[adjId]);

  } else { // odd

    adjId = center->getId() - 3;
    if(triangleIdxValid(adjId))
      out.push_back(triangles_[adjId]);

    out.push_back(triangles_[center->getId() - 1]);

    adjId = center->getId() + (size_horizontal_ * 2 - 1);
    if(triangleIdxValid(adjId))
      out.push_back(triangles_[adjId]);
  }

  return out;
}

std::list<Edge*> Grid::edgeNeighboursOfCell(Triangle* center) {
  std::list<Edge*> out;
  for(int i = 0; i < 3; i++)
    out.push_back(center->getEdges()[i]);
  return out;
}

std::list<Vertex*> Grid::vertexNeighboursOfCell(Triangle* center) {
  std::list<Vertex*> out;
  for(int i = 0; i < 3; i++)
    out.push_back(center->getVertices()[i]);
  return out;
}

std::list<Triangle*> Grid::cellNeighboursOfEdge(Edge* center) {
  std::list<Triangle*> out;
  // out.push_back(center->getFromCell());
  // out.push_back(center->getToCell());
  // return out;
  int gridIdx = center->grid_i_ + size_horizontal_ * center->grid_j_;
  int adjId;

  switch(center->color_) {
  case 0:
    adjId = gridIdx * 2 + 1;
    if(triangleIdxValid(adjId))
      out.push_back(triangles_[adjId]);
    gridIdx += size_horizontal_;
    adjId = gridIdx * 2;
    if(triangleIdxValid(adjId))
      out.push_front(triangles_[adjId]);
    break;
  case 1:
    adjId = gridIdx * 2 + 1;
    if(triangleIdxValid(adjId))
      out.push_back(triangles_[adjId]);
    if(center->grid_i_ != 0) {
      gridIdx--;
      adjId = gridIdx * 2;
      out.push_front(triangles_[adjId]);
    }
    break;
  case 2:
    adjId = gridIdx * 2;
    assert(triangleIdxValid(adjId));
    out.push_back(triangles_[adjId]);
    out.push_back(triangles_[adjId + 1]);
    break;
  }
  return out;
}
//
//
//

std::list<Vertex*> Grid::vertexNeighboursOfEdge(Edge* center) {
  std::list<Vertex*> out;
  for(int i = 0; i < 2; i++)
    out.push_back(center->getVertices()[i]);
  return out;
}

CellData::CellData(Grid& grid) : Data(grid) {
  for(const auto& cell : grid.getTriangles()) {
    data_.emplace(cell, 0);
  }
}

EdgeData::EdgeData(Grid& grid) : Data(grid) {
  for(const auto& edge : grid.getEdges()) {
    data_.emplace(edge, 0);
  }
}