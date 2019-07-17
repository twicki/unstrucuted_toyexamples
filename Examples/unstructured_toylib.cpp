#include "unstructured_toylib.h"

Edge::Edge(int grid_i, int grid_j, int color, int idGen, Grid& grid)
   : grid_i_(grid_i), grid_j_(grid_j), color_(color), id_(idGen) { 
    assert(color>=0 && color<3); 
    v_[0] = color==0 ? grid.getVertex(grid_i+1, grid_j+1) : (color==1 ? grid.getVertex(grid_i, grid_j+1) : grid.getVertex(grid_i, grid_j) );
    v_[1] = color==0 ? grid.getVertex(grid_i, grid_j+1) : (color==1 ? grid.getVertex(grid_i, grid_j) : grid.getVertex(grid_i+1, grid_j+1) );
  };

Triangle::Triangle(int id, Grid& grid) : id_(id) {
    int gridIdx = id / 2;
    int triangleColor = id % 2; // 0 for upper, 1 for lower
    int gridI = gridIdx % grid.size_horizontal_;
    int gridJ = gridIdx / grid.size_horizontal_;

    v_[0] = grid.getVertex(gridI, gridJ);
    v_[1] = triangleColor == 0 ? grid.getVertex(gridI + 1, gridJ) : grid.getVertex(gridI, gridJ + 1);
    v_[2] = grid.getVertex(gridI + 1, gridJ + 1);

    e_[0] = triangleColor == 0 ? grid.getEdge(gridI+1, gridJ, 1) : grid.getEdge(gridI, gridJ, 0);
    e_[1] = triangleColor == 0 ? grid.getEdge(gridI, gridJ-1, 0) : grid.getEdge(gridI, gridJ, 1);
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

Data::Data(Grid& grid) : grid_(grid) {
  for(const auto& vtx : grid.getVertices()) {
    vertexData_.emplace(vtx, 0);
  }
  for(const auto& cell : grid.getTriangles()) {
    triangleData_.emplace(cell, 0);
  }
}