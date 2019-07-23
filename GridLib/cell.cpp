#include "cell.h"
#include "grid.h"

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