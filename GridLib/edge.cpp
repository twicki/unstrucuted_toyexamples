#include "edge.h"
#include "grid.h"
#include <assert.h>

Edge::Edge(int grid_i, int grid_j, int color, int idGen, Grid& grid)
    : grid_i_(grid_i), grid_j_(grid_j), color_(color), id_(idGen) {
  assert(color >= 0 && color < 3);
  v_[0] = color == 0
              ? grid.getVertex(grid_i + 1, grid_j + 1)
              : (color == 1 ? grid.getVertex(grid_i, grid_j + 1) : grid.getVertex(grid_i, grid_j));
  v_[1] = color == 0 ? grid.getVertex(grid_i, grid_j + 1)
                     : (color == 1 ? grid.getVertex(grid_i, grid_j)
                                   : grid.getVertex(grid_i + 1, grid_j + 1));
}
