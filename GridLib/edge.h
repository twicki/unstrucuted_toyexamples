#ifndef EDGE_H
#define EDGE_H

class Grid;
class Triangle;
class Vertex;

class Edge {
  friend class Grid;
  Vertex* v_[2];
  Triangle* from_ = nullptr;
  Triangle* to_ = nullptr;

  Edge(int grid_i, int grid_j, int color, int idGen, Grid& grid);

public:
  const int grid_i_, grid_j_;
  const int color_; // color 0: edge below of grid cell, 1: left edge of grid
                    // cell, 2: diagonal edge of grid cell
  const int id_;
  inline Vertex** getVertices() { return v_; }
  inline Triangle* getFromCell() { return from_; }
  inline Triangle* getToCell() { return to_; }
  inline int getColor() { return color_; }
};

struct EdgeCompare {
  bool operator()(Edge* const& lhs, Edge* const& rhs) { return lhs->id_ < rhs->id_; }
};

#endif // EDGE_H
