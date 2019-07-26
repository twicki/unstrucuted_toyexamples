#ifndef VERTEX_H
#define VERTEX_H
class Grid;

class Vertex {
  friend class Grid;
  Vertex(int i, int j, int idGen) : i_(i), j_(j), id_(idGen) {}

public:
  const int i_, j_;
  const int id_;
};

struct VertexCompare {
  bool operator()(Vertex* const& lhs, Vertex* const& rhs) { return lhs->id_ < rhs->id_; }
};

#endif // VERTEX_H
