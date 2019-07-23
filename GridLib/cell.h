#ifndef CELL_H
#define CELL_H
class Vertex;
class Edge;
class Grid;

class Triangle {
  Vertex* v_[3];
  Edge* e_[3];
  int id_;

public:
  Triangle(int id, Grid& grid);
  int getId() const { return id_; }
  inline Vertex** getVertices() { return v_; }
  inline Edge** getEdges() { return e_; }
  inline int getColor() { return id_ % 2; }
};

#endif // CELL_H