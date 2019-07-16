#ifndef UNSTRUCTURED_GRID_LIB_H
#define UNSTRUCTURED_GRID_LIB_H
#include <cmath>
#include <iostream>
#include <list>
#include <ostream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#define DOM_SIZE_HORIZONTAL 3
#define DOM_SIZE_VERTICAL 3
#define NUM_TRIANGLES DOM_SIZE_HORIZONTAL* DOM_SIZE_VERTICAL * 2

typedef std::tuple<int, int> coord_t;

struct key_hash : public std::unary_function<coord_t, std::size_t> {
  std::size_t operator()(const coord_t& k) const { return std::get<0>(k) ^ std::get<1>(k); }
};

struct key_equal : public std::binary_function<coord_t, coord_t, bool> {
  bool operator()(const coord_t& v0, const coord_t& v1) const {
    return (std::get<0>(v0) == std::get<0>(v1) && std::get<1>(v0) == std::get<1>(v1));
  }
};

class Grid;
class Vertex;
typedef std::unordered_map<const coord_t, Vertex, key_hash, key_equal> vmap_t;

class Vertex {

  friend class Grid;

  static vmap_t vmap_;
  static int idGen_;
  Vertex(int i, int j) : i_(i), j_(j), id_(idGen_) {}

public:
  static Vertex* get(int i, int j) {
    auto pair = vmap_.emplace(std::make_tuple(i, j), std::move(Vertex(i, j)));
    if(pair.second == true)
      idGen_++;
    return &(pair.first->second);
  }
  int i_, j_;
  int id_;
};

struct VertexCompare {
  bool operator()(Vertex* const& lhs, Vertex* const& rhs) { return lhs->id_ < rhs->id_; }
};

class Triangle {
  int id_;
  Vertex* v[3];

public:
  Triangle(int id);
  int getId() const { return id_; }
  Vertex** getVertices() { return v; }
  double data_;
};
template <typename collection>
std::ostream& printVertices(std::ostream& os, const collection& list) {
  os << "POINTS " << list.size() << " float" << std::endl;
  for(auto v : list)
    os << v->i_ << " " << v->j_ << " "
       << "0" << std::endl;
  return os;
}
template <typename collection>
std::ostream& printTriangles(std::ostream& os, const collection& list) {
  os << "CELLS " << list.size() << " " << list.size() * 4 << std::endl;
  for(auto tri : list) {
    auto v = tri->getVertices();
    os << "3 " << v[0]->id_ << " " << v[1]->id_ << " " << v[2]->id_ << std::endl;
  }
  os << "CELL_TYPES " << list.size() << std::endl;
  for(auto tri : list) {
    os << "5" << std::endl;
  }
  os << "CELL_DATA " << list.size() << "\nSCALARS temperature  float 1\nLOOKUP_TABLE default\n";
  for(auto tri : list) {
    os << tri->data_ << std::endl;
  }
  return os;
}

class Grid {

  std::vector<Triangle*> triangles_;
  std::set<Vertex*, VertexCompare> vertices_;

  void init() {
    triangles_.reserve(NUM_TRIANGLES);
    for(int i = 0; i < NUM_TRIANGLES; ++i)
      triangles_.push_back(new Triangle(i));
    for(auto& tri : triangles_)
      tri->data_ = pow(double(tri->getVertices()[1]->i_ - 2), 2) +
                   pow(double(tri->getVertices()[1]->j_ - 2), 2);
    for(auto& pair : Vertex::vmap_)
      vertices_.insert(&(pair.second));
  }

public:
  Grid() { init(); }
  std::vector<Triangle*>& getTriangles() { return triangles_; }
  std::set<Vertex*, VertexCompare>& getVertices() { return vertices_; }
  void toVtk();

  bool triangleIdxValid(int idx) { return idx >= 0 && idx < NUM_TRIANGLES; }

  std::list<Triangle*> adjacencyList(Triangle& center) {
    std::list<Triangle*> out;
    int adjId;
    if(center.getId() % 2 == 0) { // even

      out.push_back(triangles_[center.getId() + 1]);

      adjId = center.getId() - (DOM_SIZE_HORIZONTAL * 2 - 1);
      if(triangleIdxValid(adjId))
        out.push_back(triangles_[adjId]);
      adjId = center.getId() + 3;
      if(triangleIdxValid(adjId))
        out.push_back(triangles_[adjId]);

    } else { // odd

      adjId = center.getId() - 3;
      if(triangleIdxValid(adjId))
        out.push_back(triangles_[adjId]);

      out.push_back(triangles_[center.getId() - 1]);

      adjId = center.getId() + (DOM_SIZE_HORIZONTAL * 2 - 1);
      if(triangleIdxValid(adjId))
        out.push_back(triangles_[adjId]);
    }

    return out;
  }
};

#endif // UNSTRUCTURED_GRID_LIB_H