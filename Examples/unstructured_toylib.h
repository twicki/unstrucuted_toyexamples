#ifndef UNSTRUCTURED_GRID_LIB_H
#define UNSTRUCTURED_GRID_LIB_H
#include <cmath>
#include <iostream>
#include <list>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

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
  Triangle(int id, int dom_size_horizontal);
  int getId() const { return id_; }
  Vertex** getVertices() { return v; }
  double data_;
};

class Grid {
  const int size_horizontal_, size_vertical_, num_triangles_;
  std::vector<Triangle*> triangles_;
  std::set<Vertex*, VertexCompare> vertices_;

  void init() {
    triangles_.reserve(num_triangles_);
    for(int i = 0; i < num_triangles_; ++i)
      triangles_.push_back(new Triangle(i, size_horizontal_));
    for(auto& tri : triangles_)
      tri->data_ = pow(double(tri->getVertices()[1]->i_ - 2), 2) +
                   pow(double(tri->getVertices()[1]->j_ - 2), 2);
    for(auto& pair : Vertex::vmap_)
      vertices_.insert(&(pair.second));
  }

public:
  Grid(int size_horizontal, int size_vertical) : 
    size_horizontal_(size_horizontal), size_vertical_(size_vertical), num_triangles_(size_horizontal * size_vertical * 2) { init(); }
  std::vector<Triangle*>& getTriangles() { return triangles_; }
  std::set<Vertex*, VertexCompare>& getVertices() { return vertices_; }
  std::string toVtk();

  std::string printVertices();
  std::string printTriangles();

  bool triangleIdxValid(int idx) { return idx >= 0 && idx < num_triangles_; }

  std::list<Triangle*> adjacencyList(Triangle& center) {
    std::list<Triangle*> out;
    int adjId;
    if(center.getId() % 2 == 0) { // even

      out.push_back(triangles_[center.getId() + 1]);

      adjId = center.getId() - (size_horizontal_ * 2 - 1);
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

      adjId = center.getId() + (size_horizontal_ * 2 - 1);
      if(triangleIdxValid(adjId))
        out.push_back(triangles_[adjId]);
    }

    return out;
  }
};

#endif // UNSTRUCTURED_GRID_LIB_H