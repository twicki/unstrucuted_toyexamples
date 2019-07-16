#ifndef UNSTRUCTURED_GRID_LIB_H
#define UNSTRUCTURED_GRID_LIB_H
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
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
  Vertex* v[3];
  int id_;

public:
  Triangle(int id, int dom_size_horizontal) : id_(id) {
    int gridIdx = id / 2;
    int subGridIdx = id % 2;
    int gridI = gridIdx % dom_size_horizontal;
    int gridJ = gridIdx / dom_size_horizontal;

    v[0] = Vertex::get(gridI, gridJ);
    v[1] = subGridIdx == 0 ? Vertex::get(gridI + 1, gridJ) : Vertex::get(gridI, gridJ + 1);
    v[2] = Vertex::get(gridI + 1, gridJ + 1);
  }
  int getId() const { return id_; }
  Vertex** getVertices() { return v; }
};

class Grid {
  std::vector<Triangle*> triangles_;
  std::set<Vertex*, VertexCompare> vertices_;

public:
  Grid(int size_horizontal, int size_vertical)
      : size_horizontal_(size_horizontal), size_vertical_(size_vertical),
        num_triangles_(size_horizontal * size_vertical * 2) {
    triangles_.reserve(num_triangles_);
    for(int i = 0; i < num_triangles_; ++i)
      triangles_.push_back(new Triangle(i, size_horizontal_));
    for(auto& pair : Vertex::vmap_)
      vertices_.insert(&(pair.second));
  }

  const int size_horizontal_, size_vertical_, num_triangles_;
  std::vector<Triangle*>& getTriangles() { return triangles_; }
  std::set<Vertex*, VertexCompare>& getVertices() { return vertices_; }

  std::string toVtk() {
    std::string output;
    output += "# vtk DataFile Version 3.0\n2D scalar data\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    output += printVertices();
    output += printTriangles();

    return output;
  }

  // void initGauss() {
  //   for(auto& tri : triangles_)
  //     tri->data_ = pow(double(tri->getVertices()[1]->i_ - 2), 2) +
  //                  pow(double(tri->getVertices()[1]->j_ - 2), 2);
  // }

  std::string printVertices() {
    std::string output = "POINTS " + std::to_string(vertices_.size()) + " float\n";
    for(auto v : vertices_)
      output += std::to_string(v->i_) + " " + std::to_string(v->j_) + " " + "0" + "\n";
    return output;
  }
  std::string printTriangles() {
    std::ostringstream os;
    os << "CELLS " << triangles_.size() << " " << triangles_.size() * 4 << std::endl;
    for(auto tri : triangles_) {
      auto v = tri->getVertices();
      os << "3 " << v[0]->id_ << " " << v[1]->id_ << " " << v[2]->id_ << std::endl;
    }
    os << "CELL_TYPES " << triangles_.size() << std::endl;
    for(auto tri : triangles_) {
      os << "5" << std::endl;
    }

    return os.str();
  }

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

class Data {
public:
  Data(Grid& grid);

  std::map<Vertex*, double> vertexData_;
  std::map<Triangle*, double> triangleData_;
  Grid& grid_;
  std::list<double*> data() {
    std::list<double*> output;
    for(auto tri : grid_.getTriangles()) {
      output.push_back(&(triangleData_[tri]));
    }
    return output;
  }

  std::string toVtk() {
    std::ostringstream os;
    os << "CELL_DATA " << grid_.getTriangles().size()
       << "\nSCALARS temperature  float 1\nLOOKUP_TABLE default\n";
    for(auto tri : grid_.getTriangles()) {
      os << triangleData_[tri] << std::endl;
    }
    return os.str();
  }

  void initGauss() {
    for(auto& tri : grid_.getTriangles()) {
      auto vertices = tri->getVertices();
      double center_i = double(vertices[0]->i_ + vertices[1]->i_ + vertices[2]->i_) / 3.0;
      double center_j = double(vertices[0]->j_ + vertices[1]->j_ + vertices[2]->j_) / 3.0;

      triangleData_[tri] =
          exp(-(pow(center_i - double(grid_.size_horizontal_) / 2.0, 2) +
          pow(center_j - double(grid_.size_vertical_) / 2.0, 2)));
    }

  }
};

#endif // UNSTRUCTURED_GRID_LIB_H