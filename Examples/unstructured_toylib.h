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
#include <assert.h>

typedef std::tuple<int, int> coord_t;
typedef std::tuple<int, int, int> edge_coord_t;

struct vertex_key_hash : public std::unary_function<coord_t, std::size_t> {
  std::size_t operator()(const coord_t& k) const { return std::get<0>(k) ^ std::get<1>(k); }
};
struct vertex_key_equal : public std::binary_function<coord_t, coord_t, bool> {
  bool operator()(const coord_t& v0, const coord_t& v1) const {
    return (std::get<0>(v0) == std::get<0>(v1) && std::get<1>(v0) == std::get<1>(v1));
  }
};
struct edge_key_hash : public std::unary_function<edge_coord_t, std::size_t> {
  std::size_t operator()(const edge_coord_t& k) const { return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k); }
};
struct edge_key_equal : public std::binary_function<edge_coord_t, edge_coord_t, bool> {
  bool operator()(const edge_coord_t& v0, const edge_coord_t& v1) const {
    return (std::get<0>(v0) == std::get<0>(v1) && std::get<1>(v0) == std::get<1>(v1) && std::get<2>(v0) == std::get<2>(v1));
  }
};

class Grid;
class Vertex;
class Edge;
class Triangle;
typedef std::unordered_map<const coord_t, Vertex, vertex_key_hash, vertex_key_equal> vmap_t;
typedef std::unordered_map<const edge_coord_t, Edge, edge_key_hash, edge_key_equal> emap_t;

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

class Edge {
  friend class Grid;
  Vertex * v_[2];
  Edge(int grid_i, int grid_j, int color, int idGen, Grid& grid);

public:
  const int grid_i_, grid_j_, color_; //color 0: edge below of grid cell, 1: left edge of grid cell, 2: diagonal edge of grid cell
  const int id_;
  inline Vertex** getVertices() { return v_; }
};

struct EdgeCompare {
  bool operator()(Edge* const& lhs, Edge* const& rhs) { return lhs->id_ < rhs->id_; }
};

class Triangle {
  Vertex* v_[3];
  Edge* e_[3];
  int id_;

public:
  Triangle(int id, Grid& grid);
  int getId() const { return id_; }
  inline Vertex** getVertices() { return v_; }
  inline Edge** getEdges() { return e_; }
};

class Grid {
  friend class Triangle;
  friend class Edge;
  std::vector<Triangle*> triangles_;
  std::set<Vertex*, VertexCompare> vertices_;
  std::set<Edge*, EdgeCompare> edges_;
  vmap_t vmap_;
  emap_t emap_;
  int vertexIdGen_=0, edgeIdGen_=0;

  Vertex* getVertex(int i, int j) {
    auto pair = vmap_.emplace(std::make_tuple(i, j), std::move(Vertex(i, j, vertexIdGen_)));
    if(pair.second == true)
      vertexIdGen_++;
    return &(pair.first->second);
  }

  Edge* getEdge(int grid_i, int grid_j, int color) {
    auto pair = emap_.emplace(std::make_tuple(grid_i, grid_j, color), std::move(Edge(grid_i, grid_j, color, edgeIdGen_, *this)));
    if(pair.second == true)
      edgeIdGen_++;
    return &(pair.first->second);
  }

public:

  const int size_horizontal_, size_vertical_, num_triangles_;

  Grid(int size_horizontal, int size_vertical)
      : size_horizontal_(size_horizontal), size_vertical_(size_vertical),
        num_triangles_(size_horizontal * size_vertical * 2) {
    triangles_.reserve(num_triangles_);
    for(int i = 0; i < num_triangles_; ++i)
      triangles_.push_back(new Triangle(i, *this));
    for(auto& pair : vmap_)
      vertices_.insert(&(pair.second));
    for(auto& pair : emap_)
      edges_.insert(&(pair.second));
  }

  inline std::vector<Triangle*>& getTriangles() { return triangles_; }
  inline std::set<Vertex*, VertexCompare>& getVertices() { return vertices_; }
  inline std::set<Edge*, EdgeCompare>& getEdges() { return edges_; }

  std::string toVtk() {
    std::string output;
    output += "# vtk DataFile Version 3.0\n2D scalar data\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    output += printVertices();
    output += printTriangles();

    return output;
  }

  std::string printVertices() {
    std::string output = "POINTS " + std::to_string(vertices_.size()) + " float\n";
    for(auto v : vertices_)
      output += std::to_string(v->i_) + " " + std::to_string(v->j_) + " " + "0" + "\n";
    return output;
  }
  std::string printEdges() {
    //! not vtk format
    std::string output = "EDGES\n";
    for(auto e : edges_)
      output += "id:" + std::to_string(e->id_) + " ijc:{" + std::to_string(e->grid_i_) + " " + std::to_string(e->grid_j_) + " " + std::to_string(e->color_) + "}" + "\n";
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

  std::list<Triangle*> cellNeighboursOfCell(Triangle* center);
  std::list<Edge*> edgeNeighboursOfCell(Triangle* center);
  std::list<Vertex*> vertexNeighboursOfCell(Triangle* center);
  std::list<Triangle*> cellNeighboursOfEdge(Edge* center); // returns first an upper triangle then a lower triangle
  // std::list<Edge*> edgeNeighboursOfEdge(Edge* center);
  std::list<Vertex*> vertexNeighboursOfEdge(Edge* center);
  // std::list<Triangle*> cellNeighboursOfVertex(Vertex* center);
  // std::list<Edge*> edgeNeighboursOfVertex(Vertex* center);
  // std::list<Vertex*> vertexNeighboursOfVertex(Vertex* center);
  
};

template<typename T>
class Data {
protected:
  std::map<T*, double> data_;
  Grid& grid_;
  Data(Grid& grid) : grid_(grid) {}
public:
  virtual std::string toVtk() = 0;
  std::map<T*, double>& getData() {
    return data_;
  }
};

class CellData : public Data<Triangle> {
public:
  CellData(Grid& grid);

  std::string toVtk() {
    std::ostringstream os;
    os << "CELL_DATA " << grid_.getTriangles().size()
       << "\nSCALARS temperature  float 1\nLOOKUP_TABLE default\n";
    for(auto tri : grid_.getTriangles()) {
      os << data_[tri] << std::endl;
    }
    return os.str();
  }

  void initGauss() {
    for(auto& tri : grid_.getTriangles()) {
      auto vertices = tri->getVertices();
      double center_i = double(vertices[0]->i_ + vertices[1]->i_ + vertices[2]->i_) / 3.0;
      double center_j = double(vertices[0]->j_ + vertices[1]->j_ + vertices[2]->j_) / 3.0;

      data_[tri] =
          exp(-0.1*(pow(center_i - double(grid_.size_horizontal_) / 2.0, 2) +
          pow(center_j - double(grid_.size_vertical_) / 2.0, 2)));
    }

  }
};

#endif // UNSTRUCTURED_GRID_LIB_H