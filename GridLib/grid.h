#ifndef UNSTRUCTURED_GRID_LIB_H
#define UNSTRUCTURED_GRID_LIB_H
#include "cell.h"
#include "edge.h"
#include "vertex.h"
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

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
  std::size_t operator()(const edge_coord_t& k) const {
    return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k);
  }
};
struct edge_key_equal : public std::binary_function<edge_coord_t, edge_coord_t, bool> {
  bool operator()(const edge_coord_t& v0, const edge_coord_t& v1) const {
    return (std::get<0>(v0) == std::get<0>(v1) && std::get<1>(v0) == std::get<1>(v1) &&
            std::get<2>(v0) == std::get<2>(v1));
  }
};

typedef std::unordered_map<const coord_t, Vertex, vertex_key_hash, vertex_key_equal> vmap_t;
typedef std::unordered_map<const edge_coord_t, Edge, edge_key_hash, edge_key_equal> emap_t;

class Grid {
  friend class Triangle;
  friend class Edge;
  std::vector<Triangle*> triangles_;
  std::set<Vertex*, VertexCompare> vertices_;
  std::set<Edge*, EdgeCompare> edges_;
  vmap_t vmap_;
  emap_t emap_;
  int vertexIdGen_ = 0;
  int edgeIdGen_ = 0;

public:
  Grid(int size_horizontal, int size_vertical);

  const int size_horizontal_, size_vertical_, num_triangles_;

  std::vector<Triangle*>& getTriangles() { return triangles_; }
  std::set<Vertex*, VertexCompare>& getVertices() { return vertices_; }
  std::set<Edge*, EdgeCompare>& getEdges() { return edges_; }
  Vertex* getVertex(int i, int j);
  Edge* getEdge(int grid_i, int grid_j, int color);
  Triangle* getTriangle(int grid_i, int grid_j, int color);

  std::string toVtk();

  std::list<Triangle*> cellNeighboursOfCell(Triangle* center);
  std::list<Edge*> edgeNeighboursOfCell(Triangle* center);
  std::list<Vertex*> vertexNeighboursOfCell(Triangle* center);
  std::list<Triangle*> cellNeighboursOfEdge(Edge* center);
  // std::list<Edge*> edgeNeighboursOfEdge(Edge* center);
  std::list<Vertex*> vertexNeighboursOfEdge(Edge* center);
  // std::list<Triangle*> cellNeighboursOfVertex(Vertex* center);
  // std::list<Edge*> edgeNeighboursOfVertex(Vertex* center);
  // std::list<Vertex*> vertexNeighboursOfVertex(Vertex* center);
private:
  std::string printVertices();
  std::string printEdges();
  std::string printTriangles();
};

template <typename T>
class Data {
protected:
  std::map<T*, double> data_;
  Grid& grid_;
  Data(Grid& grid) : grid_(grid) {}

public:
  virtual std::string toVtk() = 0;
  std::map<T*, double>& getData() { return data_; }
  double& getData(T* key) { return data_[key]; }
};

class CellData : public Data<Triangle> {
public:
  CellData(Grid& grid);

  std::string toVtk();

  void initGauss(double width);
};

class EdgeData : public Data<Edge> {
public:
  EdgeData(Grid& grid);

  std::string toVtk() {
    return ""; // there's no edge data in vtk format
  }
};

#endif // UNSTRUCTURED_GRID_LIB_H