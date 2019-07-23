#include "grid.h"
#include <cmath>
#include <iostream>
#include <sstream>

Grid::Grid(int size_horizontal, int size_vertical)
    : size_horizontal_(size_horizontal), size_vertical_(size_vertical),
      num_triangles_(size_horizontal * size_vertical * 2) {
  triangles_.reserve(num_triangles_);
  for(int i = 0; i < num_triangles_; ++i)
    triangles_.push_back(new Triangle(i, *this));
  for(auto& pair : vmap_)
    vertices_.insert(&(pair.second));
  for(auto& pair : emap_)
    edges_.insert(&(pair.second));
  for(auto& edge : edges_) {
    auto neighbours = cellNeighboursOfEdge(edge);
    for(auto& triangle : neighbours) {
      if(triangle->getColor() == 0)
        edge->from_ = triangle;
      else
        edge->to_ = triangle;
    }
  }
}

std::list<Triangle*> Grid::cellNeighboursOfCell(Triangle* center) {
  std::list<Triangle*> out;
  int i = ((center->getId() - center->getColor()) / 2) % size_horizontal_;
  int j = ((center->getId() - center->getColor()) / 2) / size_horizontal_;

  switch(center->getColor()) {
  case 0:
    out.push_back(getTriangle(i, j, 1));
    if(j > 0)
      out.push_back(getTriangle(i, j - 1, 1));
    if(i < size_horizontal_ - 1)
      out.push_back(getTriangle(i + 1, j, 1));
    break;
  case 1:
    out.push_back(getTriangle(i, j, 0));
    if(i > 0)
      out.push_back(getTriangle(i - 1, j, 0));
    if(j < size_horizontal_ - 1)
      out.push_back(getTriangle(i, j + 1, 0));
    break;
  }
  return out;
}

std::list<Edge*> Grid::edgeNeighboursOfCell(Triangle* center) {
  std::list<Edge*> out;
  for(int i = 0; i < 3; i++)
    out.push_back(center->getEdges()[i]);
  return out;
}

std::list<Vertex*> Grid::vertexNeighboursOfCell(Triangle* center) {
  std::list<Vertex*> out;
  for(int i = 0; i < 3; i++)
    out.push_back(center->getVertices()[i]);
  return out;
}

std::list<Triangle*> Grid::cellNeighboursOfEdge(Edge* center) {
  std::list<Triangle*> out;
  int gridIdx = center->grid_i_ + size_horizontal_ * center->grid_j_;
  int adjId;
  int i = center->grid_i_;
  int j = center->grid_j_;

  switch(center->color_) {
  case 0:
    if(j >= 0)
      out.push_back(getTriangle(i, j, 1));
    if(j < size_vertical_ - 1)
      out.push_back(getTriangle(i, j + 1, 0));
    break;
  case 1:
    if(i < size_horizontal_)
      out.push_back(getTriangle(i, j, 1));
    if(i > 0)
      out.push_back(getTriangle(i - 1, j, 0));
    break;
  case 2:
    out.push_back(getTriangle(i, j, 0));
    out.push_back(getTriangle(i, j, 1));
    break;
  }
  return out;
}
std::list<Vertex*> Grid::vertexNeighboursOfEdge(Edge* center) {
  std::list<Vertex*> out;
  for(int i = 0; i < 2; i++)
    out.push_back(center->getVertices()[i]);
  return out;
}

Vertex* Grid::getVertex(int i, int j) {
  auto pair = vmap_.emplace(std::make_tuple(i, j), std::move(Vertex(i, j, vertexIdGen_)));
  if(pair.second == true)
    vertexIdGen_++;
  return &(pair.first->second);
}

Edge* Grid::getEdge(int grid_i, int grid_j, int color) {
  auto pair = emap_.emplace(std::make_tuple(grid_i, grid_j, color),
                            std::move(Edge(grid_i, grid_j, color, edgeIdGen_, *this)));
  if(pair.second == true)
    edgeIdGen_++;
  return &(pair.first->second);
}

Triangle* Grid::getTriangle(int grid_i, int grid_j, int color) {
  return triangles_[(grid_i + size_horizontal_ * grid_j) * 2 + color];
}

std::string Grid::toVtk() {
  std::string output;
  output += "# vtk DataFile Version 3.0\n2D scalar data\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  output += printVertices();
  output += printTriangles();

  return output;
}

std::string Grid::printVertices() {
  std::string output = "POINTS " + std::to_string(vertices_.size()) + " float\n";
  for(auto v : vertices_)
    output += std::to_string(v->i_) + " " + std::to_string(v->j_) + " " + "0" + "\n";
  return output;
}
std::string Grid::printEdges() {
  //! not vtk format
  std::string output = "EDGES\n";
  for(auto e : edges_)
    output += "id:" + std::to_string(e->id_) + " ijc:{" + std::to_string(e->grid_i_) + " " +
              std::to_string(e->grid_j_) + " " + std::to_string(e->color_) + "}" + "\n";
  return output;
}
std::string Grid::printTriangles() {
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

CellData::CellData(Grid& grid) : Data(grid) {
  for(const auto& cell : grid.getTriangles()) {
    data_.emplace(cell, 0);
  }
}

void CellData::initGauss(double width) {
  for(auto& tri : grid_.getTriangles()) {
    auto vertices = tri->getVertices();
    double center_i = double(vertices[0]->i_ + vertices[1]->i_ + vertices[2]->i_) / 3.0;
    double center_j = double(vertices[0]->j_ + vertices[1]->j_ + vertices[2]->j_) / 3.0;

    data_[tri] = exp(-width * (pow(center_i - double(grid_.size_horizontal_) / 2.0, 2) +
                               pow(center_j - double(grid_.size_vertical_) / 2.0, 2)));
  }
}

std::string CellData::toVtk() {
  std::ostringstream os;
  os << "CELL_DATA " << grid_.getTriangles().size()
     << "\nSCALARS temperature  float 1\nLOOKUP_TABLE default\n";
  for(auto tri : grid_.getTriangles()) {
    os << data_[tri] << std::endl;
  }
  return os.str();
}

EdgeData::EdgeData(Grid& grid) : Data(grid) {
  for(const auto& edge : grid.getEdges()) {
    data_.emplace(edge, 0);
  }
}