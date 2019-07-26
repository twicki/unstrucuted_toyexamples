#include "grid.h"
#include <cmath>
#include <fstream>
#include <iostream>

int main(int argc, char const* argv[]) {
  if(argc != 5) {
    std::cout << "intended use is\n"
              << argv[0] << " <X> <Y> <# Timesteps> <outputfile base>" << std::endl;
    return -1;
  }
  int x = atoi(argv[1]);
  int y = atoi(argv[2]);
  int timeSteps = atoi(argv[3]);
  std::string filename(argv[4]);
  Grid grid(x, y);

  // setup of one grid per timestep
  std::vector<EdgeData> temperatures;
  std::vector<CellData> gradients;
  std::vector<VertexData> curls;
  for(int i = 0; i < timeSteps; i++) {
    temperatures.emplace_back(grid);
    gradients.emplace_back(grid);
    curls.emplace_back(grid);
  }

  // initial conditions
  temperatures[0].initGauss(0.01);
  double area = 0.5;
  double dual_area = 1;
  double dual_length = 0.5;

  for(int time = 0; time < timeSteps - 1; time++) {
    //==--------------------------------------------------------------------------------------------
    // computation of the gradient
    //==--------------------------------------------------------------------------------------------
    for(auto cell : grid.getTriangles()) {
      auto edgeNeighbours = grid.edgeNeighboursOfCell(cell);
      gradients[time].getData(cell) = 0;
      for(auto edge : edgeNeighbours) {
        double edgelength = (edge->getColor() == 2) ? std::sqrt(2) : 1;
        gradients[time].getData(cell) += temperatures[time].getData(edge) *
                                         (edge->getFromCell() == cell ? -1.0 : 1.0) * edgelength /
                                         area;
      }
    }
    //==--------------------------------------------------------------------------------------------
    // computation of the curl
    //==--------------------------------------------------------------------------------------------
    for(auto vertex : grid.getVertices()) {
      double curl = 0;
      auto edgeNeighbours = grid.edgeNeighboursOfVertex(vertex);
      auto ptr = edgeNeighbours.front();
      // 0 - 3 (north - south)
      // 1 - 4 (east - west)
      // 2 - 5 (SE - NW)
      for(int i = 0; i < 3; ++i) {
        curl += temperatures[time].getData(ptr) - temperatures[time].getData(ptr + 3);
        ptr++;
      }
      curl *= (dual_length / dual_area);
      curls[time].getData(vertex) = curl;
    }
    //==--------------------------------------------------------------------------------------------
    // computation of the laplacian
    //==--------------------------------------------------------------------------------------------
    for(auto edge : grid.getEdges()) {
      double grad, curl;
      double edgelength = (edge->getColor() == 2) ? std::sqrt(2) : 1;
      grad = curl = 0;
      //==------------------------------------------------------------------------------------------
      // divergence of the gradient
      //==------------------------------------------------------------------------------------------
      auto cellNeighbours = grid.cellNeighboursOfEdge(edge);
      for(auto cell : cellNeighbours) {
        int sign = edge->getFromCell() == cell ? -1 : 1;
        grad += gradients[time].getData(cell) * sign;
      }
      grad /= dual_length;
      //==------------------------------------------------------------------------------------------
      // divergence of the curl
      //==------------------------------------------------------------------------------------------
      auto vertexNeighbours = grid.vertexNeighboursOfEdge(edge);
      auto front = vertexNeighbours.front();
      auto back = vertexNeighbours.back();
      curl = (curls[time].getData(front) - curls[time].getData(back));
      curl /= edgelength;
      temperatures[time + 1].getData(edge) = curl - grad;
    }
  }

  // write to files
  for(int step = 0; step < temperatures.size(); ++step) {
    std::ofstream ofile;
    ofile.open(filename + "_" + std::to_string(step) + ".vtk");
    ofile << grid.toVtk() << temperatures[step].toVtk() << std::endl;
    ofile.close();
  }
  return 0;
}
