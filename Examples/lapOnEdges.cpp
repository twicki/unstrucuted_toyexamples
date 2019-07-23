#include "grid.h"
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
  std::vector<CellData> temperatures;
  std::vector<EdgeData> flux;
  for(int i = 0; i < timeSteps; i++) {
    temperatures.emplace_back(grid);
    flux.emplace_back(grid);
  }

  // initial conditions
  temperatures[0].initGauss(0.01);

  // time stepping
  for(int step = 0; step < timeSteps - 1; step++) {
    for(auto edge : grid.getEdges()) {
      flux[step].getData(edge) = grid.cellNeighboursOfEdge(edge).size() > 1
                                     ? temperatures[step].getData(edge->getFromCell()) -
                                           temperatures[step].getData(edge->getToCell())
                                     : 0.0;
    }
    for(auto triangle : grid.getTriangles()) {
      auto neighbours = grid.edgeNeighboursOfCell(triangle);
      temperatures[step + 1].getData(triangle) = temperatures[step].getData(triangle);

      for(auto edge : neighbours) {
        double scaling = 1;
        //     edge->getColor() == 2 ? 3.0 * sqrt(2) / (2.0 + sqrt(2)) : 3.0 / (2.0 + sqrt(2));
        temperatures[step + 1].getData(triangle) += (edge->getFromCell() == triangle ? -1.0 : 1.0) *
                                                    (0.1) * scaling * flux[step].getData(edge);
      }
    }
  }

  // write to file
  for(int step = 0; step < timeSteps; step++) {
    std::ofstream ofile;
    ofile.open(filename + "_" + std::to_string(step) + ".vtk");
    ofile << grid.toVtk() << temperatures[step].toVtk() << std::endl;
    ofile.close();
  }
  return 0;
}
