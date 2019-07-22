#include "unstructured_toylib.h"
#include <fstream>
#include <iostream>

int main(int argc, char const* argv[]) {
  int x = atoi(argv[1]);
  int y = atoi(argv[2]);
  std::string filename(argv[3]);
  Grid grid(x, y);

  const int timeSteps = 1000;
  std::vector<CellData> temperatures;
  std::vector<EdgeData> flux;
  for(int i = 0; i < timeSteps; i++) {
    temperatures.emplace_back(grid);
    flux.emplace_back(grid);
  }

  temperatures[0].initGauss(0.01);

  // TODO: calculate gradient of divergence
  for(int step = 0; step < timeSteps - 1; step++) {
    for(auto edge : grid.getEdges()) {
      flux[step].getData()[edge] = grid.cellNeighboursOfEdge(edge).size() > 1
                                       ? temperatures[step].getData()[edge->getFromCell()] -
                                             temperatures[step].getData()[edge->getToCell()]
                                       : 0.0;
    }
    for(auto triangle : grid.getTriangles()) {
      auto neighbours = grid.edgeNeighboursOfCell(triangle);
      temperatures[step + 1].getData()[triangle] = temperatures[step].getData()[triangle];

      for(auto edge : neighbours)
        temperatures[step + 1].getData()[triangle] +=
            (edge->getFromCell() == triangle ? -1.0 : 1.0) * (0.1) * flux[step].getData()[edge];
    }
  }

  // for(auto triangle : grid.getTriangles()) {
  //   auto neighbours = grid.cellNeighboursOfCell(triangle);
  //   temperatures[step+1].getData()[triangle] = -0.1 * 3.0 *
  //   temperatures[step].getData()[triangle]; for(auto neigh : neighbours)
  //     temperatures[step+1].getData()[triangle] +=
  //     0.1*temperatures[step].getData()[neigh];
  //   temperatures[step+1].getData()[triangle] +=
  //   temperatures[step].getData()[triangle];
  // }

  // write to file
  for(int step = 0; step < timeSteps - 1; step++) {
    std::ofstream ofile;
    ofile.open(filename + "_" + std::to_string(step) + ".vtk");
    ofile << grid.toVtk() << temperatures[step].toVtk() << std::endl;
    ofile.close();
  }
  return 0;
}
