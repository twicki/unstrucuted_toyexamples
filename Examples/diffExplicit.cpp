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
  for(int i = 0; i < timeSteps; i++)
    temperatures.emplace_back(grid);

  // initial conditions
  temperatures[0].initGauss(0.01);

  // time stepping
  for(int step = 0; step < timeSteps - 1; step++)
    for(auto triangle : grid.getTriangles()) {
      auto neighbours = grid.cellNeighboursOfCell(triangle);
      temperatures[step + 1].getData(triangle) =
          -0.1 * double(neighbours.size()) * temperatures[step].getData(triangle);
      for(auto neigh : neighbours)
        temperatures[step + 1].getData(triangle) += 0.1 * temperatures[step].getData(neigh);
      temperatures[step + 1].getData(triangle) += temperatures[step].getData(triangle);
    }

  // write to files
  int step = 0;
  for(auto& temperature : temperatures) {
    std::ofstream ofile;
    ofile.open(filename + "_" + std::to_string(step) + ".vtk");
    ofile << grid.toVtk() << temperatures[step].toVtk() << std::endl;
    ofile.close();
    step++;
  }
  return 0;
}
