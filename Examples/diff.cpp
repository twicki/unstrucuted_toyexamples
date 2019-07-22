#include "unstructured_toylib.h"
#include <iostream>
#include <fstream>

int main(int argc, char const* argv[]) {
  int x = atoi(argv[1]);
  int y = atoi(argv[2]);
  std::string filename(argv[3]);
  Grid grid(x, y);

  const int timeSteps = 100;
  std::vector<CellData> temperatures;
  for(int i=0; i<timeSteps; i++)
    temperatures.emplace_back(grid);

  temperatures[0].initGauss();

  for(int step = 0; step < timeSteps-1; step++)
    for(auto triangle : grid.getTriangles()) {
      auto neighbours = grid.cellNeighboursOfCell(triangle);
      temperatures[step+1].getData()[triangle] = -0.1 * 3.0 * temperatures[step].getData()[triangle];
      for(auto neigh : neighbours)
        temperatures[step+1].getData()[triangle] += 0.1*temperatures[step].getData()[neigh];
      temperatures[step+1].getData()[triangle] += temperatures[step].getData()[triangle];
    }


  // write to files
  int step = 0;
  for(auto& temperature : temperatures) {
    std::ofstream ofile;
    ofile.open (filename + "_" + std::to_string(step) + ".vtk");
    ofile << grid.toVtk() << temperatures[step].toVtk() << std::endl;
    ofile.close();
    step++;
  }
  return 0;
}
