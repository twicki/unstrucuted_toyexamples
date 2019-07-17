#include "unstructured_toylib.h"
#include <iostream>
#include <fstream>

int main(int argc, char const* argv[]) {
  /* code */
  // std::cout << "hello" << std::endl;
  int x = atoi(argv[1]);
  int y = atoi(argv[2]);
  std::string filename(argv[3]);
  Grid grid(x, y);

  const int timeSteps = 2;
  std::vector<Data> winds;
  for(int i=0; i<timeSteps; i++)
    winds.emplace_back(grid);

  winds[0].initGauss();


  for(int step = 0; step < timeSteps-1; step++)
    for(auto triangle : grid.getTriangles()) {
      auto neighbours = grid.cellNeighboursOfCell(triangle);
      //grid.vertexNeighbourOfCell
      winds[step+1].triangleData_[triangle] = -0.1 * 3.0 * winds[step].triangleData_[triangle];
      for(auto neigh : neighbours)
        winds[step+1].triangleData_[triangle] += 0.1*winds[step].triangleData_[neigh];
      winds[step+1].triangleData_[triangle] += winds[step].triangleData_[triangle];
    }


  // write to files
  int step = 0;
  for(auto& wind : winds) {
    std::ofstream ofile;
    ofile.open (filename + "_" + std::to_string(step) + ".vtk");
    ofile << grid.toVtk() << winds[step].toVtk() << std::endl;
    ofile.close();
    step++;
  }
  return 0;
}
