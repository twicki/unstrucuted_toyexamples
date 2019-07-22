#include "unstructured_toylib.h"
#include <iostream>
#include <fstream>

int main(int argc, char const* argv[]) {
  int x = atoi(argv[1]);
  int y = atoi(argv[2]);
  std::string filename(argv[3]);
  Grid grid(x, y);

  CellData field_on_cells(grid);
  EdgeData field_on_edges(grid);

  //TODO: calculate gradient of divergence

  // write to file
  // for(auto& wind : winds) {
  //   std::ofstream ofile;
  //   ofile.open (filename + "_" + std::to_string(step) + ".vtk");
  //   ofile << grid.toVtk() << winds[step].toVtk() << std::endl;
  //   ofile.close();
  //   step++;
  // }
  return 0;
}
