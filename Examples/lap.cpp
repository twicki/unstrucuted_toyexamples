#include "unstructured_toylib.h"
#include <iostream>

int main(int argc, char const* argv[]) {
  /* code */
  // std::cout << "hello" << std::endl;
  int x = atoi(argv[1]);
  int y = atoi(argv[2]);
  Grid grid(x, y);

  std::vector<Data> winds;
  winds.emplace_back(grid);
  winds.emplace_back(grid);

  // winds[0].initGauss();

  for(auto triangle : grid.getTriangles()) {
    if(triangle->getId() == 200) {
      winds[0].triangleData_[triangle] = 1;
    } else {
      winds[0].triangleData_[triangle] = 0;
    }
  }

  // for(auto triangle : grid.getTriangles()) {
  //   auto neighbours = grid.adjacencyList(*triangle);
  //   winds[1].triangleData_[triangle] = -3.0 * winds[0].triangleData_[triangle];
  //   for(auto neigh : neighbours)
  //     winds[1].triangleData_[triangle] += winds[0].triangleData_[neigh];
  // }

  std::cout << grid.toVtk() << winds[0].toVtk() << std::endl;
  return 0;
}
