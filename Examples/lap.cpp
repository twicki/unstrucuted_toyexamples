#include "unstructured_toylib.h"
#include <iostream>

int main(int argc, char const* argv[]) {
  /* code */
  std::cout << "hello" << std::endl;
  Grid grid;
  grid.toVtk();
  return 0;
}