#include "unstructured_toylib.h"
#include <iostream>

int main(int argc, char const* argv[]) {
  /* code */
  std::cout << "hello" << std::endl;
  Grid grid;
  std::cout << grid.toVtk() << std::endl;
  return 0;
}