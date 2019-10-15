#include "grid.h"
#include <cmath>
#include <iostream>

int main(int argc, char const* argv[]) {
  std::cout << "hello" << std::endl;
  int x = 20;
  Grid g(x, 2);
  double width = 0.01;
  for(int i = 0; i < x; ++i) {
    for(int j = 0; j < x; ++j) {
      g(i, j) = exp(-width * (pow(i - double(x) / 2.0, 2) + pow(j - double(x) / 2.0, 2)));
    }
  }
  GridPrinter p("halo.vtk");
  p.printWithHalo(g);
  GridPrinter p2("nohalo.vtk");
  p2.printInterior(g);
  return 0;
}