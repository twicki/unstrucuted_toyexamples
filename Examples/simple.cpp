#include "grid.h"
#include <cmath>
#include <iostream>

int main(int argc, char const* argv[]) {
  std::cout << "hello" << std::endl;
  int horizontal = 20;
  int nHalo = 2;
  Grid g(horizontal, nHalo);
  double width = 0.01;
  for(int i = 0; i < horizontal; ++i) {
    for(int j = 0; j < horizontal; ++j) {
      g(i, j) = exp(-width *
                    (pow(i - double(horizontal) / 2.0, 2) + pow(j - double(horizontal) / 2.0, 2)));
    }
  }
  GridPrinter p("halo.vtk");
  p.printWithHalo(g);
  GridPrinter p2("nohalo.vtk");
  p2.printInterior(g);
  return 0;
}