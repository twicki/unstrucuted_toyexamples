#ifndef STRUCTURED_LIB_GRID_H
#define STRUCTURED_LIB_GRID_H
#include <string>
#include <vector>

class Grid;

struct GridData {
  GridData(int size) : data(size) {}
  std::vector<double> data;
  double operator()(int i, int j, int stride) const;
  double& operator()(int i, int j, int stride);
};

class GridPrinter {
  std::string filename_;

public:
  GridPrinter(std::string filename) : filename_(filename) {}
  void printInterior(Grid& g);
  void printWithHalo(Grid& g);
};
std::string printHaloData(Grid g, int nCells);

class Grid {

  int iSize;
  int jSize;
  int haloSize;
  GridData data;

public:
  Grid(int horizontal, int halo = 3)
      : iSize(horizontal), jSize(horizontal), haloSize(halo),
        data((2 * halo + horizontal) * (2 * halo + horizontal)) {}
  double operator()(int i, int j) const;
  double& operator()(int i, int j);

  friend void GridPrinter::printInterior(Grid& g);
  friend void GridPrinter::printWithHalo(Grid& g);
  friend std::string printHaloData(Grid g, int nCells);
};

#endif // STRUCTURED_LIB_GRID_H