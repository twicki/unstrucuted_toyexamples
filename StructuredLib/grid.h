#ifndef STRUCTURED_LIB_GRID_H
#define STRUCTURED_LIB_GRID_H
#include <vector>
struct GridData {
  GridData(int size) : data(size) {}
  std::vector<double> data;
  double operator()(int i, int j, int stride);
};

class Grid {

  int iSize;
  int jSize;
  int haloSize;
  GridData data;

public:
  Grid(int horizontal, int halo = 3)
      : iSize(horizontal), jSize(horizontal), haloSize(halo),
        data((2 * halo + horizontal) * (2 * halo + horizontal)) {}
  double operator()(int i, int j);
};

#endif // STRUCTURED_LIB_GRID_H