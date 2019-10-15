#include "grid.h"

double GridData::operator()(int i, int j, int stride) { return data[i + j * stride]; }

double Grid::operator()(int i, int j) {
  return data(i + haloSize, j + haloSize, jSize + 2 * haloSize);
}