#include "grid.h"
#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>
#include <string>

double GridData::operator()(int i, int j, int stride) const { return data[i + j * stride]; }

double Grid::operator()(int i, int j) const {
  return data(i + haloSize, j + haloSize, jSize + 2 * haloSize);
}

double& GridData::operator()(int i, int j, int stride) { return data[i + j * stride]; }

double& Grid::operator()(int i, int j) {
  return data(i + haloSize, j + haloSize, jSize + 2 * haloSize);
}

std::string printInteriorData(Grid& g, int iSize, int jSize) {
  std::ostringstream os;
  for(int i = 0; i < iSize - 1; ++i) {
    for(int j = 0; j < jSize - 1; ++j) {
      os << g(i, j) << std::endl;
    }
  }
  return os.str();
}

std::string printHaloData(Grid g, int nCells) {
  std::ostringstream os;
  for(int i = 0; i < nCells; ++i) {
    os << std::to_string(g.data.data[i]) << std::endl;
  }
  return os.str();
}

static void printImpl(Grid& g, int iSize, int jSize, std::string filename, bool printHalo = false) {
  std::ostringstream os;
  int nPoints = iSize * jSize;
  os << "# vtk DataFile Version 3.0\n2D scalar data\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
     << nPoints << " double\n";
  for(int i = 0, iMax = iSize; i < iMax; ++i) {
    for(int j = 0, jMax = jSize; j < jMax; ++j) {
      os << i << " " << j << " 0" << std::endl;
    }
  }

  int nCells = (iSize - 1) * (jSize - 1);
  os << "CELLS " << nCells << " " << 5 * nCells << std::endl;
  for(int i = 0; i < iSize - 1; ++i) {
    for(int j = 0; j < jSize - 1; ++j) {
      os << "4 " << i + j * iSize << " " << (i + 1) + j * iSize << " " << i + (j + 1) * iSize << " "
         << (i + 1) + (j + 1) * iSize << std::endl;
    }
  }
  os << "CELL_TYPES " << nCells << std::endl;
  for(int i = 0; i < nCells; ++i) {
    os << "8" << std::endl;
  }
  os << "CELL_DATA " << nCells << "\nSCALARS temperature double 1\nLOOKUP_TABLE default\n";
  if(printHalo)
    os << printHaloData(g, nCells);
  else
    os << printInteriorData(g, iSize, jSize);

  std::ofstream outfile;
  outfile.open(filename, std::ios::trunc);
  outfile << os.str() << std::endl;
  outfile.close();
}

void GridPrinter::printWithHalo(Grid& g) {
  int iSize = (g.iSize + 2 * g.haloSize + 1);
  int jSize = (g.jSize + 2 * g.haloSize + 1);
  printImpl(g, iSize, jSize, filename_, true);
}

void GridPrinter::printInterior(Grid& g) {
  int iSize = (g.iSize + 1);
  int jSize = (g.jSize + 1);
  printImpl(g, iSize, jSize, filename_);
}
