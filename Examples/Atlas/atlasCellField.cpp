#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/UnstructuredGrid.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/meshgenerator.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/detail/DelaunayMeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"
#include "eckit/types/FloatCompare.h"
#include <array/ArrayShape.h>
#include <array/DataType.h>
#include <array_fwd.h>
#include <iostream>
#include <sstream>

static std::vector<int>
cellNeighborOfCells(int idx, const atlas::Mesh& mesh,
                    const atlas::mesh::HybridElements::Connectivity& connectitvity) {
  std::vector<int> neighs;
  for(int n = 0; n < connectitvity.cols(idx); ++n) {
    int initialEdge = connectitvity(idx, n);
    for(int c1 = 0; c1 < mesh.cells().size(); ++c1) {
      for(int n1 = 0; n1 < connectitvity.cols(c1); ++n1) {
        int compareEdge = connectitvity(c1, n1);
        if(initialEdge == compareEdge && c1 != idx) {
          neighs.push_back(c1);
        }
      }
    }
  }
  return neighs;
}

int main(int argc, char const* argv[]) {
  if(argc != 5) {
    std::cout << "intended use is\n"
              << argv[0] << " <X> <Y> <# Timesteps> <outputfile base>" << std::endl;
    return -1;
  }
  int x = atoi(argv[1]);
  int y = atoi(argv[2]);
  int timeSteps = atoi(argv[3]);
  std::string filename(argv[4]);
  ////////////////////// Regular Grid setup
  // Grid setup: 10 x 10 lat lon points
  atlas::StructuredGrid structuredGrid = atlas::Grid("L3x4");
  std::cout << "grid.npts() = " << structuredGrid.size() << std::endl;

  // Mesh setup to this grid (quads)
  atlas::StructuredMeshGenerator generator;
  auto mesh = generator.generate(structuredGrid);
  atlas::output::Gmsh gmesh("mymesh.msh");
  gmesh.write(mesh);

  // Field
  int nb_levels = 1;
  atlas::Field field{"tmp", atlas::array::make_datatype<double>(),
                     atlas::array::make_shape(mesh.cells().size(), nb_levels)};
  auto temp2 = atlas::array::make_view<double, 2>(field);

  // Cell to cell connectivity
  atlas::mesh::actions::build_edges(mesh);
  atlas::mesh::actions::build_node_to_edge_connectivity(mesh);
  const auto& conn = mesh.cells().edge_connectivity();
  for(int c = 0; c < mesh.cells().size(); ++c) {
    std::cout << "neighbours of cell " << c << ":" << std::endl;
    auto neighs = cellNeighborOfCells(c, mesh, conn);
    for(auto neigh : neighs) {
      std::cout << "\t"
                << "cell: " << neigh << std::endl;
    }
  }

  // Initialize everything with a nice gaussian
  const double rpi = 2.0 * asin(1.0);
  const double deg2rad = rpi / 180.;
  double width = 0.1;
  auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
  const auto& nodeCellConn = mesh.cells().node_connectivity();
  for(int jCell = 0, size = mesh.cells().size(); jCell < size; ++jCell) {
    double iCoordinate = lonlat(nodeCellConn(jCell, 0), 0) * deg2rad;
    double jCoordinate = lonlat(nodeCellConn(jCell, 0), 1) * deg2rad + 1;
    temp2(jCell, 0) =
        exp(-width * (pow(iCoordinate - 6.0 / 2.0, 2) + pow(jCoordinate - 2.0 / 2.0, 2)));
  }
  //   field.metadata().set("step", 0);
  gmesh.write(field);

  return 0;
}