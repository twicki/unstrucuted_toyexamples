#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/UnstructuredGrid.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildNode2CellConnectivity.h"
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

int main() {
  int timeSteps = 10;
  ////////////////////// Regular Grid setup
  // Grid setup: 10 x 10 lat lon points
  atlas::StructuredGrid structuredGrid = atlas::Grid("L3x4");
  std::cout << "grid.npts() = " << structuredGrid.size() << std::endl;

  // Mesh setup to this grid (quads)
  atlas::StructuredMeshGenerator generator;
  auto mesh = generator.generate(structuredGrid);
  atlas::mesh::actions::BuildNode2CellConnectivity{mesh}();
  const auto& node_to_cell = mesh.nodes().cell_connectivity();
  for(int c = 0; c < mesh.nodes().size(); ++c) {
    std::cout << "neighbours of cell " << c << ":" << std::endl;
    for(int n = 0; n < node_to_cell.cols(c); ++n) {
      int node = node_to_cell(c, n);
      std::cout << "\t"
                << "cell: " << node << std::endl;
    }
  }
  atlas::output::Gmsh gmesh("mymesh.msh");
  gmesh.write(mesh);

  std::cout << "connectivity_table_size(): " << mesh.nodes().cell_connectivity().size()
            << std::endl;
  // Field
  atlas::functionspace::NodeColumns fs_nodes(mesh, atlas::option::levels(1));
  auto field_temp2 = fs_nodes.createField<double>(atlas::option::name("temperature"));
  auto temp2 = atlas::array::make_view<double, 2>(field_temp2);

  // Initialize everything with a nice gaussian
  auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
  const double rpi = 2.0 * asin(1.0);
  const double deg2rad = rpi / 180.;
  double width = 0.1;
  for(int jNode = 0, size = mesh.nodes().size(); jNode < size; ++jNode) {
    double iCoordinate = lonlat(jNode, 0) * deg2rad;
    double jCoordinate = lonlat(jNode, 1) * deg2rad + 1;
    std::cout << "node: " << jNode << " (" << lonlat(jNode, 0) << "," << lonlat(jNode, 1) << ")"
              << std::endl;
    temp2(jNode, 0) =
        exp(-width * (pow(iCoordinate - 6.0 / 2.0, 2) + pow(jCoordinate - 2.0 / 2.0, 2)));
  }
  field_temp2.metadata().set("step", 0);
  gmesh.write(field_temp2);
  return 0;
}