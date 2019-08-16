#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/UnstructuredGrid.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
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
  ////////////////////// Regular Grid setup
  // Grid setup: 10 x 10 lat lon points
  atlas::StructuredGrid structuredGrid = atlas::Grid("L9x10");
  // Mesh setup to this grid (quads)
  atlas::StructuredMeshGenerator generator;
  auto mesh = generator.generate(structuredGrid);
  atlas::output::Gmsh gmesh("mymesh.msh");
  gmesh.write(mesh);

  // Field based on the mesh
  atlas::functionspace::NodeColumns fs_nodes(mesh, atlas::option::levels(1));
  auto field_temp2 = fs_nodes.createField<double>(atlas::option::name("newtemp"));

  // covert radial coords to x-y
  const double rpi = 2.0 * asin(1.0);
  const double deg2rad = rpi / 180.;

  auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
  auto temp2 = atlas::array::make_view<double, 2>(field_temp2);
  for(int jNode = 0, size = mesh.nodes().size(); jNode < size; ++jNode) {
    int iCoordinate = lonlat(jNode, 0) * deg2rad;
    int jCoordinate = lonlat(jNode, 1) * deg2rad + 1;
    // std::cout << iCoordinate << ", " << jCoordinate << std::endl;
    temp2(jNode, 0) = iCoordinate + 10 * jCoordinate;
    auto& node = mesh.nodes();
    auto& connectivity = node.cell_connectivity();
  }
  gmesh.write(field_temp2);
}