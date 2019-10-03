//cell diffusion example using MPI to explore compute domain
//
//  simply run with mpirun -np <N> cellDiffusionMPI <nx> <ny> <n_timestep> <outfname>
//


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

#include "eckit/types/Types.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/library/Library.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/mpi/mpi.h"

#include "eckit/config/Resource.h"
#include "eckit/runtime/Tool.h"

#include <stdio.h>

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

using namespace atlas;

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

static std::vector<std::vector<int>> buildCellConnectivitySlow(const atlas::Mesh& mesh) {
  std::vector<std::vector<int>> cellConn;
  for(int jCell = 0, size = mesh.cells().size(); jCell < size; ++jCell) {
    cellConn.push_back(cellNeighborOfCells(jCell, mesh, mesh.cells().edge_connectivity()));
  }
  return cellConn;
}

static std::vector<std::vector<int>> buildCellConnectivityFast(const atlas::Mesh& mesh) {
  std::vector<std::vector<int>> cellCellConn;
  
  const mesh::HybridElements::Connectivity& cellEdgeConnectivity = mesh.cells().edge_connectivity();
  const mesh::HybridElements::Connectivity& edgeCellConnectivity = mesh.edges().cell_connectivity();

  for(int idxCell = 0, size = mesh.cells().size(); idxCell < size; ++idxCell) {
    const int nbEdgesPerElem = cellEdgeConnectivity.cols(idxCell);
    std::vector<int> cellNeighbors;
    for (int iterEdge = 0; iterEdge < nbEdgesPerElem; iterEdge++) {
      const int idxEdge = cellEdgeConnectivity(idxCell, iterEdge);
      const int nbCellsPerEdge = edgeCellConnectivity.cols(idxEdge);
      for (int iterCell = 0; iterCell < nbCellsPerEdge; iterCell++) {
        const int idxCellNeighbor = edgeCellConnectivity(idxEdge, iterCell);
        if (idxCellNeighbor != idxCell && idxCellNeighbor != edgeCellConnectivity.missing_value()) {
          cellNeighbors.push_back(idxCellNeighbor);
        }
      }
    }
    cellCellConn.push_back(cellNeighbors);
  }

  return cellCellConn;
}

void buildCellConnectivityFastAtlas(atlas::Mesh& mesh) { 
  const mesh::HybridElements::Connectivity& cellEdgeConnectivity = mesh.cells().edge_connectivity();
  const mesh::HybridElements::Connectivity& edgeCellConnectivity = mesh.edges().cell_connectivity();
  mesh::HybridElements::Connectivity&       cellCellConnectivity = mesh.cells().cell_connectivity();

  for(int idxCell = 0, size = mesh.cells().size(); idxCell < size; ++idxCell) {
    const int nbEdgesPerElem = cellEdgeConnectivity.cols(idxCell);
    std::vector<int> cellNeighbors;

    for (int iterEdge = 0; iterEdge < nbEdgesPerElem; iterEdge++) {
      const int idxEdge = cellEdgeConnectivity(idxCell, iterEdge);
      const int nbCellsPerEdge = edgeCellConnectivity.cols(idxEdge);
      for (int iterCell = 0; iterCell < nbCellsPerEdge; iterCell++) {
        const int idxCellNeighbor = edgeCellConnectivity(idxEdge, iterCell);
        if (idxCellNeighbor != idxCell && idxCellNeighbor != edgeCellConnectivity.missing_value()) {
          cellNeighbors.push_back(idxCellNeighbor);
        }
      }
    }

    std::vector<int> initData(cellNeighbors.size(), cellCellConnectivity.missing_value());
    cellCellConnectivity.add(1, cellNeighbors.size(), initData.data());
    for (int copyIter = 0; copyIter < cellNeighbors.size(); copyIter++) {
      cellCellConnectivity.set(idxCell, copyIter, cellNeighbors[copyIter]);
    }
  }
}

static void compareCellConns(std::vector<std::vector<int>> left, std::vector<std::vector<int>> right) {
  assert(left.size() == right.size());
  for (int i = 0; i < left.size(); i++) {
    assert(left[i].size() == right[i].size());
    for (int j = 0; j < left[i].size(); j++) {
      assert(left[i][j] == right[i][j]);
    }
  }
}

static void verifyAtlasCellConns(const atlas::Mesh& mesh) { 
  auto refCellCellConnectivity = buildCellConnectivityFast(mesh);
  const mesh::HybridElements::Connectivity& cellCellConnectivity = mesh.cells().cell_connectivity();
  for(int jCell = 0, size = mesh.cells().size(); jCell < size; ++jCell) {
    int numNbh = cellCellConnectivity.cols(jCell);
    for(int nbhCelliter = 0; nbhCelliter < numNbh; nbhCelliter++) {
      int nbhCell = cellCellConnectivity(jCell, nbhCelliter);
      assert(refCellCellConnectivity[jCell][nbhCelliter] == nbhCell);
    }
  }
}

class myTool : public eckit::Tool {
  void run() override;

  int x_ = 0, y_ = 0;
  int timeSteps_ = 0;
  std::string filename_;

public:
  myTool(int argc, char** argv) : eckit::Tool(argc, argv) {
    x_ = atoi(argv[1]);
    y_ = atoi(argv[2]);
    timeSteps_ = atoi(argv[3]);
    filename_ = std::string(argv[4]);
  }
};

void myTool::run() {
  //mesh generation
  std::string x_s = std::to_string(x_);
  std::string y_s = std::to_string(y_);
  atlas::StructuredGrid grid = atlas::Grid("L" + x_s + "x" + y_s);
  MeshGenerator meshgenerator("structured");
  Mesh mesh = meshgenerator.generate(grid);

  //gmesh for output
  atlas::output::Gmsh gmesh(filename_ + ".msh");
  gmesh.write(mesh);  //write complete mesh in the beginning for easy visualization

  //vertical levels
  int nb_levels = 1;

  //halos (1 cell for cell to cell diffusion)
  const int halo = 1;
  functionspace::CellColumns cells_fs(mesh, option::halo(halo));

  // Cell to cell connectivity
  //      IMPORTANT: build this only AFTER halo has been set!
  atlas::mesh::actions::build_edges(mesh);
  atlas::mesh::actions::build_node_to_edge_connectivity(mesh);
  atlas::mesh::actions::build_element_to_edge_connectivity(mesh);

  //two temp fields for double buffering
  Field in_field = cells_fs.createField<double>(option::name("temp_in") | option::levels(nb_levels));
  Field out_field = cells_fs.createField<double>(option::name("temp_out") | option::levels(nb_levels));

  // Initialize everything with a nice gaussian
  const double rpi = 2.0 * asin(1.0);
  const double deg2rad = rpi / 180.;
  double width = 0.1;
  auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());

  auto temp_in = array::make_view<double, 2>(in_field);
  auto temp_out = array::make_view<double, 2>(out_field);

  const auto& nodeCellConn = mesh.cells().node_connectivity();
  for(int jCell = 0, size = mesh.cells().size(); jCell < size; ++jCell) {
    double iCoordinate = lonlat(nodeCellConn(jCell, 0), 0) * deg2rad;
    double jCoordinate = lonlat(nodeCellConn(jCell, 0), 1) * deg2rad + 1;
    double val = exp(-width * (pow(iCoordinate - 6.0 / 2.0, 2) + pow(jCoordinate - 2.0 / 2.0, 2)));
    temp_in(jCell, 0) = val;   
  }
  
  // precompute cell to cell connectivity
  // compareCellConns(buildCellConnectivityFast(mesh), buildCellConnectivitySlow(mesh));  //checks out!
  buildCellConnectivityFastAtlas(mesh);
  // verifyAtlasCellConns(mesh);  //checks out too!

  // time stepping
  double diffusionCoefficient = 1e-1;
  for(int step = 0; step < timeSteps_; step++) {
    // exchange halos
    cells_fs.haloExchange(in_field);

    // write output data
    in_field.metadata().set("step", step);
    gmesh.write(in_field);

    // do diffusion
     mesh::HybridElements::Connectivity& cellCellConnectivity = mesh.cells().cell_connectivity();
    for(int jCell = 0, size = mesh.cells().size(); jCell < size; ++jCell) {      
      double sum = 0.;
      int numNbh = cellCellConnectivity.cols(jCell);
      for(int nbhCelliter = 0; nbhCelliter < numNbh; nbhCelliter++) {
        sum += temp_in(cellCellConnectivity(jCell, nbhCelliter), 0);
      }
      temp_out(jCell, 0) =
          temp_in(jCell, 0) + diffusionCoefficient * (sum - numNbh * temp_in(jCell, 0));
    }

    // in = out for the next step
    std::swap(temp_in, temp_out);
  }
}

int main(int argc, char* argv[]) {
  if(argc != 5) {
    std::cout << "intended use is\n"
              << argv[0] << " <X> <Y> <# Timesteps> <outputfile>" << std::endl;
    return -1;
  }

  myTool tool(argc, argv);
  tool.start();
}