#ifndef __MeshMatrices_2D_1Cells_Mock_H
#define __MeshMatrices_2D_1Cells_Mock_H

#include "MeshMatrices.hpp"

using namespace std;

namespace GedimUnitTesting
{
  /// \brief Mesh generated with Triangle with parameter 0.5
  class MeshMatrices_2D_1Cells_Mock final
  {
    public:
      Gedim::MeshMatrices Mesh;

      MeshMatrices_2D_1Cells_Mock() {
        Mesh.Dimension = 2;
        Mesh.NumberCell0D = 4;
        Mesh.Cell0DCoordinates = { 0.0, 0.0, 0.0,
                                   1.0, 0.0, 0.0,
                                   1.0, 1.0, 0.0,
                                   0.0, 1.0, 0.0 };
        Mesh.Cell0DMarkers = { 1,2,3,4 };
        Mesh.ActiveCell0D = { 1,1,1,1 };
        Mesh.NumberCell0DNeighbourCell1D = { 0,0,0,0,0 };
        Mesh.Cell0DNeighbourCell1Ds = { };
        Mesh.NumberCell0DNeighbourCell2D = { 0,0,0,0,0 };
        Mesh.Cell0DNeighbourCell2Ds = { };
        Mesh.NumberCell1D = 4;
        Mesh.Cell1DVertices = { 3,0,
                                0,1,
                                1,2,
                                2,3 };
        Mesh.Cell1DMarkers = { 8,5,6,7 };
        Mesh.ActiveCell1D = { 1,1,1,1 };
        Mesh.Cell1DAdjacency.resize(4, 4);
        Mesh.Cell1DAdjacency.reserve(4);
        Mesh.Cell1DAdjacency.insert(3, 0) = 1;
        Mesh.Cell1DAdjacency.insert(0, 1) = 2;
        Mesh.Cell1DAdjacency.insert(1, 2) = 3;
        Mesh.Cell1DAdjacency.insert(2, 3) = 4;
        Mesh.NumberCell1DNeighbourCell2D = { 0,2,4,6,8 };
        Mesh.Cell1DNeighbourCell2Ds = { 1,0,1,0,1,0,1,0 };
        Mesh.NumberCell2D = 1;
        Mesh.NumberCell2DVertices = { 0, 4 };
        Mesh.NumberCell2DEdges = { 0, 4  };
        Mesh.Cell2DVertices = { 3,0,1,2 };
        Mesh.Cell2DEdges = { 0,1,2,3 };
        Mesh.Cell2DMarkers = { 0 };
        Mesh.ActiveCell2D = { 1 };
      }
  };
}

#endif // __MeshMatrices_2D_1Cells_Mock_H
