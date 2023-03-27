#ifndef __MESHMATRICES_2D_2CELLS_MOCK_H
#define __MESHMATRICES_2D_2CELLS_MOCK_H

#include "MeshMatrices.hpp"

using namespace std;

namespace GedimUnitTesting
{
  /// \brief Mesh generated with Triangle with parameter 0.5
  class MeshMatrices_2D_2Cells_Mock final
  {
    public:
      Gedim::MeshMatrices Mesh;

      MeshMatrices_2D_2Cells_Mock() {
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
        Mesh.NumberCell1D = 5;
        Mesh.Cell1DVertices = { 3,0,
                                0,1,
                                1,3,
                                1,2,
                                2,3 };
        Mesh.Cell1DMarkers = { 8,5,0,6,7 };
        Mesh.ActiveCell1D = { 1,1,1,1,1 };
        Mesh.Cell1DOriginalCell1Ds = { 5,5,5,5,5 };
        Mesh.Cell1DAdjacency.resize(4, 4);
        Mesh.Cell1DAdjacency.reserve(5);
        Mesh.Cell1DAdjacency.insert(3, 0) = 1;
        Mesh.Cell1DAdjacency.insert(0, 1) = 2;
        Mesh.Cell1DAdjacency.insert(1, 3) = 3;
        Mesh.Cell1DAdjacency.insert(1, 2) = 4;
        Mesh.Cell1DAdjacency.insert(2, 3) = 5;
        Mesh.NumberCell1DNeighbourCell2D = { 0,2,4,6,8,10 };
        Mesh.Cell1DNeighbourCell2Ds = { 2,0,2,0,1,0,2,1,2,1 };
        Mesh.NumberCell2D = 2;
        Mesh.NumberCell2DVertices = { 0, 3, 6 };
        Mesh.NumberCell2DEdges = { 0, 3, 6  };
        Mesh.Cell2DVertices = { 3,0,1,
                                1,2,3 };
        Mesh.Cell2DEdges = { 0,1,2,
                             3,4,2 };
        Mesh.Cell2DMarkers = { 0,0 };
        Mesh.ActiveCell2D = { 1,1 };
      }
  };
}

#endif // __MESHMATRICES_2D_2CELLS_MOCK_H
