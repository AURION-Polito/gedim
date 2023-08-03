#ifndef __MESHMATRICES_2D_4CELLS_MOCK_H
#define __MESHMATRICES_2D_4CELLS_MOCK_H

#include "MeshMatrices.hpp"

using namespace std;

namespace GedimUnitTesting
{
  /// \brief Mesh generated with Triangle with parameter 0.25
  class MeshMatrices_2D_4Cells_Mock final
  {
    public:
      Gedim::MeshMatrices Mesh;

      MeshMatrices_2D_4Cells_Mock() {
        Mesh.Dimension = 2;
        Mesh.NumberCell0D = 5;
        Mesh.Cell0DCoordinates = { 0,0,0,
                                   1,0,0,
                                   1,1,0,
                                   0,1,0,
                                   0.5,0.5,0 };
        Mesh.Cell0DMarkers = { 1,2,3,4,0 };
        Mesh.ActiveCell0D = { 1,1,1,1,1 };
        Mesh.NumberCell0DNeighbourCell1D = { 0,0,0,0,0,0 };
        Mesh.Cell0DNeighbourCell1Ds = { };
        Mesh.NumberCell0DNeighbourCell2D = { 0,0,0,0,0,0 };
        Mesh.Cell0DNeighbourCell2Ds = { };
        Mesh.NumberCell1D = 8;
        Mesh.Cell1DVertices = { 1,2,
                                2,4,
                                4,1,
                                3,0,
                                0,4,
                                4,3,
                                2,3,
                                0,1 };
        Mesh.Cell1DMarkers = { 6,0,0,8,0,0,7,5 };
        Mesh.ActiveCell1D = { 1,1,1,1,1,1,1,1 };
        Mesh.Cell1DOriginalCell1Ds.resize(Mesh.NumberCell1D, std::numeric_limits<unsigned int>::max());
        Mesh.NumberCell1DNeighbourCell2D = { 0,2,4,6,8,10,12,14,16 };
        Mesh.Cell1DNeighbourCell2Ds = { 4,0,2,0,3,0,4,1,3,1,2,1,4,2,4,3 };
        std::replace(Mesh.Cell1DNeighbourCell2Ds.begin(), Mesh.Cell1DNeighbourCell2Ds.end(), static_cast<unsigned int>(4), std::numeric_limits<unsigned int>::max());
        Mesh.NumberCell2D = 4;
        Mesh.NumberCell2DVertices = { 0,3,6,9,12 };
        Mesh.NumberCell2DEdges = { 0,3,6,9,12 };
        Mesh.Cell2DVertices = { 1,2,4,
                                3,0,4,
                                4,2,3,
                                0,1,4 };
        Mesh.Cell2DEdges = { 0,1,2,
                             3,4,5,
                             1,6,5,
                             7,2,4 };
        Mesh.NumberCell2DSubdivision.resize(Mesh.NumberCell2D + 1, 0);
        Mesh.Cell2DMarkers = { 0,0,0,0 };
        Mesh.Cell2DOriginalCell2Ds.resize(Mesh.NumberCell2D, std::numeric_limits<unsigned int>::max());
        Mesh.ActiveCell2D = { 1,1,1,1 };
        Mesh.NumberCell2DNeighbourCell3D = { 0,0,0,0,0 };
        Mesh.Cell2DNeighbourCell3Ds = { };
      }
  };
}

#endif // __MESHMATRICES_2D_4CELLS_MOCK_H
