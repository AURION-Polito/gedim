#ifndef __MeshMatrices_3D_1Cells_Mock_H
#define __MeshMatrices_3D_1Cells_Mock_H

#include "MeshMatrices.hpp"

using namespace std;

namespace GedimUnitTesting
{
  /// \brief Mesh generated with Triangle with parameter 0.5
  class MeshMatrices_3D_1Cells_Mock final
  {
    public:
      Gedim::MeshMatrices Mesh;

      MeshMatrices_3D_1Cells_Mock() {
        Mesh.Dimension = 3;
        Mesh.NumberCell0D = 8;
        Mesh.Cell0DCoordinates = { 0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,1.0000000000000000e+00,1.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,1.0000000000000000e+00 };
        Mesh.Cell0DMarkers = { 1,2,3,4,5,6,7,8 };
        Mesh.ActiveCell0D = { 1,1,1,1,1,1,1,1 };
        Mesh.UpdatedCell0Ds = {};
        Mesh.NumberCell0DNeighbourCell1D = { 0,0,0,0,0,0,0,0,0 };
        Mesh.Cell0DNeighbourCell1Ds = {};
        Mesh.NumberCell0DNeighbourCell2D = { 0,0,0,0,0,0,0,0,0 };
        Mesh.Cell0DNeighbourCell2Ds = {};
        Mesh.Cell0DDoublePropertyIds = {};
        Mesh.Cell0DDoublePropertyIndices = {};
        Mesh.Cell0DDoublePropertySizes = {};
        Mesh.Cell0DDoublePropertyValues = {};
        Mesh.NumberCell1D = 12;
        Mesh.Cell1DVertices = { 0,1,1,2,2,3,3,0,4,5,5,6,6,7,7,4,0,4,1,5,2,6,3,7 };
        Mesh.Cell1DAdjacency.resize(8, 8);
        Mesh.Cell1DAdjacency.reserve(12);
        Mesh.Cell1DAdjacency.insert(3, 0) = 4;
        Mesh.Cell1DAdjacency.insert(0, 1) = 1;
        Mesh.Cell1DAdjacency.insert(1, 2) = 2;
        Mesh.Cell1DAdjacency.insert(2, 3) = 3;
        Mesh.Cell1DAdjacency.insert(0, 4) = 9;
        Mesh.Cell1DAdjacency.insert(7, 4) = 8;
        Mesh.Cell1DAdjacency.insert(1, 5) = 10;
        Mesh.Cell1DAdjacency.insert(4, 5) = 5;
        Mesh.Cell1DAdjacency.insert(2, 6) = 11;
        Mesh.Cell1DAdjacency.insert(5, 6) = 6;
        Mesh.Cell1DAdjacency.insert(3, 7) = 12;
        Mesh.Cell1DAdjacency.insert(6, 7) = 7;
        Mesh.Cell1DAdjacency.makeCompressed();
        Mesh.Cell1DMarkers = { 9,10,11,12,13,14,15,16,17,18,19,20 };
        Mesh.ActiveCell1D = { 1,1,1,1,1,1,1,1,1,1,1,1 };
        Mesh.UpdatedCell1Ds = {};
        Mesh.NumberCell1DNeighbourCell2D = { 0,0,0,0,0,0,0,0,0,0,0,0,0 };
        Mesh.Cell1DNeighbourCell2Ds = {};
        Mesh.Cell1DDoublePropertyIds = {};
        Mesh.Cell1DDoublePropertyIndices = {};
        Mesh.Cell1DDoublePropertySizes = {};
        Mesh.Cell1DDoublePropertyValues = {};
        Mesh.NumberCell2D = 6;
        Mesh.NumberCell2DVertices = { 0,4,8,12,16,20,24 };
        Mesh.Cell2DVertices = { 0,1,2,3,4,5,6,7,0,3,7,4,1,2,6,5,0,1,5,4,3,2,6,7 };
        Mesh.NumberCell2DEdges = { 0,4,8,12,16,20,24 };
        Mesh.Cell2DEdges = { 0,1,2,3,4,5,6,7,3,11,7,8,1,10,5,9,0,9,4,8,2,10,6,11 };
        Mesh.NumberCell2DSubdivision = { 0,0,0,0,0,0,0 };
        Mesh.Cell2DSubdivision = {};
        Mesh.Cell2DMarkers = { 21,22,23,24,25,26 };
        Mesh.ActiveCell2D = { 1,1,1,1,1,1 };
        Mesh.Cell2DOriginalCell2Ds = { 6,6,6,6,6,6 };
        Mesh.UpdatedCell2Ds = {};
        Mesh.Cell2DDoublePropertyIds = {};
        Mesh.Cell2DDoublePropertyIndices = {};
        Mesh.Cell2DDoublePropertySizes = {};
        Mesh.Cell2DDoublePropertyValues = {};
        Mesh.NumberCell3D = 1;
        Mesh.NumberCell3DVertices = { 0,8 };
        Mesh.Cell3DVertices = { 0,1,2,3,4,5,6,7 };
        Mesh.NumberCell3DEdges = { 0,12 };
        Mesh.Cell3DEdges = { 0,1,2,3,4,5,6,7,8,9,10,11 };
        Mesh.NumberCell3DFaces = { 0,6 };
        Mesh.Cell3DFaces = { 0,1,2,3,4,5 };
        Mesh.Cell3DMarkers = { 0 };
        Mesh.ActiveCell3D = { 1 };
        Mesh.UpdatedCell3Ds = {};
        Mesh.Cell3DDoublePropertyIds = {};
        Mesh.Cell3DDoublePropertyIndices = {};
        Mesh.Cell3DDoublePropertySizes = {};
        Mesh.Cell3DDoublePropertyValues = {};
      }
  };
}

#endif // __MeshMatrices_3D_1Cells_Mock_H
