#ifndef __MESHMATRICES_2D_26CELLS_MOCK_H
#define __MESHMATRICES_2D_26CELLS_MOCK_H

#include "MeshMatrices.hpp"

namespace GedimUnitTesting
{
  /// \brief Mesh generated with Triangle with parameter 0.033
  class MeshMatrices_2D_26Cells_Mock final
  {
    public:
      Gedim::MeshMatrices Mesh;

      MeshMatrices_2D_26Cells_Mock() {
        Mesh.Dimension = 2;
        Mesh.NumberCell0D = 26;
        Mesh.Cell0DCoordinates = { 0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,5.0000000000000000e-01,5.0000000000000000e-01,0.0000000000000000e+00,0.0000000000000000e+00,5.0000000000000000e-01,0.0000000000000000e+00,5.0000000000000000e-01,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,5.0000000000000000e-01,0.0000000000000000e+00,5.0000000000000000e-01,1.0000000000000000e+00,0.0000000000000000e+00,7.5000000000000000e-01,2.5000000000000000e-01,0.0000000000000000e+00,2.5000000000000000e-01,7.5000000000000000e-01,0.0000000000000000e+00,7.5000000000000000e-01,7.5000000000000000e-01,0.0000000000000000e+00,5.0000000000000000e-01,2.5000000000000000e-01,0.0000000000000000e+00,2.5000000000000000e-01,3.7500000000000000e-01,0.0000000000000000e+00,3.1250000000000000e-01,5.6250000000000000e-01,0.0000000000000000e+00,0.0000000000000000e+00,2.5000000000000000e-01,0.0000000000000000e+00,7.5000000000000000e-01,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,2.5000000000000000e-01,0.0000000000000000e+00,7.5000000000000000e-01,5.0000000000000000e-01,0.0000000000000000e+00,5.0000000000000000e-01,7.5000000000000000e-01,0.0000000000000000e+00,2.5000000000000000e-01,1.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,7.5000000000000000e-01,0.0000000000000000e+00,1.0000000000000000e+00,7.5000000000000000e-01,0.0000000000000000e+00,7.5000000000000000e-01,1.0000000000000000e+00,0.0000000000000000e+00,2.5000000000000000e-01,0.0000000000000000e+00,0.0000000000000000e+00,3.1250000000000000e-01,1.8750000000000000e-01,0.0000000000000000e+00 };
        Mesh.Cell0DMarkers = { 1,2,3,4,0,8,5,6,7,0,0,0,0,0,0,8,5,6,0,0,7,8,6,7,5,0 };
        Mesh.ActiveCell0D = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
        Mesh.UpdatedCell0Ds = {};
        Mesh.NumberCell0DNeighbourCell1D = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
        Mesh.Cell0DNeighbourCell1Ds = {};
        Mesh.NumberCell0DNeighbourCell2D = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
        Mesh.Cell0DNeighbourCell2Ds = {};
        Mesh.Cell0DDoublePropertyIds = {};
        Mesh.Cell0DDoublePropertyIndices = {};
        Mesh.Cell0DDoublePropertySizes = {};
        Mesh.Cell0DDoublePropertyValues = {};
        Mesh.NumberCell1D = 59;
        Mesh.Cell1DVertices = { 18,11,11,4,4,18,25,6,6,12,12,25,19,8,8,10,10,19,5,13,13,14,14,5,12,13,13,25,6,9,9,12,4,9,9,18,19,11,11,8,17,7,7,9,9,17,16,1,1,9,9,16,21,5,5,10,10,21,20,3,3,10,10,20,23,8,11,23,22,2,2,11,11,22,12,4,4,13,0,24,24,15,15,0,4,14,19,14,4,19,14,10,6,16,5,15,15,13,1,17,7,11,18,7,8,20,3,21,7,22,2,23,15,25,24,25,24,6 };
        Mesh.Cell1DAdjacency.resize(26, 26);
        Mesh.Cell1DAdjacency.reserve(59);
        Mesh.Cell1DAdjacency.insert(15, 0) = 42;
        Mesh.Cell1DAdjacency.insert(16, 1) = 24;
        Mesh.Cell1DAdjacency.insert(22, 2) = 35;
        Mesh.Cell1DAdjacency.insert(20, 3) = 30;
        Mesh.Cell1DAdjacency.insert(11, 4) = 2;
        Mesh.Cell1DAdjacency.insert(12, 4) = 38;
        Mesh.Cell1DAdjacency.insert(14, 5) = 12;
        Mesh.Cell1DAdjacency.insert(21, 5) = 27;
        Mesh.Cell1DAdjacency.insert(24, 6) = 59;
        Mesh.Cell1DAdjacency.insert(25, 6) = 4;
        Mesh.Cell1DAdjacency.insert(17, 7) = 21;
        Mesh.Cell1DAdjacency.insert(18, 7) = 52;
        Mesh.Cell1DAdjacency.insert(11, 8) = 20;
        Mesh.Cell1DAdjacency.insert(19, 8) = 7;
        Mesh.Cell1DAdjacency.insert(23, 8) = 33;
        Mesh.Cell1DAdjacency.insert(1, 9) = 25;
        Mesh.Cell1DAdjacency.insert(4, 9) = 17;
        Mesh.Cell1DAdjacency.insert(6, 9) = 15;
        Mesh.Cell1DAdjacency.insert(7, 9) = 22;
        Mesh.Cell1DAdjacency.insert(3, 10) = 31;
        Mesh.Cell1DAdjacency.insert(5, 10) = 28;
        Mesh.Cell1DAdjacency.insert(8, 10) = 8;
        Mesh.Cell1DAdjacency.insert(14, 10) = 46;
        Mesh.Cell1DAdjacency.insert(2, 11) = 36;
        Mesh.Cell1DAdjacency.insert(7, 11) = 51;
        Mesh.Cell1DAdjacency.insert(18, 11) = 1;
        Mesh.Cell1DAdjacency.insert(19, 11) = 19;
        Mesh.Cell1DAdjacency.insert(6, 12) = 5;
        Mesh.Cell1DAdjacency.insert(9, 12) = 16;
        Mesh.Cell1DAdjacency.insert(4, 13) = 39;
        Mesh.Cell1DAdjacency.insert(5, 13) = 10;
        Mesh.Cell1DAdjacency.insert(12, 13) = 13;
        Mesh.Cell1DAdjacency.insert(15, 13) = 49;
        Mesh.Cell1DAdjacency.insert(4, 14) = 43;
        Mesh.Cell1DAdjacency.insert(13, 14) = 11;
        Mesh.Cell1DAdjacency.insert(19, 14) = 44;
        Mesh.Cell1DAdjacency.insert(5, 15) = 48;
        Mesh.Cell1DAdjacency.insert(24, 15) = 41;
        Mesh.Cell1DAdjacency.insert(6, 16) = 47;
        Mesh.Cell1DAdjacency.insert(9, 16) = 26;
        Mesh.Cell1DAdjacency.insert(1, 17) = 50;
        Mesh.Cell1DAdjacency.insert(9, 17) = 23;
        Mesh.Cell1DAdjacency.insert(4, 18) = 3;
        Mesh.Cell1DAdjacency.insert(9, 18) = 18;
        Mesh.Cell1DAdjacency.insert(4, 19) = 45;
        Mesh.Cell1DAdjacency.insert(10, 19) = 9;
        Mesh.Cell1DAdjacency.insert(8, 20) = 53;
        Mesh.Cell1DAdjacency.insert(10, 20) = 32;
        Mesh.Cell1DAdjacency.insert(3, 21) = 54;
        Mesh.Cell1DAdjacency.insert(10, 21) = 29;
        Mesh.Cell1DAdjacency.insert(7, 22) = 55;
        Mesh.Cell1DAdjacency.insert(11, 22) = 37;
        Mesh.Cell1DAdjacency.insert(2, 23) = 56;
        Mesh.Cell1DAdjacency.insert(11, 23) = 34;
        Mesh.Cell1DAdjacency.insert(0, 24) = 40;
        Mesh.Cell1DAdjacency.insert(12, 25) = 6;
        Mesh.Cell1DAdjacency.insert(13, 25) = 14;
        Mesh.Cell1DAdjacency.insert(15, 25) = 57;
        Mesh.Cell1DAdjacency.insert(24, 25) = 58;
        Mesh.Cell1DAdjacency.makeCompressed();
        Mesh.Cell1DMarkers = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,0,0,5,0,0,8,0,0,7,0,0,7,0,6,0,0,0,0,5,0,8,0,0,0,0,5,8,0,6,0,0,7,8,6,7,0,0,5 };
        Mesh.ActiveCell1D = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
        Mesh.UpdatedCell1Ds = {};
        Mesh.NumberCell1DNeighbourCell2D = { 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118 };
        Mesh.Cell1DNeighbourCell2Ds = { 23,0,25,0,6,0,33,1,5,1,4,1,7,2,27,2,26,2,21,3,17,3,19,3,14,4,31,4,20,5,15,5,15,6,24,6,25,7,12,7,34,8,24,8,22,8,34,9,22,9,20,9,34,10,19,10,28,10,34,11,28,11,27,11,34,12,30,12,34,13,30,13,29,13,15,14,17,14,34,16,32,16,34,16,18,17,26,18,25,18,26,19,34,20,34,21,31,21,34,22,29,23,24,23,34,27,34,28,34,29,34,30,32,31,33,32,34,33 };
        Mesh.Cell1DDoublePropertyIds = {};
        Mesh.Cell1DDoublePropertyIndices = {};
        Mesh.Cell1DDoublePropertySizes = {};
        Mesh.Cell1DDoublePropertyValues = {};
        Mesh.NumberCell2D = 34;
        Mesh.NumberCell2DVertices = { 0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102 };
        Mesh.Cell2DVertices = { 18,11,4,25,6,12,19,8,10,5,13,14,25,12,13,6,9,12,4,9,18,19,11,8,17,7,9,16,1,9,21,5,10,20,3,10,23,8,11,22,2,11,13,12,4,9,4,12,0,24,15,14,13,4,19,14,4,5,14,10,9,6,16,13,5,15,9,1,17,7,11,18,9,7,18,4,11,19,14,19,10,10,8,20,10,3,21,11,7,22,11,2,23,15,25,13,24,25,15,6,25,24 };
        Mesh.NumberCell2DEdges = { 0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102 };
        Mesh.Cell2DEdges = { 0,1,2,3,4,5,6,7,8,9,10,11,5,12,13,14,15,4,16,17,2,18,19,6,20,21,22,23,24,25,26,27,28,29,30,31,32,19,33,34,35,36,12,37,38,16,37,15,39,40,41,10,38,42,43,42,44,11,45,27,14,46,25,9,47,48,24,49,22,50,0,51,21,51,17,1,18,44,43,8,45,7,52,31,30,53,28,50,54,36,35,55,33,56,13,48,57,56,40,3,57,58 };
        Mesh.NumberCell2DSubdivision = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
        Mesh.Cell2DSubdivision = {};
        Mesh.Cell2DMarkers = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
        Mesh.ActiveCell2D = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
        Mesh.Cell2DOriginalCell2Ds = { 34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34 };
        Mesh.UpdatedCell2Ds = {};
        Mesh.Cell2DDoublePropertyIds = {};
        Mesh.Cell2DDoublePropertyIndices = {};
        Mesh.Cell2DDoublePropertySizes = {};
        Mesh.Cell2DDoublePropertyValues = {};
        Mesh.NumberCell2DNeighbourCell3D = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
        Mesh.Cell2DNeighbourCell3Ds = {  };
        Mesh.NumberCell3D = 0;
        Mesh.NumberCell3DVertices = {};
        Mesh.Cell3DVertices = {};
        Mesh.NumberCell3DEdges = {};
        Mesh.Cell3DEdges = {};
        Mesh.NumberCell3DFaces = {};
        Mesh.Cell3DFaces = {};
        Mesh.Cell3DMarkers = {};
        Mesh.ActiveCell3D = {};
        Mesh.UpdatedCell3Ds = {};
        Mesh.Cell3DDoublePropertyIds = {};
        Mesh.Cell3DDoublePropertyIndices = {};
        Mesh.Cell3DDoublePropertySizes = {};
        Mesh.Cell3DDoublePropertyValues = {};
      }
  };
}

#endif // __MESHMATRICES_2D_26CELLS_MOCK_H
