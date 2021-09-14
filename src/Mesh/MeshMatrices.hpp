#ifndef __MeshMatrices_H
#define __MeshMatrices_H

#include "Output.hpp"
#include "Eigen"

using namespace std;

namespace Gedim
{
  struct MeshMatrices final
  {
      unsigned int Dimension = 0; ///< Geometric dimension
      unsigned int NumberCell0D = 0; ///< number of Cell0D
      vector<double> Cell0DCoordinates = {}; ///< Cell0D coordinates, size 3 x NumberCell0D (x,y,z)
      vector<unsigned int> Cell0DMarkers = {}; ///< Cell0D markers, size 1 x NumberCell0D (marker)
      vector<bool> ActiveCell0D = {}; ///< active Cell0D
      map<unsigned int, set<unsigned int>> UpdatedCell0Ds = {}; ///< for each cell0D the list to the new cell0Ds
      vector<string> Cell0DDoublePropertyIds = {}; ///< Cell0D double property id - double property index
      map<string, unsigned int> Cell0DDoublePropertyIndices = {}; ///< Cell0D double property id - double property index
      vector<vector<unsigned int>> Cell0DDoublePropertySizes = {}; ///< Cell0D double property sizes
      vector<vector<double>> Cell0DDoublePropertyValues = {}; ///< Cell0D double property values
      unsigned int NumberCell1D = 0; ///< number of Cell1D
      vector<unsigned int> Cell1DVertices = {}; ///< Cell1D vertices indices, size 2 x NumberCell1D (fromId,toId)
      vector<unsigned int> NumberCell1DNeighbourCell2D = {}; ///< Cell1D neighbour Cell2D indices per cell, size  1 x NumberCell1D + 1
      vector<unsigned int> Cell1DNeighbourCell2Ds = {}; ///< Cell1D neighbour Cell2D indices, size 1 x NumberCell1DNeighbourCell2D[NumberCell1D]
      vector<unsigned int> Cell1DMarkers = {}; ///< Cell1D propertoes, size 1 x NumberCell1D (marker)
      vector<bool> ActiveCell1D = {}; ///< active Cell1D
      map<unsigned int, set<unsigned int>> UpdatedCell1Ds = {}; ///< for each cell1D the list to the new cell1Ds
      vector<string> Cell1DDoublePropertyIds = {}; ///< Cell1D double property id - double property index
      map<string, unsigned int> Cell1DDoublePropertyIndices = {}; ///< Cell1D double property id - double property index
      vector<vector<unsigned int>> Cell1DDoublePropertySizes = {}; ///< Cell1D double property sizes
      vector<vector<double>> Cell1DDoublePropertyValues = {}; ///< Cell1D double property values
      unsigned int NumberCell2D = 0; ///< number of Cell2D
      vector<unsigned int> NumberCell2DVertices = {}; ///< number of Vertices per Cell2D, size 1 x NumberCell2D + 1
      vector<unsigned int> NumberCell2DEdges = {}; ///< number of Edges per Cell2D, size 1 x NumberCell2D + 1
      vector<unsigned int> Cell2DVertices = {}; ///< Cell2D Vertices indices, size 1 x NumberCell2DVertices[NumberCell2D]
      vector<unsigned int> Cell2DEdges = {}; ///< Cell2D Cell1D indices, size 1 x NumberCell2DEdges[NumberCell2D]
      vector<unsigned int> Cell2DMarkers = {}; ///< Cell2D markers, size 1 x NumberCell2D (marker)
      vector<bool> ActiveCell2D = {}; ///< active Cell2D
      map<unsigned int, set<unsigned int>> UpdatedCell2Ds = {}; ///< for each cell2D the list to the new cell2Ds
      vector<string> Cell2DDoublePropertyIds = {}; ///< Cell2D double property id - double property index
      map<string, unsigned int> Cell2DDoublePropertyIndices = {}; ///< Cell2D double property id - double property index
      vector<vector<unsigned int>> Cell2DDoublePropertySizes = {}; ///< Cell2D double property sizes
      vector<vector<double>> Cell2DDoublePropertyValues = {}; ///< Cell2D double property values
      unsigned int NumberCell3D = 0; ///< number of Cell3D
      vector<unsigned int> NumberCell3DVertices = {}; ///< number of Vertices per Cell3D, size 1 x NumberCell3D + 1
      vector<unsigned int> NumberCell3DEdges = {}; ///< number of Edges per Cell3D, size 1 x NumberCell3D + 1
      vector<unsigned int> NumberCell3DFaces = {}; ///< number of Faces per Cell3D, size 1 x NumberCell3D + 1
      vector<unsigned int> Cell3DVertices = {}; ///< Cell3D Cell0D indices, size 1 x NumberCell3DVertices[NumberCell3D]
      vector<unsigned int> Cell3DEdges = {}; ///< Cell3D Cell1D indices, size 1 x NumberCell3DEdges[NumberCell3D]
      vector<unsigned int> Cell3DFaces = {}; ///< Cell3D Cell2D indices, size 1 x NumberCell3DFaces[NumberCell3D]
      vector<unsigned int> Cell3DMarkers = {}; ///< Cell3D markers, size 1 x NumberCell3D (marker)
      vector<bool> ActiveCell3D = {}; ///< active Cell3D
      map<unsigned int, set<unsigned int>> UpdatedCell3Ds = {}; ///< for each cell3D the list to the new cell3Ds
      vector<string> Cell3DDoublePropertyIds = {}; ///< Cell3D double property id - double property index
      map<string, unsigned int> Cell3DDoublePropertyIndices = {}; ///< Cell3D double property id - double property index
      vector<vector<unsigned int>> Cell3DDoublePropertySizes = {}; ///< Cell3D double property sizes
      vector<vector<double>> Cell3DDoublePropertyValues = {}; ///< Cell3D double property values
  };
}

#endif
