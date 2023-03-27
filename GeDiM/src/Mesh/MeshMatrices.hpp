#ifndef __MeshMatrices_H
#define __MeshMatrices_H

#include "IOUtilities.hpp"
#include "Eigen/Eigen"

namespace Gedim
{
  struct MeshMatrices final
  {
      unsigned int Dimension = 0; ///< Geometric dimension
      unsigned int NumberCell0D = 0; ///< number of Cell0D
      std::vector<double> Cell0DCoordinates = {}; ///< Cell0D coordinates, size 3 x NumberCell0D (x,y,z)
      std::vector<unsigned int> Cell0DMarkers = {}; ///< Cell0D markers, size 1 x NumberCell0D (marker)
      std::vector<unsigned int> NumberCell0DNeighbourCell1D = {}; ///< Cell0D neighbour Cell1D indices per cell, size  1 x NumberCell0D + 1
      std::vector<unsigned int> Cell0DNeighbourCell1Ds = {}; ///< Cell0D neighbour Cell1D indices, size 1 x NumberCell0DNeighbourCell1D[NumberCell0D]
      std::vector<unsigned int> NumberCell0DNeighbourCell2D = {}; ///< Cell0D neighbour Cell2D indices per cell, size  1 x NumberCell0D + 1
      std::vector<unsigned int> Cell0DNeighbourCell2Ds = {}; ///< Cell0D neighbour Cell2D indices, size 1 x NumberCell0DNeighbourCell2D[NumberCell0D]
      std::vector<bool> ActiveCell0D = {}; ///< active Cell0D
      std::map<unsigned int, std::set<unsigned int>> UpdatedCell0Ds = {}; ///< for each cell0D the list to the new cell0Ds
      std::vector<std::string> Cell0DDoublePropertyIds = {}; ///< Cell0D double property id - double property index
      std::map<std::string, unsigned int> Cell0DDoublePropertyIndices = {}; ///< Cell0D double property id - double property index
      std::vector<std::vector<unsigned int>> Cell0DDoublePropertySizes = {}; ///< Cell0D double property sizes
      std::vector<std::vector<double>> Cell0DDoublePropertyValues = {}; ///< Cell0D double property values
      unsigned int NumberCell1D = 0; ///< number of Cell1D
      std::vector<unsigned int> Cell1DVertices = {}; ///< Cell1D vertices indices, size 2 x NumberCell1D (fromId,toId)
      Eigen::SparseMatrix<unsigned int> Cell1DAdjacency; ///< Cell1D adjacency, each element is the Cell1D index, size NumberCell0D x NumberCell0D
      std::vector<unsigned int> NumberCell1DNeighbourCell2D = {}; ///< Cell1D neighbour Cell2D indices per cell, size  1 x NumberCell1D + 1
      std::vector<unsigned int> Cell1DNeighbourCell2Ds = {}; ///< Cell1D neighbour Cell2D indices, size 1 x NumberCell1DNeighbourCell2D[NumberCell1D]
      std::vector<unsigned int> Cell1DMarkers = {}; ///< Cell1D propertoes, size 1 x NumberCell1D (marker)
      std::vector<bool> ActiveCell1D = {}; ///< active Cell1D
      std::vector<unsigned int> Cell1DOriginalCell1Ds = {}; ///< for each cell1D the index of original cell1D, NumberCell1D is the default value (no original cell), size 1 x NumberCell1D
      std::map<unsigned int, std::set<unsigned int>> UpdatedCell1Ds = {}; ///< for each cell1D the list to the new cell1Ds
      std::vector<std::string> Cell1DDoublePropertyIds = {}; ///< Cell1D double property id - double property index
      std::map<std::string, unsigned int> Cell1DDoublePropertyIndices = {}; ///< Cell1D double property id - double property index
      std::vector<std::vector<unsigned int>> Cell1DDoublePropertySizes = {}; ///< Cell1D double property sizes
      std::vector<std::vector<double>> Cell1DDoublePropertyValues = {}; ///< Cell1D double property values
      unsigned int NumberCell2D = 0; ///< number of Cell2D
      std::vector<unsigned int> NumberCell2DVertices = {}; ///< number of Vertices per Cell2D, size 1 x NumberCell2D + 1
      std::vector<unsigned int> NumberCell2DEdges = {}; ///< number of Edges per Cell2D, size 1 x NumberCell2D + 1
      std::vector<unsigned int> Cell2DVertices = {}; ///< Cell2D Vertices indices, size 1 x NumberCell2DVertices[NumberCell2D]
      std::vector<unsigned int> Cell2DEdges = {}; ///< Cell2D Cell1D indices, size 1 x NumberCell2DEdges[NumberCell2D]
      std::vector<unsigned int> Cell2DMarkers = {}; ///< Cell2D markers, size 1 x NumberCell2D (marker)
      std::vector<bool> ActiveCell2D = {}; ///< active Cell2D
      std::vector<unsigned int> NumberCell2DSubdivision = {}; ///< number of sub-division per Cell2D, size 1 x NumberCell2D + 1
      std::vector<unsigned int> Cell2DSubdivision = {}; ///< Sub-division of Cell2Ds, used for Concave polygons, size 1 x NumberCell2DSubdivision[NumberCell2D]
      std::vector<unsigned int> Cell2DOriginalCell2Ds = {}; ///< for each cell2D the index of original cell2D, NumberCell2D is the default value (no original cell), size 1 x NumberCell2D
      std::map<unsigned int, std::set<unsigned int>> UpdatedCell2Ds = {}; ///< for each cell2D the list to the new cell2Ds
      std::vector<std::string> Cell2DDoublePropertyIds = {}; ///< Cell2D double property id - double property index
      std::map<std::string, unsigned int> Cell2DDoublePropertyIndices = {}; ///< Cell2D double property id - double property index
      std::vector<std::vector<unsigned int>> Cell2DDoublePropertySizes = {}; ///< Cell2D double property sizes
      std::vector<std::vector<double>> Cell2DDoublePropertyValues = {}; ///< Cell2D double property values
      unsigned int NumberCell3D = 0; ///< number of Cell3D
      std::vector<unsigned int> NumberCell3DVertices = {}; ///< number of Vertices per Cell3D, size 1 x NumberCell3D + 1
      std::vector<unsigned int> NumberCell3DEdges = {}; ///< number of Edges per Cell3D, size 1 x NumberCell3D + 1
      std::vector<unsigned int> NumberCell3DFaces = {}; ///< number of Faces per Cell3D, size 1 x NumberCell3D + 1
      std::vector<unsigned int> Cell3DVertices = {}; ///< Cell3D Cell0D indices, size 1 x NumberCell3DVertices[NumberCell3D]
      std::vector<unsigned int> Cell3DEdges = {}; ///< Cell3D Cell1D indices, size 1 x NumberCell3DEdges[NumberCell3D]
      std::vector<unsigned int> Cell3DFaces = {}; ///< Cell3D Cell2D indices, size 1 x NumberCell3DFaces[NumberCell3D]
      std::vector<unsigned int> Cell3DMarkers = {}; ///< Cell3D markers, size 1 x NumberCell3D (marker)
      std::vector<bool> ActiveCell3D = {}; ///< active Cell3D
      std::map<unsigned int, std::set<unsigned int>> UpdatedCell3Ds = {}; ///< for each cell3D the list to the new cell3Ds
      std::vector<std::string> Cell3DDoublePropertyIds = {}; ///< Cell3D double property id - double property index
      std::map<std::string, unsigned int> Cell3DDoublePropertyIndices = {}; ///< Cell3D double property id - double property index
      std::vector<std::vector<unsigned int>> Cell3DDoublePropertySizes = {}; ///< Cell3D double property sizes
      std::vector<std::vector<double>> Cell3DDoublePropertyValues = {}; ///< Cell3D double property values
  };
}

#endif
