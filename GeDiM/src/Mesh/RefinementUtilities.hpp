#ifndef __RefinementUtilities_H
#define __RefinementUtilities_H

#include "IMeshDAO.hpp"
#include "GeometryUtilities.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
  /// \brief RefinementUtilities
  /// \copyright See top level LICENSE file for details.
  class RefinementUtilities final
  {
    private:
      const GeometryUtilities& geometryUtilities;
      const MeshUtilities& meshUtilities;

      /// \brief Split Triangle From Four Vertices
      /// \param cell2DIndex the cell2D to split
      /// \param splitCell1DsIndex the new cell1Ds index
      /// \param cell1DDirection the original cell1D direction
      /// \param oppositeVertexIndex the vertex opposite to cell1D
      /// \param cell0DOppositeIndex the mesh cell0D index opposite to cell1D
      /// \param newCell0DIndex the new cell0D index
      /// \param mesh the mesh
      void SplitTriangle_FromFourVertices(const unsigned int& cell2DIndex,
                                          const std::vector<unsigned int>& splitCell1DsIndex,
                                          const bool& cell1DDirection,
                                          const unsigned int& oppositeVertexIndex,
                                          const unsigned int& cell0DOppositeIndex,
                                          const unsigned int& newCell0DIndex,
                                          IMeshDAO& mesh) const;

    public:
      RefinementUtilities(const GeometryUtilities& geometryUtilities,
                          const MeshUtilities& meshUtilities);
      ~RefinementUtilities();

      /// \brief Refine Triangle Cell2D By Max Edge
      /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
      /// \param cell2DEdgeLengths the cell2D edge lenghts
      /// \param mesh the mesh to be updated
      void RefineTriangleCellByMaxEdge(const unsigned int& cell2DIndex,
                                       const Eigen::VectorXd& cell2DEdgesLength,
                                       const std::vector<bool>& cell2DEdgesDirection,
                                       IMeshDAO& mesh) const;
  };

}

#endif // __RefinementUtilities_H
