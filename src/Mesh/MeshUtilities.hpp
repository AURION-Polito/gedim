#ifndef __MeshUtilities_H
#define __MeshUtilities_H

#include "IMeshDAO.hpp"
#include "GeometryUtilities.hpp"

using namespace std;

namespace Gedim
{
  /// \brief MeshUtilities
  /// \copyright See top level LICENSE file for details.
  class MeshUtilities final {
    public:
      struct ExtractActiveMeshData final
      {
          map<unsigned int, unsigned int> OldCell0DToNewCell0D; ///< each pair is {old Cell0D index, new Cell0D index}
          map<unsigned int, unsigned int> OldCell1DToNewCell1D; ///< each pair is {old Cell1D index, new Cell1D index}
          map<unsigned int, unsigned int> OldCell2DToNewCell2D; ///< each pair is {old Cell2D index, new Cell2D index}
          map<unsigned int, unsigned int> OldCell3DToNewCell3D; ///< each pair is {old Cell3D index, new Cell3D index}
          map<unsigned int, unsigned int> NewCell0DToOldCell0D; ///< each pair is {new Cell0D index, old Cell0D index}
          map<unsigned int, unsigned int> NewCell1DToOldCell1D; ///< each pair is {new Cell1D index, old Cell1D index}
          map<unsigned int, unsigned int> NewCell2DToOldCell2D; ///< each pair is {new Cell2D index, old Cell2D index}
          map<unsigned int, unsigned int> NewCell3DToOldCell3D; ///< each pair is {new Cell3D index, old Cell3D index}
      };

    public:
      MeshUtilities();
      ~MeshUtilities();

      /// \brief Extract Active Cells from mesh
      /// \note the data are duplied from mesh to extractedMesh
      void ExtractActiveMesh(IMeshDAO& mesh,
                             ExtractActiveMeshData& extractionData) const;

      /// \brief Fill Mesh 1D From segment Coordinates
      /// \param segmentOrigin the segment origin
      /// \param segmentTangent the segment tangent vector
      /// \param coordinates relative coordinates between [0.0, 1.0]
      /// \param mesh the resulting mesh
      void FillMesh1D(const GeometryUtilities& geometryUtilities,
                      const Vector3d& segmentOrigin,
                      const Vector3d& segmentTangent,
                      const vector<double>& coordinates,
                      IMeshDAO& mesh) const;

      /// \brief Fill a Mesh 2D with vertices and polygons
      /// \param cell0Ds the coordinates as Eigen MatrixXd of cell0Ds, size 3xCell0DTotalNumber()
      /// \param cell1Ds the origin and end as Eigen MatrixXd of cell1Ds, size 2xCell1DTotalNumber()
      /// \param cell2Ds the vertices and edges indices of the cell2Ds ordered counterclockwise, size Cell2DTotalNumber()x2xCell2DNumberVertices()
      void FillMesh2D(const Eigen::MatrixXd& cell0Ds,
                      const Eigen::MatrixXi& cell1Ds,
                      const vector<Eigen::MatrixXi>& cell2Ds,
                      IMeshDAO& mesh) const;
  };

}

#endif // __MeshUtilities_H
