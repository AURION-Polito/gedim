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

      struct MeshGeometricData final
      {
          vector<Eigen::MatrixXd> Cell2DsVertices; ///< cell2D vertices coordinates
          vector<vector<Eigen::Matrix3d>> Cell2DsTriangulations; ///< cell2D triangulations
          vector<double> Cell2DsAreas; ///< cell2D areas
          vector<Eigen::Vector3d> Cell2DsCentroids; ///< cell2D centroids
          vector<double> Cell2DsDiameters; ///< cell2D diameters
          vector<vector<bool>> Cell2DsEdgeDirections; ///< cell2D edge directions
          vector<Eigen::VectorXd> Cell2DsEdgeLengths; ///< cell2D edge lenghts
          vector<Eigen::MatrixXd> Cell2DsEdgeTangents; ///< cell2D edge tangents
          vector<Eigen::MatrixXd> Cell2DsEdgeNormals; ///< cell2D edge normals
      };

    public:
      MeshUtilities();
      ~MeshUtilities();

      /// \brief Extract Active Cells from mesh
      /// \note the resulting mesh has no inactive elements
      void ExtractActiveMesh(IMeshDAO& mesh,
                             ExtractActiveMeshData& extractionData) const;

      /// \brief Fill Mesh 1D From segment Coordinates
      /// \param segmentOrigin the segment origin
      /// \param segmentTangent the segment tangent vector
      /// \param coordinates relative coordinates between [0.0, 1.0]
      /// \param mesh the resulting mesh
      void FillMesh1D(const GeometryUtilities& geometryUtilities,
                      const Eigen::Vector3d& segmentOrigin,
                      const Eigen::Vector3d& segmentTangent,
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

      /// \brief Create a Mesh 2D with a polygon
      /// \param polygonVertices the polygon coordinates, size 3xNumPolygonVertices()
      /// \param vertexMarkers mesh markers of vertices, size 1xNumPolygonVertices()
      /// \param edgeMarkers mesh markers of edges, size 1xNumPolygonVertices()
      void Mesh2DFromPolygon(const Eigen::MatrixXd& polygonVertices,
                             const vector<unsigned int> vertexMarkers,
                             const vector<unsigned int> edgeMarkers,
                             IMeshDAO& mesh) const;

      /// \brief Extract the mesh Cell2D Roots
      /// \param mesh the mesh
      /// \return the root cell for each cell2D, size 1xCell2DTotalNumber()
      vector<unsigned int> MeshCell2DRoots(const IMeshDAO& mesh) const;

      /// \brief Fill Mesh2D Geometric Data given a mesh with convex mesh cells
      /// \param convexMesh the convex mesh
      /// \return the MeshGeometricData computed
      MeshGeometricData FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                const IMeshDAO& convexMesh) const;

      /// \brief Fill Mesh2D Geometric Data starting given a mesh with non convex mesh cells and its convex sub-mesh cells
      /// \param mesh the mesh
      /// \param convexMesh the convex mesh cells of mesh
      /// \param meshCell2DToConvexCell2DIndices the collection of convex cell2Ds for each mesh cell2D
      /// \return the MeshGeometricData computed
      MeshGeometricData FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                const IMeshDAO& mesh,
                                                const IMeshDAO& convexMesh,
                                                const vector<vector<unsigned int>>& meshCell2DToConvexCell2DIndices) const;
  };

}

#endif // __MeshUtilities_H
