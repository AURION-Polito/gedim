#ifndef __CONFORMERMESHSEGMENT_H
#define __CONFORMERMESHSEGMENT_H

#include "Eigen/Eigen"
#include "IMeshDAO.hpp"
#include "GeometryUtilities.hpp"
#include "IntersectorMesh2DSegment.hpp"
#include "UnionMeshSegment.hpp"
#include "MeshMatrices.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
  class ConformerMeshSegment final
  {
    public:
      struct ConformMesh final
      {
          struct ConformMeshPoint final
          {
              enum Types
              {
                Unknown = 0,
                Original = 1, ///< point belong to original intersection mesh
                Inherited = 2, ///< point belong to other intersection mesh
                External = 3 ///< other points not inherith from intersection mesh
              };

              std::set<unsigned int> Cell2DIds = {};
              std::set<unsigned int> Edge2DIds = {};
              std::set<unsigned int> Vertex2DIds = {};
              Types Type;
          };

          struct ConformMeshSegment final
          {
              std::vector<double> Points = {};
              std::set<unsigned int> Cell2DIds = {};
              std::set<unsigned int> Edge2DIds = {};
          };

          std::map<double, ConformMeshPoint> Points;
          std::vector<ConformMeshSegment> Segments;
      };

      /// \brief convert IntersectionMesh to Curvilinear Coordinates vector
      static void ToCurvilinearCoordinates(const ConformMesh& conformMesh,
                                           std::vector<double>& curvilinearCoordinates);

      static std::string ToString(const ConformMesh& conformMesh);

      static void CreateConformSegments(ConformMesh& result);

      static void Serialize(std::ostream &os,
                            const ConformMesh& mesh);

      static void Deserialize(std::istream &is,
                              ConformMesh& mesh);

    private:
      const Gedim::GeometryUtilities& _geometryUtilities;

      ConformMesh::ConformMeshPoint& InsertNewIntersection(const double& curvilinearCoordinate,
                                                           ConformMesh& result,
                                                           bool& found) const;
      void CreateConformPoints(const Gedim::IntersectorMesh2DSegment::IntersectionMesh& meshIntersection,
                               const Gedim::UnionMeshSegment::UnionMesh& meshUnion,
                               const unsigned int& meshIntersectionPosition,
                               ConformMesh& result);

    public:
      ConformerMeshSegment(const Gedim::GeometryUtilities& geometryUtilities);
      ~ConformerMeshSegment();

      /// \brief Create ConformMesh on segment starting from an intersection mesh and the corresponding union mesh
      /// \param meshIntersection the segment intersection mesh
      /// \param meshUnion the segment mesh union of meshIntersection and an other meshIntersection
      /// \param meshIntersectionPosition the position of meshIntersection inside meshUnion, starting from 0
      /// \param result the resulting conform mesh
      void CreateConformMesh(const Gedim::IntersectorMesh2DSegment::IntersectionMesh& meshIntersection,
                             const Gedim::UnionMeshSegment::UnionMesh& meshUnion,
                             const unsigned int& meshIntersectionPosition,
                             ConformMesh& result);

      /// \brief Insert an external point on conform mesh
      /// \param mesh2D the bidimensional mesh where to add the new point
      /// \param curvilinearCoordinate the curvilinear coordinate of the new point
      /// \param result the resulting conform mesh
      void InsertExternalPoint(const Eigen::Vector3d& segmentOrigin,
                               const Eigen::Vector3d& segmentEnd,
                               const Gedim::IMeshDAO& mesh2D,
                               const double& curvilinearCoordinate,
                               ConformMesh& result);

      /// \brief Update the conformed 1D mesh with updated mesh 2D data
      /// \param mesh2D the updated mesh data
      /// \param conformedMesh the resulting conformed mesh
      void UpdateWithUpdatedMesh2D(const Gedim::IMeshDAO& mesh2D,
                                   ConformMesh& conformedMesh) const;

      /// \brief Update the conformed 1D mesh with active mesh 2D data
      /// \param activeMesh2DData the active mesh data
      /// \param conformedMesh the resulting conformed mesh
      void UpdateWithActiveMesh2D(const Gedim::MeshUtilities::ExtractActiveMeshData& activeMesh2DData,
                                  ConformMesh& conformedMesh) const;

      /// \brief Add Missing Mesh2D Cell0Ds generated from the intersections with other interfaces
      /// \param mesh2D the mesh data
      /// \param conformedMesh the resulting conformed mesh
      /// \note the mesh2D shall be made by only active elements
      /// \return the list of coordinates added
      std::vector<double> AddMissingMesh2DCell0Ds(const Eigen::Vector3d& segmentOrigin,
                                                  const Eigen::Vector3d& segmentTangent,
                                                  const double& segmentSquaredLength,
                                                  const Gedim::IMeshDAO& mesh2D,
                                                  ConformMesh& conformedMesh) const;
  };
}

#endif
