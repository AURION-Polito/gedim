#ifndef __UNIONMESHSEGMENT_H
#define __UNIONMESHSEGMENT_H

#include "GeometryUtilities.hpp"

namespace Gedim
{
  class UnionMeshSegment final
  {
    public:
      struct UnionMesh final {
          struct UnionMeshPoint final
          {
              enum struct Types
              {
                Unknown = 0,
                First = 1,
                Second = 2,
                Both = 3
              };

              Types Type = Types::Unknown;
              std::vector<unsigned int> MeshIndices = {}; ///< vector of size 2 containing in each i the indices in mesh_i
          };

          struct UnionMeshSegment final
          {
              enum struct Types
              {
                Unknown = 0,
                First = 1,
                Second = 2,
                Both = 3
              };

              Types Type = Types::Unknown;
              std::vector<double> Points = {};
              std::vector<unsigned int> MeshIndices = {}; ///< vector of size 2 containing in each i the indices in mesh_i
          };

          std::map<double, UnionMeshPoint> Points = {};
          std::vector<UnionMeshSegment> Segments = {};
      };

      /// \brief convert UnionMesh to Curvilinear Coordinates vector
      static void ToCurvilinearCoordinates(const UnionMesh& unionMesh,
                                           std::vector<double>& curvilinearCoordinates);

      static void ToString(const UnionMesh& unionMesh);

    private:
      const Gedim::GeometryUtilities& _geometryUtilities;

      UnionMesh::UnionMeshPoint& InsertNewIntersection(const double& curvilinearCoordinate,
                                                       UnionMesh& result,
                                                       bool& found);

      void CreateUnionPoints(const std::vector<double>& curvilinearCoordinatesMeshOne,
                             const std::vector<double>& curvilinearCoordinatesMeshTwo,
                             UnionMesh& result);
      void CreateUnionSegments(const std::vector<double>& curvilinearCoordinatesMeshOne,
                               const std::vector<double>& curvilinearCoordinatesMeshTwo,
                               UnionMesh& result);

    public:
      UnionMeshSegment(const Gedim::GeometryUtilities& geometryUtilities);
      ~UnionMeshSegment();

      void CreateUnionMesh(const std::vector<double>& curvilinearCoordinatesMeshOne,
                           const std::vector<double>& curvilinearCoordinatesMeshTwo,
                           UnionMesh& result);
  };
}

#endif
