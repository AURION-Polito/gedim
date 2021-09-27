#ifndef __UNIONMESHSEGMENT_H
#define __UNIONMESHSEGMENT_H

#include "GeometryUtilities.hpp"

using namespace std;

namespace Gedim
{
  class UnionMeshSegment final
  {
    public:
      struct UnionMesh final {
          struct UnionMeshPoint final {
              enum Types {
                Unknown = 0,
                First = 1,
                Second = 2,
                Both = 3
              };

              Types Type;
              vector<unsigned int> MeshIndices; ///< vector of size 2 containing in each i the indices in mesh_i
          };

          map<double, UnionMeshPoint> Points;
      };

      /// \brief convert UnionMesh to Curvilinear Coordinates vector
      static void ToCurvilinearCoordinates(const UnionMesh& unionMesh,
                                           vector<double>& curvilinearCoordinates);

      static void ToString(const UnionMesh& unionMesh);

    private:
      const Gedim::GeometryUtilities& _geometryUtilities;

      UnionMesh::UnionMeshPoint& InsertNewIntersection(const double& curvilinearCoordinate,
                                                       UnionMesh& result,
                                                       bool& found);

    public:
      UnionMeshSegment(const Gedim::GeometryUtilities& geometryUtilities);
      ~UnionMeshSegment();

      void CreateUnionMesh(const vector<double>& curvilinearCoordinatesMeshOne,
                           const vector<double>& curvilinearCoordinatesMeshTwo,
                           UnionMesh& result);
  };
}

#endif
