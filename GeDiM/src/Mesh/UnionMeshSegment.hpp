// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __UNIONMESHSEGMENT_H
#define __UNIONMESHSEGMENT_H

#include "GeometryUtilities.hpp"

namespace Gedim
{
class UnionMeshSegment final
{
  public:
    struct UnionMesh final
    {
        struct UnionMeshPoint final
        {
            enum struct Types
            {
                Unknown = 0,
                First = 1,
                Second = 2,
                Both = 3
            };

            Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::Types Type = Types::Unknown;
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

            Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types Type = Types::Unknown;
            std::vector<double> Points = {};
            std::vector<unsigned int> MeshIndices = {}; ///< vector of size 2 containing in each i the indices in mesh_i
        };

        std::map<double, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint> Points = {};
        std::vector<Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment> Segments = {};
    };

    /// \brief convert UnionMesh to Curvilinear Coordinates vector
    static void ToCurvilinearCoordinates(const Gedim::UnionMeshSegment::UnionMesh &unionMesh, std::vector<double> &curvilinearCoordinates);

    static void ToString(const Gedim::UnionMeshSegment::UnionMesh &unionMesh);

  private:
    const Gedim::GeometryUtilities &_geometryUtilities;

    UnionMesh::UnionMeshPoint &InsertNewIntersection(const double &curvilinearCoordinate, UnionMesh &result, bool &found);

    void CreateUnionPoints(const std::vector<double> &curvilinearCoordinatesMeshOne,
                           const std::vector<double> &curvilinearCoordinatesMeshTwo,
                           Gedim::UnionMeshSegment::UnionMesh &result);
    void CreateUnionSegments(const std::vector<double> &curvilinearCoordinatesMeshOne,
                             const std::vector<double> &curvilinearCoordinatesMeshTwo,
                             Gedim::UnionMeshSegment::UnionMesh &result);

  public:
    UnionMeshSegment(const Gedim::GeometryUtilities &geometryUtilities);
    ~UnionMeshSegment();

    void CreateUnionMesh(const std::vector<double> &curvilinearCoordinatesMeshOne,
                         const std::vector<double> &curvilinearCoordinatesMeshTwo,
                         Gedim::UnionMeshSegment::UnionMesh &result);
};
} // namespace Gedim

#endif
