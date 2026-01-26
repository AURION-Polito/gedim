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

#ifndef __IntersectorMesh2DSegment_H
#define __IntersectorMesh2DSegment_H

#include "Eigen/Eigen"
#include "GeometryUtilities.hpp"
#include "IMeshDAO.hpp"

namespace Gedim
{
class IntersectorMesh2DSegment final
{
  public:
    struct IntersectionMesh final
    {
        struct IntersectionMeshPoint final
        {
            std::set<unsigned int> Cell2DIds = {};
            std::set<unsigned int> Cell1DIds = {};
            std::set<unsigned int> Cell0DIds = {};
        };

        struct IntersectionMeshSegment final
        {
            std::vector<double> Points = {};
            std::set<unsigned int> Cell2DIds = {};
            std::set<unsigned int> Cell1DIds = {};
        };

        std::map<double, IntersectionMeshPoint> Points;
        std::vector<IntersectionMeshSegment> Segments;
    };

    /// \brief convert IntersectionMesh to Curvilinear Coordinates vector
    static void ToCurvilinearCoordinates(const IntersectionMesh &intersectingMesh, std::vector<double> &curvilinearCoordinates);

    static void ToString(const IntersectionMesh &intersectingMesh);

  private:
    const Gedim::IMeshDAO &_mesh;
    const Gedim::GeometryUtilities &_geometryUtilities;

    IntersectionMesh::IntersectionMeshPoint &InsertNewIntersection(const double &curvilinearCoordinate,
                                                                   IntersectionMesh &result,
                                                                   bool &found);

    void CheckOriginAndEndSegmentPosition(const Eigen::Vector3d &segmentOrigin,
                                          const Eigen::Vector3d &segmentEnd,
                                          IntersectionMesh &result,
                                          const bool concave = false);

    void CreateIntersectionPoints(const Eigen::Vector3d &segmentOrigin,
                                  const Eigen::Vector3d &segmentEnd,
                                  const Eigen::Vector3d &segmentTangent,
                                  const Eigen::Vector3d &segmentBarycenter,
                                  const double &segmentLength,
                                  IntersectionMesh &result);
    void CreateIntersectionSegments(IntersectionMesh &result);

    void CheckOriginAndEndPointExistence(IntersectionMesh &result);

    void SmoothIntersections(const Eigen::Vector3d &segmentOrigin, const Eigen::Vector3d &segmentEnd, IntersectionMesh &result);

  public:
    IntersectorMesh2DSegment(const Gedim::IMeshDAO &mesh, const Gedim::GeometryUtilities &geometryUtilities);
    ~IntersectorMesh2DSegment();

    void CreateIntersectionMesh(const Eigen::Vector3d &segmentOrigin,
                                const Eigen::Vector3d &segmentEnd,
                                const Eigen::Vector3d &segmentTangent,
                                const Eigen::Vector3d &segmentBarycenter,
                                const double &segmentLength,
                                IntersectionMesh &result,
                                const bool concave = false);
};
} // namespace Gedim

#endif
