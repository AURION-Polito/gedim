#ifndef __IntersectorMesh3DSegment_H
#define __IntersectorMesh3DSegment_H

#include "Eigen/Eigen"
#include "GeometryUtilities.hpp"
#include "IMeshDAO.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
class IntersectorMesh3DSegment final
{
  public:
    struct IntersectionMesh final
    {
        struct IntersectionMeshPoint final
        {
            double CurvilinearCoordinate;
            std::vector<unsigned int> Cell3DIds;
            std::vector<GeometryUtilities::PointPolyhedronPositionResult> Positions;
        };

        struct IntersectionMeshSegment final
        {
            std::array<unsigned int, 2> PointsIndex;
            std::vector<unsigned int> Cell3DIds;
            std::vector<GeometryUtilities::SegmentPolyhedronPositionResult> Positions;
        };

        std::vector<IntersectionMeshPoint> Points;
        std::vector<IntersectionMeshSegment> Segments;
    };

    struct IntersectionPoint final
    {
        std::set<unsigned int> Cell3DIds;
    };

    struct FindSegmentStartingCell3DResult final
    {
        unsigned int StartingCell3DIndex;
        bool SegmentFullInsideCell3D;
        GeometryUtilities::PointPolyhedronPositionResult Position;
    };

  private:
    const Gedim::GeometryUtilities &geometryUtilities;
    const Gedim::MeshUtilities &meshUtilities;

    IntersectionPoint &CreateOrFindNewIntersection(const double &curvilinearCoordinate,
                                                   std::map<double, IntersectionPoint> &points,
                                                   bool &found) const;

    bool CheckSegmentPolyhedronIntersection(const double &segment_intersection_coordinate,
                                            const Eigen::Vector3d &segment_vertex,
                                            const IMeshDAO &mesh3D,
                                            const unsigned int cell3D_index,
                                            const Eigen::MatrixXd &cell3D_BoundingBox,
                                            const std::vector<Eigen::MatrixXi> &cell3D_faces,
                                            const std::vector<Eigen::MatrixXd> &cell3D_faces_3D_vertices,
                                            const std::vector<Eigen::MatrixXd> &cell3D_faces_2D_vertices,
                                            const std::vector<Eigen::Vector3d> &cell3D_faces_normal,
                                            const std::vector<bool> &cell3D_faces_normal_direction,
                                            const std::vector<Eigen::Vector3d> &cell3D_faces_translation,
                                            const std::vector<Eigen::Matrix3d> &cell3D_faces_rotation,
                                            const std::vector<Eigen::MatrixXd> &cell3D_tetrahedrons,
                                            std::map<double, IntersectionPoint> &mesh1D_intersections,
                                            std::list<unsigned int> &cell3Ds_index) const;

    bool CheckSegmentFaceIntersection(const Eigen::Vector3d &segment_origin,
                                      const Eigen::Vector3d &segment_tangent,
                                      const GeometryUtilities::PointSegmentPositionTypes segment_intersection_position,
                                      const double &segment_intersection_coordinate,
                                      const Eigen::Matrix3d &face_rotation,
                                      const Eigen::Vector3d &face_translation,
                                      const Eigen::MatrixXd &face_2D_vertices,
                                      const IMeshDAO &mesh3D,
                                      const unsigned int cell3D_index,
                                      const unsigned int cell2D_index,
                                      std::map<double, IntersectionPoint> &mesh1D_intersections,
                                      std::list<unsigned int> &cell3Ds_index) const;

    bool CheckSegmentEdgeIntersection(const Gedim::GeometryUtilities::PointSegmentPositionTypes segment_intersection_position,
                                      const double &segment_intersection_coordinate,
                                      const Gedim::IMeshDAO &mesh3D,
                                      const unsigned int cell3D_index,
                                      const unsigned int cell2D_index,
                                      std::map<double, IntersectionPoint> &mesh1D_intersections,
                                      std::list<unsigned int> &cell3Ds_index) const;

    void InsertNewIntersection(const double &segment_intersection_coordinate,
                               const IMeshDAO &mesh3D,
                               const unsigned int cell3D_index,
                               const unsigned int cell2D_index,
                               std::map<double, IntersectionPoint> &mesh1D_intersections,
                               std::list<unsigned int> &cell3Ds_index) const;

    std::vector<IntersectionMesh::IntersectionMeshSegment> CreateIntersectionSegments(
        const Gedim::IMeshDAO &mesh3D,
        const std::vector<IntersectionMesh::IntersectionMeshPoint> &mesh1D_points) const;

    std::vector<IntersectionMesh::IntersectionMeshPoint> FindSegmentIntersectionPoints(
        const Eigen::Vector3d &segmentOrigin,
        const Eigen::Vector3d &segmentEnd,
        const Eigen::Vector3d &segmentTangent,
        const Gedim::IMeshDAO &mesh3D,
        const Gedim::MeshUtilities::MeshGeometricData3D &mesh3D_geometricData,
        const unsigned int starting_cell3D_index) const;

    void SegmentCell3DIntersection(const Eigen::Vector3d &segmentOrigin,
                                   const Eigen::Vector3d &segmentEnd,
                                   const Eigen::Vector3d &segmentTangent,
                                   const Gedim::IMeshDAO &mesh3D,
                                   const Gedim::MeshUtilities::MeshGeometricData3D &mesh3D_geometricData,
                                   const unsigned int cell3D_index,
                                   std::map<double, IntersectionPoint> &mesh1D_intersections,
                                   std::list<unsigned int> &cell3Ds_index) const;

  public:
    IntersectorMesh3DSegment(const Gedim::GeometryUtilities &geometryUtilities, const Gedim::MeshUtilities &meshUtilities);
    ~IntersectorMesh3DSegment();

    static std::vector<double> ToCurvilinearCoordinates(const IntersectorMesh3DSegment::IntersectionMesh &intersectingMesh);
    static std::vector<std::vector<unsigned int>> MeshSegmentsCell3Ds(const IntersectorMesh3DSegment::IntersectionMesh &intersectingMesh);
    static std::string ToString(const IntersectorMesh3DSegment::IntersectionMesh &intersectingMesh);

    IntersectionMesh CreateIntersectionMesh(const Eigen::Vector3d &segmentOrigin,
                                            const Eigen::Vector3d &segmentEnd,
                                            const Eigen::Vector3d &segmentTangent,
                                            const Gedim::IMeshDAO &mesh3D,
                                            const Gedim::MeshUtilities::MeshGeometricData3D &mesh3D_geometricData) const;

    FindSegmentStartingCell3DResult FindSegmentStartingCell3D(const Eigen::Vector3d &segmentOrigin,
                                                              const Eigen::Vector3d &segmentEnd,
                                                              const Gedim::IMeshDAO &mesh3D,
                                                              const Gedim::MeshUtilities::MeshGeometricData3D &mesh3D_geometricData) const;
};
} // namespace Gedim

#endif
