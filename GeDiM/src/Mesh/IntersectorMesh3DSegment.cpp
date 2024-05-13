#include "IntersectorMesh3DSegment.hpp"
#include "IOUtilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  IntersectorMesh3DSegment::IntersectorMesh3DSegment(const Gedim::GeometryUtilities& geometryUtilities,
                                                     const Gedim::MeshUtilities& meshUtilities) :
    geometryUtilities(geometryUtilities),
    meshUtilities(meshUtilities)
  {
  }
  IntersectorMesh3DSegment::~IntersectorMesh3DSegment()
  {
  }
  // ***************************************************************************
  IntersectorMesh3DSegment::IntersectionMesh::IntersectionMeshPoint& IntersectorMesh3DSegment::InsertNewIntersection(const double& curvilinearCoordinate,
                                                                                                                     std::map<double, IntersectionMesh::IntersectionMeshPoint>& points,
                                                                                                                     bool& found) const
  {
    double foundCoordinate = -1.0;
    for (auto& it : points)
    {
      if (!geometryUtilities.IsValuePositive(abs(it.first - curvilinearCoordinate),
                                             geometryUtilities.Tolerance1D()))
      {
        foundCoordinate = it.first;
        break;
      }
    }

    if (foundCoordinate != -1.0)
    {
      found = true;
      return points[foundCoordinate];
    }

    points.insert(std::make_pair(curvilinearCoordinate,
                                 IntersectionMesh::IntersectionMeshPoint()));
    found = false;
    return points[curvilinearCoordinate];
  }
  // ***************************************************************************
  std::vector<IntersectorMesh3DSegment::IntersectionMesh::IntersectionMeshSegment> IntersectorMesh3DSegment::CreateIntersectionSegments(const std::vector<IntersectionMesh::IntersectionMeshPoint>& mesh1D_points) const
  {
    std::vector<IntersectorMesh3DSegment::IntersectionMesh::IntersectionMeshSegment> mesh1D_segments;

    mesh1D_segments.resize(mesh1D_points.size() - 1);

    std::vector<std::set<unsigned int>> points_cell3Ds(mesh1D_points.size());
    for (unsigned int s = 0; s < mesh1D_points.size(); s++)
    {
      const std::vector<unsigned int>& cell3Ds = mesh1D_points[s].Cell3DIds;
      points_cell3Ds[s] = std::set<unsigned int>(cell3Ds.begin(),
                                                 cell3Ds.end());
    }

    for (unsigned int s = 0; s < mesh1D_segments.size(); s++)
    {
      const std::vector<unsigned int>& origin_cell3Ds = mesh1D_points[s].Cell3DIds;
      const std::vector<unsigned int>& end_cell3Ds = mesh1D_points[s + 1].Cell3DIds;

      std::vector<unsigned int> segment_cell3Ds;
      std::set_intersection(origin_cell3Ds.begin(),
                            origin_cell3Ds.end(),
                            end_cell3Ds.begin(),
                            end_cell3Ds.end(),
                            std::back_inserter(segment_cell3Ds));

      mesh1D_segments[s] =
      {
        { s, s + 1 },
        segment_cell3Ds
      };
    }

    return mesh1D_segments;
  }
  // ***************************************************************************
  unsigned int IntersectorMesh3DSegment::FindSegmentVertexCell3D(const Eigen::Vector3d& vertex,
                                                                 const Gedim::IMeshDAO& mesh3D,
                                                                 const MeshUtilities::MeshGeometricData3D& mesh3D_geometricData) const
  {
    unsigned int cell3D_index_found = mesh3D.Cell3DTotalNumber();

    for (unsigned int c3D_index = 0; c3D_index < mesh3D.Cell3DTotalNumber(); c3D_index++)
    {
      if (!mesh3D.Cell3DIsActive(c3D_index))
        continue;

      const auto& cell3D_vertices = mesh3D_geometricData.Cell3DsVertices.at(c3D_index);
      const auto& cell3D_faces = mesh3D_geometricData.Cell3DsFaces.at(c3D_index);
      const auto& cell3D_faces_3D_vertices = mesh3D_geometricData.Cell3DsFaces3DVertices.at(c3D_index);
      const auto& cell3D_faces_2D_vertices = mesh3D_geometricData.Cell3DsFaces2DVertices.at(c3D_index);
      const auto& cell3D_faces_normal = mesh3D_geometricData.Cell3DsFacesNormals.at(c3D_index);
      const auto& cell3D_faces_normal_direction = mesh3D_geometricData.Cell3DsFacesNormalDirections.at(c3D_index);
      const auto& cell3D_faces_translation = mesh3D_geometricData.Cell3DsFacesTranslations.at(c3D_index);
      const auto& cell3D_faces_rotation = mesh3D_geometricData.Cell3DsFacesRotationMatrices.at(c3D_index);

      const Eigen::MatrixXd cell3D_BoundingBox = geometryUtilities.PointsBoundingBox(cell3D_vertices);

      if (!geometryUtilities.IsPointInBoundingBox(vertex,
                                                  cell3D_BoundingBox))
        continue;

      const GeometryUtilities::PointPolyhedronPositionResult vertex_position =
          geometryUtilities.PointPolyhedronPosition(vertex,
                                                    cell3D_faces,
                                                    cell3D_faces_3D_vertices,
                                                    cell3D_faces_2D_vertices,
                                                    cell3D_faces_normal,
                                                    cell3D_faces_normal_direction,
                                                    cell3D_faces_translation,
                                                    cell3D_faces_rotation);

      switch (vertex_position.Type)
      {
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Outside:
          continue;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace:
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex:
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Inside:
          cell3D_index_found = c3D_index;
          break;
        default:
          throw std::runtime_error("Vertex position not supported");
      }
    }

    Gedim::Output::Assert(cell3D_index_found < mesh3D.Cell3DTotalNumber());

    return cell3D_index_found;
  }
  // ***************************************************************************
  IntersectorMesh3DSegment::IntersectionMesh IntersectorMesh3DSegment::CreateIntersectionMesh(const Eigen::Vector3d& segmentOrigin,
                                                                                              const Eigen::Vector3d& segmentEnd,
                                                                                              const Eigen::Vector3d& segmentTangent,
                                                                                              const Eigen::Vector3d& segmentBarycenter,
                                                                                              const double& segmentLength,
                                                                                              const Gedim::IMeshDAO& mesh3D,
                                                                                              const Gedim::MeshUtilities::MeshGeometricData3D& mesh3D_geometricData) const
  {
    const unsigned int segment_origin_cell3D_index = FindSegmentVertexCell3D(segmentOrigin,
                                                                             mesh3D,
                                                                             mesh3D_geometricData);
    const unsigned int segment_end_cell3D_index = FindSegmentVertexCell3D(segmentEnd,
                                                                          mesh3D,
                                                                          mesh3D_geometricData);

    if (segment_origin_cell3D_index ==
        segment_end_cell3D_index)
    {
      IntersectionMesh mesh1D;

      mesh1D.Points =
      {
        {
          0.0,
          { segment_origin_cell3D_index }
        },
        {
          1.0,
          { segment_origin_cell3D_index }
        }
      };
      mesh1D.Segments = CreateIntersectionSegments(mesh1D.Points);

      return mesh1D;
    }

    IntersectionMesh mesh1D;

    return mesh1D;
  }
  // ***************************************************************************
}
