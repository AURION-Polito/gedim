#include "IntersectorMesh3DSegment.hpp"
#include "IOUtilities.hpp"

namespace Gedim
{
// ***************************************************************************
IntersectorMesh3DSegment::IntersectorMesh3DSegment(const Gedim::GeometryUtilities &geometryUtilities,
                                                   const Gedim::MeshUtilities &meshUtilities)
    : geometryUtilities(geometryUtilities), meshUtilities(meshUtilities)
{
}
IntersectorMesh3DSegment::~IntersectorMesh3DSegment()
{
}
// ***************************************************************************
std::vector<double> IntersectorMesh3DSegment::ToCurvilinearCoordinates(const IntersectorMesh3DSegment::IntersectionMesh &intersectingMesh)
{
    std::vector<double> curvilinearCoordinates(intersectingMesh.Points.size());

    for (unsigned int i = 0; i < intersectingMesh.Points.size(); i++)
        curvilinearCoordinates[i] = intersectingMesh.Points[i].CurvilinearCoordinate;

    return curvilinearCoordinates;
}
// ***************************************************************************
std::vector<std::vector<unsigned int>> IntersectorMesh3DSegment::MeshSegmentsCell3Ds(const IntersectionMesh &intersectingMesh)
{
    std::vector<std::vector<unsigned int>> segment_Cell1Ds_Cell3Ds(intersectingMesh.Segments.size());

    for (unsigned int s = 0; s < intersectingMesh.Segments.size(); s++)
        segment_Cell1Ds_Cell3Ds[s] = intersectingMesh.Segments[s].Cell3DIds;

    return segment_Cell1Ds_Cell3Ds;
}
// ***************************************************************************
std::string IntersectorMesh3DSegment::ToString(const IntersectorMesh3DSegment::IntersectionMesh &intersectingMesh)
{
    std::ostringstream stream;

    stream.precision(16);
    stream << "Points:\n";
    for (const auto &intersection_point : intersectingMesh.Points)
    {
        stream << std::scientific << "{ ";
        stream << "Coordinate: " << intersection_point.CurvilinearCoordinate << "; ";
        stream << "Value: C3D: " << intersection_point.Cell3DIds << " }\n";
    };

    stream << "Segments:\n";
    for (const auto &intersection_segment : intersectingMesh.Segments)
    {
        stream << std::scientific << "{ ";
        stream << "Points: [" << intersection_segment.PointsIndex.at(0) << ", ";
        stream << intersection_segment.PointsIndex.at(1) << "]; ";
        stream << "C3D: " << intersection_segment.Cell3DIds << " }\n";
    };

    return stream.str();
}
// ***************************************************************************
IntersectorMesh3DSegment::IntersectionPoint &IntersectorMesh3DSegment::CreateOrFindNewIntersection(const double &curvilinearCoordinate,
                                                                                                   std::map<double, IntersectionPoint> &points,
                                                                                                   bool &found) const
{
    double foundCoordinate = -1.0;
    for (auto &it : points)
    {
        if (!geometryUtilities.IsValuePositive(std::abs(it.first - curvilinearCoordinate), geometryUtilities.Tolerance1D()))
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

    points.insert(std::make_pair(curvilinearCoordinate, IntersectionPoint()));
    found = false;
    return points[curvilinearCoordinate];
}
// ***************************************************************************
bool IntersectorMesh3DSegment::CheckSegmentPolyhedronIntersection(const double &segment_intersection_coordinate,
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
                                                                  std::list<unsigned int> &cell3Ds_index) const
{
    if (!geometryUtilities.IsPointInBoundingBox(segment_vertex, cell3D_BoundingBox))
        return false;

    const GeometryUtilities::PointPolyhedronPositionResult segment_vertex_position =
        geometryUtilities.PointPolyhedronPosition(segment_vertex,
                                                  cell3D_faces,
                                                  cell3D_faces_3D_vertices,
                                                  cell3D_faces_2D_vertices,
                                                  cell3D_faces_normal,
                                                  cell3D_faces_normal_direction,
                                                  cell3D_faces_translation,
                                                  cell3D_faces_rotation,
                                                  cell3D_tetrahedrons);

    switch (segment_vertex_position.Type)
    {
    case GeometryUtilities::PointPolyhedronPositionResult::Types::Outside:
    case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace:
    case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
    case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex:
        return false;
    case GeometryUtilities::PointPolyhedronPositionResult::Types::Inside: {
        InsertNewIntersection(segment_intersection_coordinate, mesh3D, cell3D_index, mesh3D.Cell2DTotalNumber(), mesh1D_intersections, cell3Ds_index);
        return true;
    }
    default:
        throw std::runtime_error("Segment origin position not supported");
    }
}
// ***************************************************************************
bool IntersectorMesh3DSegment::CheckSegmentFaceIntersection(const Eigen::Vector3d &segment_origin,
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
                                                            std::list<unsigned int> &cell3Ds_index) const
{
    switch (segment_intersection_position)
    {
    case GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin:
    case GeometryUtilities::PointSegmentPositionTypes::InsideSegment:
    case GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd: {

        const Eigen::Vector3d segment_intersection = segment_origin + segment_intersection_coordinate * segment_tangent;
        const Eigen::Vector3d segment_intersection_face_2D =
            geometryUtilities.RotatePointsFrom3DTo2D(segment_intersection, face_rotation.transpose(), face_translation);

        if (geometryUtilities.IsPointInsidePolygon_RayCasting(segment_intersection_face_2D, face_2D_vertices))
        {
            InsertNewIntersection(segment_intersection_coordinate, mesh3D, cell3D_index, cell2D_index, mesh1D_intersections, cell3Ds_index);
        }
    }
        return true;
    case GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin:
    case GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineAfterEnd:
    case GeometryUtilities::PointSegmentPositionTypes::LeftTheSegment:
    case GeometryUtilities::PointSegmentPositionTypes::RightTheSegment:
        return false;
    default:
        throw std::runtime_error("segment face single intersection not supported");
    }
}
// ***************************************************************************
bool IntersectorMesh3DSegment::CheckSegmentEdgeIntersection(const Gedim::GeometryUtilities::PointSegmentPositionTypes segment_intersection_position,
                                                            const double &segment_intersection_coordinate,
                                                            const Gedim::IMeshDAO &mesh3D,
                                                            const unsigned int cell3D_index,
                                                            const unsigned int cell2D_index,
                                                            std::map<double, IntersectionPoint> &mesh1D_intersections,
                                                            std::list<unsigned int> &cell3Ds_index) const
{
    switch (segment_intersection_position)
    {
    case GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin:
    case GeometryUtilities::PointSegmentPositionTypes::InsideSegment:
    case GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd: {
        InsertNewIntersection(segment_intersection_coordinate, mesh3D, cell3D_index, cell2D_index, mesh1D_intersections, cell3Ds_index);
    }
        return true;
    case GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin:
    case GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineAfterEnd:
    case GeometryUtilities::PointSegmentPositionTypes::LeftTheSegment:
    case GeometryUtilities::PointSegmentPositionTypes::RightTheSegment:
        return false;
    default:
        throw std::runtime_error("segment face single intersection not supported");
    }
}
// ***************************************************************************
void IntersectorMesh3DSegment::InsertNewIntersection(const double &segment_intersection_coordinate,
                                                     const IMeshDAO &mesh3D,
                                                     const unsigned int cell3D_index,
                                                     const unsigned int cell2D_index,
                                                     std::map<double, IntersectionPoint> &mesh1D_intersections,
                                                     std::list<unsigned int> &cell3Ds_index) const
{
    bool found = false;
    auto &mesh1D_intersection = CreateOrFindNewIntersection(segment_intersection_coordinate, mesh1D_intersections, found);

    if (mesh1D_intersection.Cell3DIds.find(cell3D_index) == mesh1D_intersection.Cell3DIds.end())
    {
        mesh1D_intersection.Cell3DIds.insert(cell3D_index);
    }

    if (cell2D_index < mesh3D.Cell2DTotalNumber())
    {
        for (unsigned int face_3D_neigh = 0; face_3D_neigh < mesh3D.Cell2DNumberNeighbourCell3D(cell2D_index); face_3D_neigh++)
        {
            if (!mesh3D.Cell2DHasNeighbourCell3D(cell2D_index, face_3D_neigh))
                continue;

            cell3Ds_index.push_back(mesh3D.Cell2DNeighbourCell3D(cell2D_index, face_3D_neigh));
        }
    }
}
// ***************************************************************************
std::vector<IntersectorMesh3DSegment::IntersectionMesh::IntersectionMeshSegment> IntersectorMesh3DSegment::CreateIntersectionSegments(
    const Gedim::IMeshDAO &mesh3D,
    const std::vector<IntersectionMesh::IntersectionMeshPoint> &mesh1D_points) const
{
    if (mesh1D_points.size() == 0)
        return {};

    Gedim::Output::Assert(mesh1D_points.size() >= 2);

    std::vector<IntersectorMesh3DSegment::IntersectionMesh::IntersectionMeshSegment> mesh1D_segments;

    mesh1D_segments.resize(mesh1D_points.size() - 1);

    std::vector<std::set<unsigned int>> points_cell3Ds(mesh1D_points.size());
    for (unsigned int s = 0; s < mesh1D_points.size(); s++)
    {
        const std::vector<unsigned int> &cell3Ds = mesh1D_points[s].Cell3DIds;
        points_cell3Ds[s] = std::set<unsigned int>(cell3Ds.begin(), cell3Ds.end());
    }

    for (unsigned int s = 0; s < mesh1D_segments.size(); s++)
    {
        const std::vector<unsigned int> &origin_cell3Ds = mesh1D_points[s].Cell3DIds;
        const std::vector<unsigned int> &end_cell3Ds = mesh1D_points[s + 1].Cell3DIds;

        std::vector<unsigned int> segment_cell3Ds;
        std::set_intersection(origin_cell3Ds.begin(), origin_cell3Ds.end(), end_cell3Ds.begin(), end_cell3Ds.end(), std::back_inserter(segment_cell3Ds));

        std::vector<GeometryUtilities::SegmentPolyhedronPositionResult> Positions(segment_cell3Ds.size());

        for (unsigned int c = 0; c < segment_cell3Ds.size(); c++)
        {
            const auto it_origin = std::find(origin_cell3Ds.begin(), origin_cell3Ds.end(), segment_cell3Ds[c]);
            assert(it_origin != origin_cell3Ds.end());
            const unsigned int local_id_cell3D_origin = it_origin - origin_cell3Ds.begin();

            const auto it_end = std::find(end_cell3Ds.begin(), end_cell3Ds.end(), segment_cell3Ds[c]);
            assert(it_end != end_cell3Ds.end());
            const unsigned int local_id_cell3D_end = it_end - end_cell3Ds.begin();

            if (mesh1D_points[s].Positions[local_id_cell3D_origin].Type ==
                    GeometryUtilities::PointPolyhedronPositionResult::Types::Inside ||
                mesh1D_points[s + 1].Positions[local_id_cell3D_end].Type ==
                    GeometryUtilities::PointPolyhedronPositionResult::Types::Inside)
            {
                Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::Inside;
                Positions[c].BorderIndex = 0;
                continue;
            }

            if (mesh1D_points[s].Positions[local_id_cell3D_origin].Type ==
                    GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace ||
                mesh1D_points[s + 1].Positions[local_id_cell3D_end].Type ==
                    GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace)
            {
                GeometryUtilities::PointPolyhedronPositionResult first_position;
                GeometryUtilities::PointPolyhedronPositionResult second_position;

                if (mesh1D_points[s].Positions[local_id_cell3D_origin].Type ==
                    GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace)
                {
                    first_position = mesh1D_points[s].Positions[local_id_cell3D_origin];
                    second_position = mesh1D_points[s + 1].Positions[local_id_cell3D_end];
                }
                else
                {
                    first_position = mesh1D_points[s + 1].Positions[local_id_cell3D_end];
                    second_position = mesh1D_points[s].Positions[local_id_cell3D_origin];
                }

                if (second_position.Type == GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace)
                {
                    if (second_position.BorderIndex == first_position.BorderIndex)
                    {
                        Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::BorderFace;
                        Positions[c].BorderIndex = first_position.BorderIndex;
                    }
                    else
                    {
                        Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::Inside;
                        Positions[c].BorderIndex = 0;
                    }
                }
                else if (second_position.Type == GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge)
                {
                    const unsigned int id_cell_2D = mesh3D.Cell3DFace(segment_cell3Ds[c], first_position.BorderIndex);
                    const unsigned int id_cell_1D = mesh3D.Cell3DEdge(segment_cell3Ds[c], second_position.BorderIndex);
                    const unsigned int cell_2D_id_edge = mesh3D.Cell2DFindEdge(id_cell_2D, id_cell_1D);
                    if (cell_2D_id_edge == mesh3D.Cell2DNumberEdges(id_cell_2D))
                    {
                        Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::Inside;
                        Positions[c].BorderIndex = 0;
                    }
                    else
                    {
                        Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::BorderFace;
                        Positions[c].BorderIndex = first_position.BorderIndex;
                    }
                }
                else if (second_position.Type == GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex)
                {
                    const unsigned int id_cell_2D = mesh3D.Cell3DFace(segment_cell3Ds[c], first_position.BorderIndex);
                    const unsigned int id_cell_0D = mesh3D.Cell3DEdge(segment_cell3Ds[c], second_position.BorderIndex);
                    const unsigned int cell_2D_id_vertex = mesh3D.Cell2DFindVertex(id_cell_2D, id_cell_0D);
                    if (cell_2D_id_vertex == mesh3D.Cell2DNumberVertices(id_cell_2D))
                    {
                        Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::Inside;
                        Positions[c].BorderIndex = 0;
                    }
                    else
                    {
                        Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::BorderFace;
                        Positions[c].BorderIndex = first_position.BorderIndex;
                    }
                }
                else
                    throw std::runtime_error("Not valid Point-Polyhedron Position");
            }
            else if (mesh1D_points[s].Positions[local_id_cell3D_origin].Type ==
                         GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge ||
                     mesh1D_points[s + 1].Positions[local_id_cell3D_end].Type ==
                         GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge)
            {
                GeometryUtilities::PointPolyhedronPositionResult first_position;
                GeometryUtilities::PointPolyhedronPositionResult second_position;

                if (mesh1D_points[s].Positions[local_id_cell3D_origin].Type ==
                    GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge)
                {
                    first_position = mesh1D_points[s].Positions[local_id_cell3D_origin];
                    second_position = mesh1D_points[s + 1].Positions[local_id_cell3D_end];
                }
                else
                {
                    first_position = mesh1D_points[s + 1].Positions[local_id_cell3D_end];
                    second_position = mesh1D_points[s].Positions[local_id_cell3D_origin];
                }

                if (second_position.Type == GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge)
                {
                    if (second_position.BorderIndex == first_position.BorderIndex)
                    {
                        Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::BorderEdge;
                        Positions[c].BorderIndex = first_position.BorderIndex;
                    }
                    else
                    {
                        const unsigned int id_cell_1D_first = mesh3D.Cell3DEdge(segment_cell3Ds[c], first_position.BorderIndex);
                        const unsigned int id_cell_1D_second = mesh3D.Cell3DEdge(segment_cell3Ds[c], second_position.BorderIndex);
                        std::vector<unsigned int> neighbours_first = mesh3D.Cell1DNeighbourCell2Ds(id_cell_1D_first);
                        std::vector<unsigned int> neighbours_second = mesh3D.Cell1DNeighbourCell2Ds(id_cell_1D_second);

                        std::vector<unsigned int> all_segment_cell2Ds;
                        std::set_intersection(neighbours_first.begin(),
                                              neighbours_first.end(),
                                              neighbours_second.begin(),
                                              neighbours_second.end(),
                                              std::back_inserter(all_segment_cell2Ds));

                        std::vector<unsigned int> segment_cell2Ds;
                        std::vector<unsigned int> faces = mesh3D.Cell3DFaces(segment_cell3Ds[c]);
                        std::sort(faces.begin(), faces.end());
                        std::set_intersection(faces.begin(),
                                              faces.end(),
                                              all_segment_cell2Ds.begin(),
                                              all_segment_cell2Ds.end(),
                                              std::back_inserter(segment_cell2Ds));

                        assert(segment_cell2Ds.size() <= 1);

                        if (segment_cell2Ds.size() == 1)
                        {
                            const unsigned int local_id_cell_2D = mesh3D.Cell3DFindFace(segment_cell3Ds[c], segment_cell2Ds[0]);
                            assert(local_id_cell_2D != mesh3D.Cell3DNumberFaces(segment_cell3Ds[c]));
                            Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::BorderFace;
                            Positions[c].BorderIndex = local_id_cell_2D;
                        }
                        else
                        {
                            Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::Inside;
                            Positions[c].BorderIndex = 0;
                        }
                    }
                }
                else if (second_position.Type == GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex)
                {
                    const unsigned int id_cell_1D = mesh3D.Cell3DEdge(segment_cell3Ds[c], first_position.BorderIndex);
                    const unsigned int id_cell_0D = mesh3D.Cell3DVertex(segment_cell3Ds[c], second_position.BorderIndex);
                    const unsigned int cell_1D_id_vertex = mesh3D.Cell1DFindExtreme(id_cell_1D, id_cell_0D);

                    if (cell_1D_id_vertex == 2)
                    {
                        std::vector<unsigned int> neighbours_first = mesh3D.Cell1DNeighbourCell2Ds(id_cell_1D);
                        std::vector<unsigned int> neighbours_second = mesh3D.Cell0DNeighbourCell2Ds(id_cell_0D);

                        std::vector<unsigned int> all_segment_cell2Ds;
                        std::set_intersection(neighbours_first.begin(),
                                              neighbours_first.end(),
                                              neighbours_second.begin(),
                                              neighbours_second.end(),
                                              std::back_inserter(all_segment_cell2Ds));

                        std::vector<unsigned int> segment_cell2Ds;
                        std::vector<unsigned int> faces = mesh3D.Cell3DFaces(segment_cell3Ds[c]);
                        std::sort(faces.begin(), faces.end());
                        std::set_intersection(faces.begin(),
                                              faces.end(),
                                              all_segment_cell2Ds.begin(),
                                              all_segment_cell2Ds.end(),
                                              std::back_inserter(segment_cell2Ds));

                        // Possibilmente falso nel caso di facce allineate: not considered case
                        assert(segment_cell2Ds.size() <= 1);

                        if (segment_cell2Ds.size() == 1)
                        {
                            const unsigned int local_id_cell_2D = mesh3D.Cell3DFindFace(segment_cell3Ds[c], segment_cell2Ds[0]);
                            assert(local_id_cell_2D != mesh3D.Cell3DNumberFaces(segment_cell3Ds[c]));
                            Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::BorderFace;
                            Positions[c].BorderIndex = local_id_cell_2D;
                        }
                        else
                        {
                            Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::Inside;
                            Positions[c].BorderIndex = 0;
                        }
                    }
                    else
                    {
                        Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::BorderEdge;
                        Positions[c].BorderIndex = first_position.BorderIndex;
                    }
                }
                else
                    throw std::runtime_error("Not valid Point-Polyhedron Position");
            }
            else if (mesh1D_points[s].Positions[local_id_cell3D_origin].Type ==
                         GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex &&
                     mesh1D_points[s + 1].Positions[local_id_cell3D_end].Type ==
                         GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex)
            {

                GeometryUtilities::PointPolyhedronPositionResult first_position = mesh1D_points[s].Positions[local_id_cell3D_origin];
                GeometryUtilities::PointPolyhedronPositionResult second_position =
                    mesh1D_points[s + 1].Positions[local_id_cell3D_end];

                const unsigned int id_cell_0D_first = mesh3D.Cell3DVertex(segment_cell3Ds[c], first_position.BorderIndex);
                const unsigned int id_cell_0D_second = mesh3D.Cell3DVertex(segment_cell3Ds[c], second_position.BorderIndex);
                std::vector<unsigned int> neighbours_1D_first = mesh3D.Cell0DNeighbourCell1Ds(id_cell_0D_first);
                std::vector<unsigned int> neighbours_1D_second = mesh3D.Cell0DNeighbourCell1Ds(id_cell_0D_second);

                std::vector<unsigned int> all_segment_cell1Ds;
                std::set_intersection(neighbours_1D_first.begin(),
                                      neighbours_1D_first.end(),
                                      neighbours_1D_second.begin(),
                                      neighbours_1D_second.end(),
                                      std::back_inserter(all_segment_cell1Ds));

                std::vector<unsigned int> segment_cell1Ds;
                std::vector<unsigned int> edges = mesh3D.Cell3DEdges(segment_cell3Ds[c]);
                std::sort(edges.begin(), edges.end());
                std::set_intersection(edges.begin(),
                                      edges.end(),
                                      all_segment_cell1Ds.begin(),
                                      all_segment_cell1Ds.end(),
                                      std::back_inserter(segment_cell1Ds));

                assert(segment_cell1Ds.size() <= 1);

                if (segment_cell1Ds.size() == 0)
                {
                    std::vector<unsigned int> neighbours_2D_first = mesh3D.Cell0DNeighbourCell2Ds(id_cell_0D_first);
                    std::vector<unsigned int> neighbours_2D_second = mesh3D.Cell0DNeighbourCell2Ds(id_cell_0D_second);

                    std::vector<unsigned int> all_segment_cell2Ds;
                    std::set_intersection(neighbours_2D_first.begin(),
                                          neighbours_2D_first.end(),
                                          neighbours_2D_second.begin(),
                                          neighbours_2D_second.end(),
                                          std::back_inserter(all_segment_cell2Ds));

                    std::vector<unsigned int> segment_cell2Ds;
                    std::vector<unsigned int> faces = mesh3D.Cell3DFaces(segment_cell3Ds[c]);
                    std::sort(faces.begin(), faces.end());
                    std::set_intersection(faces.begin(),
                                          faces.end(),
                                          all_segment_cell2Ds.begin(),
                                          all_segment_cell2Ds.end(),
                                          std::back_inserter(segment_cell2Ds));

                    assert(segment_cell2Ds.size() <= 1);

                    if (segment_cell2Ds.size() == 1)
                    {
                        const unsigned int local_id_cell_2D = mesh3D.Cell3DFindFace(segment_cell3Ds[c], segment_cell2Ds[0]);
                        Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::BorderFace;
                        Positions[c].BorderIndex = local_id_cell_2D;
                    }
                    else
                    {
                        Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::Inside;
                        Positions[c].BorderIndex = 0;
                    }
                }
                else
                {
                    const unsigned int local_id_cell_1D = mesh3D.Cell3DFindEdge(segment_cell3Ds[c], segment_cell1Ds[0]);
                    Positions[c].Type = GeometryUtilities::SegmentPolyhedronPositionResult::Types::BorderEdge;
                    Positions[c].BorderIndex = local_id_cell_1D;
                }
            }
            else
                throw std::runtime_error("Not valid Point-Polyhedron Position");
        }

        mesh1D_segments[s] = {{s, s + 1}, segment_cell3Ds, Positions};
    }

    return mesh1D_segments;
}
// ***************************************************************************
IntersectorMesh3DSegment::FindSegmentStartingCell3DResult IntersectorMesh3DSegment::FindSegmentStartingCell3D(
    const Eigen::Vector3d &segmentOrigin,
    const Eigen::Vector3d &segmentEnd,
    const Gedim::IMeshDAO &mesh3D,
    const MeshUtilities::MeshGeometricData3D &mesh3D_geometricData) const
{
    FindSegmentStartingCell3DResult result;
    result.StartingCell3DIndex = mesh3D.Cell3DTotalNumber();
    result.SegmentFullInsideCell3D = false;

    for (unsigned int c3D_index = 0; c3D_index < mesh3D.Cell3DTotalNumber(); c3D_index++)
    {
        if (result.StartingCell3DIndex < mesh3D.Cell3DTotalNumber())
            break;

        if (!mesh3D.Cell3DIsActive(c3D_index))
            continue;

        const auto &cell3D_vertices = mesh3D_geometricData.Cell3DsVertices.at(c3D_index);
        const auto &cell3D_faces = mesh3D_geometricData.Cell3DsFaces.at(c3D_index);
        const auto &cell3D_faces_3D_vertices = mesh3D_geometricData.Cell3DsFaces3DVertices.at(c3D_index);
        const auto &cell3D_faces_2D_vertices = mesh3D_geometricData.Cell3DsFaces2DVertices.at(c3D_index);
        const auto &cell3D_faces_normal = mesh3D_geometricData.Cell3DsFacesNormals.at(c3D_index);
        const auto &cell3D_faces_normal_direction = mesh3D_geometricData.Cell3DsFacesNormalDirections.at(c3D_index);
        const auto &cell3D_faces_translation = mesh3D_geometricData.Cell3DsFacesTranslations.at(c3D_index);
        const auto &cell3D_faces_rotation = mesh3D_geometricData.Cell3DsFacesRotationMatrices.at(c3D_index);
        const auto &cell3D_tetrahedrons = mesh3D_geometricData.Cell3DsTetrahedronPoints.at(c3D_index);

        const Eigen::MatrixXd cell3D_BoundingBox = geometryUtilities.PointsBoundingBox(cell3D_vertices);

        if (!geometryUtilities.IsPointInBoundingBox(segmentOrigin, cell3D_BoundingBox))
            continue;

        const GeometryUtilities::PointPolyhedronPositionResult segment_origin_position =
            geometryUtilities.PointPolyhedronPosition(segmentOrigin,
                                                      cell3D_faces,
                                                      cell3D_faces_3D_vertices,
                                                      cell3D_faces_2D_vertices,
                                                      cell3D_faces_normal,
                                                      cell3D_faces_normal_direction,
                                                      cell3D_faces_translation,
                                                      cell3D_faces_rotation,
                                                      cell3D_tetrahedrons);

        switch (segment_origin_position.Type)
        {
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Outside:
            continue;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Inside: {
            result.StartingCell3DIndex = c3D_index;
            result.Position = segment_origin_position;

            if (!geometryUtilities.IsPointInBoundingBox(segmentEnd, cell3D_BoundingBox))
                break;

            const GeometryUtilities::PointPolyhedronPositionResult segment_end_position =
                geometryUtilities.PointPolyhedronPosition(segmentEnd,
                                                          cell3D_faces,
                                                          cell3D_faces_3D_vertices,
                                                          cell3D_faces_2D_vertices,
                                                          cell3D_faces_normal,
                                                          cell3D_faces_normal_direction,
                                                          cell3D_faces_translation,
                                                          cell3D_faces_rotation,
                                                          cell3D_tetrahedrons);

            switch (segment_end_position.Type)
            {
            case GeometryUtilities::PointPolyhedronPositionResult::Types::Inside: {
                result.SegmentFullInsideCell3D = true;
                break;
            }
            break;
            case GeometryUtilities::PointPolyhedronPositionResult::Types::Outside:
            case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace:
            case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
            case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex:
                break;
            default:
                throw std::runtime_error("Segment end position not supported");
            }
        }
        break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace:
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex:
            result.StartingCell3DIndex = c3D_index;
            result.Position = segment_origin_position;
            break;
        default:
            throw std::runtime_error("Segment origin position not supported");
        }
    }

    Gedim::Output::Assert(result.StartingCell3DIndex < mesh3D.Cell3DTotalNumber());

    return result;
}
// ***************************************************************************
std::vector<IntersectorMesh3DSegment::IntersectionMesh::IntersectionMeshPoint> IntersectorMesh3DSegment::FindSegmentIntersectionPoints(
    const Eigen::Vector3d &segmentOrigin,
    const Eigen::Vector3d &segmentEnd,
    const Eigen::Vector3d &segmentTangent,
    const Gedim::IMeshDAO &mesh3D,
    const Gedim::MeshUtilities::MeshGeometricData3D &mesh3D_geometricData,
    const unsigned int starting_cell3D_index) const
{
    std::map<double, IntersectionPoint> mesh1D_intersections;

    bool found = false;
    auto &starting_intersection = CreateOrFindNewIntersection(0.0, mesh1D_intersections, found);
    starting_intersection.Cell3DIds.insert(starting_cell3D_index);

    std::list<unsigned int> cell3Ds_index;
    cell3Ds_index.push_back(starting_cell3D_index);
    std::vector<bool> visited(mesh3D.Cell3DTotalNumber(), false);

    while (!cell3Ds_index.empty())
    {
        const unsigned int cell3D_index = cell3Ds_index.front();
        cell3Ds_index.pop_front();

        if (visited.at(cell3D_index))
            continue;

        visited[cell3D_index] = true;

        SegmentCell3DIntersection(segmentOrigin, segmentEnd, segmentTangent, mesh3D, mesh3D_geometricData, cell3D_index, mesh1D_intersections, cell3Ds_index);
    }

    std::vector<IntersectorMesh3DSegment::IntersectionMesh::IntersectionMeshPoint> result(mesh1D_intersections.size());

    unsigned int intersection_index = 0;
    for (const auto &intersection : mesh1D_intersections)
    {
        unsigned int c = 0;
        std::vector<GeometryUtilities::PointPolyhedronPositionResult> segment_origin_position(
            intersection.second.Cell3DIds.size());
        for (unsigned int c3D_index : intersection.second.Cell3DIds)
        {
            const auto &cell3D_faces = mesh3D_geometricData.Cell3DsFaces.at(c3D_index);
            const auto &cell3D_faces_3D_vertices = mesh3D_geometricData.Cell3DsFaces3DVertices.at(c3D_index);
            const auto &cell3D_faces_2D_vertices = mesh3D_geometricData.Cell3DsFaces2DVertices.at(c3D_index);
            const auto &cell3D_faces_normal = mesh3D_geometricData.Cell3DsFacesNormals.at(c3D_index);
            const auto &cell3D_faces_normal_direction = mesh3D_geometricData.Cell3DsFacesNormalDirections.at(c3D_index);
            const auto &cell3D_faces_translation = mesh3D_geometricData.Cell3DsFacesTranslations.at(c3D_index);
            const auto &cell3D_faces_rotation = mesh3D_geometricData.Cell3DsFacesRotationMatrices.at(c3D_index);
            const auto &cell3D_tetrahedrons = mesh3D_geometricData.Cell3DsTetrahedronPoints.at(c3D_index);

            segment_origin_position[c++] =
                geometryUtilities.PointPolyhedronPosition(segmentOrigin + intersection.first * (segmentEnd - segmentOrigin),
                                                          cell3D_faces,
                                                          cell3D_faces_3D_vertices,
                                                          cell3D_faces_2D_vertices,
                                                          cell3D_faces_normal,
                                                          cell3D_faces_normal_direction,
                                                          cell3D_faces_translation,
                                                          cell3D_faces_rotation,
                                                          cell3D_tetrahedrons);
        }

        result[intersection_index] = {
            intersection.first,
            std::vector<unsigned int>(intersection.second.Cell3DIds.begin(), intersection.second.Cell3DIds.end()),
            segment_origin_position};

        intersection_index++;
    }

    return result;
}
// ***************************************************************************
void IntersectorMesh3DSegment::SegmentCell3DIntersection(const Eigen::Vector3d &segmentOrigin,
                                                         const Eigen::Vector3d &segmentEnd,
                                                         const Eigen::Vector3d &segmentTangent,
                                                         const Gedim::IMeshDAO &mesh3D,
                                                         const Gedim::MeshUtilities::MeshGeometricData3D &mesh3D_geometricData,
                                                         const unsigned int cell3D_index,
                                                         std::map<double, IntersectionPoint> &mesh1D_intersections,
                                                         std::list<unsigned int> &cell3Ds_index) const
{
    const auto &cell3D_vertices = mesh3D_geometricData.Cell3DsVertices.at(cell3D_index);
    const auto &cell3D_faces = mesh3D_geometricData.Cell3DsFaces.at(cell3D_index);
    const auto &cell3D_faces_3D_vertices = mesh3D_geometricData.Cell3DsFaces3DVertices.at(cell3D_index);
    const auto &cell3D_faces_normal = mesh3D_geometricData.Cell3DsFacesNormals.at(cell3D_index);
    const auto &cell3D_faces_normal_direction = mesh3D_geometricData.Cell3DsFacesNormalDirections.at(cell3D_index);
    const auto &cell3D_faces_3D_edges_tangent = mesh3D_geometricData.Cell3DsFacesEdge3DTangents.at(cell3D_index);
    const auto &cell3D_faces_edges_direction = mesh3D_geometricData.Cell3DsFacesEdgeDirections.at(cell3D_index);
    const auto &cell3D_faces_2D_vertices = mesh3D_geometricData.Cell3DsFaces2DVertices.at(cell3D_index);
    const auto &cell3D_faces_translation = mesh3D_geometricData.Cell3DsFacesTranslations.at(cell3D_index);
    const auto &cell3D_faces_rotation = mesh3D_geometricData.Cell3DsFacesRotationMatrices.at(cell3D_index);
    const auto &cell3D_tetrahedrons = mesh3D_geometricData.Cell3DsTetrahedronPoints.at(cell3D_index);

    const Eigen::MatrixXd cell3D_BoundingBox = geometryUtilities.PointsBoundingBox(cell3D_vertices);

    CheckSegmentPolyhedronIntersection(0.0,
                                       segmentOrigin,
                                       mesh3D,
                                       cell3D_index,
                                       cell3D_BoundingBox,
                                       cell3D_faces,
                                       cell3D_faces_3D_vertices,
                                       cell3D_faces_2D_vertices,
                                       cell3D_faces_normal,
                                       cell3D_faces_normal_direction,
                                       cell3D_faces_translation,
                                       cell3D_faces_rotation,
                                       cell3D_tetrahedrons,
                                       mesh1D_intersections,
                                       cell3Ds_index);

    CheckSegmentPolyhedronIntersection(1.0,
                                       segmentEnd,
                                       mesh3D,
                                       cell3D_index,
                                       cell3D_BoundingBox,
                                       cell3D_faces,
                                       cell3D_faces_3D_vertices,
                                       cell3D_faces_2D_vertices,
                                       cell3D_faces_normal,
                                       cell3D_faces_normal_direction,
                                       cell3D_faces_translation,
                                       cell3D_faces_rotation,
                                       cell3D_tetrahedrons,
                                       mesh1D_intersections,
                                       cell3Ds_index);

    for (unsigned int f = 0; f < cell3D_faces.size(); f++)
    {
        const unsigned int cell2D_index = mesh3D.Cell3DFace(cell3D_index, f);

        const Eigen::MatrixXd &face_3D_vertices = cell3D_faces_3D_vertices.at(f);
        const unsigned int face_num_edges = face_3D_vertices.cols();
        const Eigen::Vector3d &face_normal = cell3D_faces_normal.at(f);
        const Eigen::MatrixXd &face_3D_edges_tangent = cell3D_faces_3D_edges_tangent.at(f);
        const std::vector<bool> &face_edges_direction = cell3D_faces_edges_direction.at(f);
        const Eigen::MatrixXd &face_2D_vertices = cell3D_faces_2D_vertices.at(f);
        const Eigen::Vector3d &face_translation = cell3D_faces_translation.at(f);
        const Eigen::Matrix3d &face_rotation = cell3D_faces_rotation.at(f);

        const GeometryUtilities::IntersectionSegmentPlaneResult segment_face_plane_intersection =
            geometryUtilities.IntersectionSegmentPlane(segmentOrigin, segmentEnd, face_normal, face_3D_vertices.col(0));

        switch (segment_face_plane_intersection.Type)
        {
        case GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection:
            continue;
        case GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection: {
            const auto &segment_intersection = segment_face_plane_intersection.SingleIntersection;

            CheckSegmentFaceIntersection(segmentOrigin,
                                         segmentTangent,
                                         segment_intersection.Type,
                                         segment_intersection.CurvilinearCoordinate,
                                         face_rotation,
                                         face_translation,
                                         face_2D_vertices,
                                         mesh3D,
                                         cell3D_index,
                                         cell2D_index,
                                         mesh1D_intersections,
                                         cell3Ds_index);
        }
        break;
        case GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections: {
            CheckSegmentFaceIntersection(segmentOrigin,
                                         segmentTangent,
                                         Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin,
                                         0.0,
                                         face_rotation,
                                         face_translation,
                                         face_2D_vertices,
                                         mesh3D,
                                         cell3D_index,
                                         cell2D_index,
                                         mesh1D_intersections,
                                         cell3Ds_index);

            CheckSegmentFaceIntersection(segmentOrigin,
                                         segmentTangent,
                                         Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd,
                                         1.0,
                                         face_rotation,
                                         face_translation,
                                         face_2D_vertices,
                                         mesh3D,
                                         cell3D_index,
                                         cell2D_index,
                                         mesh1D_intersections,
                                         cell3Ds_index);

            for (unsigned int e = 0; e < face_num_edges; e++)
            {
                const double edge_direction = face_edges_direction.at(e) ? +1.0 : -1.0;
                const Eigen::Vector3d edge_origin =
                    face_edges_direction.at(e) ? face_3D_vertices.col(e) : face_3D_vertices.col((e + 1) % face_num_edges);
                const Eigen::Vector3d edge_tangent = edge_direction * face_3D_edges_tangent.col(e);

                const auto segment_edge_intersection =
                    geometryUtilities.IntersectionSegmentSegment(segmentOrigin, segmentEnd, edge_origin, edge_origin + edge_tangent);
                switch (segment_edge_intersection.IntersectionLinesType)
                {
                case GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarParallel:
                    continue;
                case GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarIntersecting: {
                    switch (segment_edge_intersection.IntersectionSegmentsType)
                    {
                    case GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection:
                        continue;
                    case GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection: {
                        const auto &segment_intersection = segment_edge_intersection.FirstSegmentIntersections.at(0);

                        CheckSegmentEdgeIntersection(segment_intersection.Type,
                                                     segment_intersection.CurvilinearCoordinate,
                                                     mesh3D,
                                                     cell3D_index,
                                                     cell2D_index,
                                                     mesh1D_intersections,
                                                     cell3Ds_index);
                    }
                    break;
                    case GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections: {
                        const auto &segment_first_intersection = segment_edge_intersection.FirstSegmentIntersections.at(0);
                        const auto &segment_second_intersection = segment_edge_intersection.FirstSegmentIntersections.at(1);

                        CheckSegmentEdgeIntersection(segment_first_intersection.Type,
                                                     segment_first_intersection.CurvilinearCoordinate,
                                                     mesh3D,
                                                     cell3D_index,
                                                     cell2D_index,
                                                     mesh1D_intersections,
                                                     cell3Ds_index);
                        CheckSegmentEdgeIntersection(segment_second_intersection.Type,
                                                     segment_second_intersection.CurvilinearCoordinate,
                                                     mesh3D,
                                                     cell3D_index,
                                                     cell2D_index,
                                                     mesh1D_intersections,
                                                     cell3Ds_index);
                    }
                    break;
                    default:
                        throw std::runtime_error("segment edge segment intersection not supported");
                    }
                }
                break;
                default:
                    throw std::runtime_error("segment edge line intersection not supported");
                }
            }
        }
            continue;
        default:
            throw std::runtime_error("segment face intersection not supported");
        }
    }
}
// ***************************************************************************
IntersectorMesh3DSegment::IntersectionMesh IntersectorMesh3DSegment::CreateIntersectionMesh(
    const Eigen::Vector3d &segmentOrigin,
    const Eigen::Vector3d &segmentEnd,
    const Eigen::Vector3d &segmentTangent,
    const Gedim::IMeshDAO &mesh3D,
    const Gedim::MeshUtilities::MeshGeometricData3D &mesh3D_geometricData) const
{
    const IntersectorMesh3DSegment::FindSegmentStartingCell3DResult segment_starting_cell3D =
        FindSegmentStartingCell3D(segmentOrigin, segmentEnd, mesh3D, mesh3D_geometricData);

    if (segment_starting_cell3D.SegmentFullInsideCell3D)
    {
        IntersectionMesh mesh1D;

        mesh1D.Points = {{0.0, {segment_starting_cell3D.StartingCell3DIndex}, {segment_starting_cell3D.Position}},
                         {1.0, {segment_starting_cell3D.StartingCell3DIndex}, {segment_starting_cell3D.Position}}};
        mesh1D.Segments = CreateIntersectionSegments(mesh3D, mesh1D.Points);

        return mesh1D;
    }

    IntersectionMesh mesh1D;
    mesh1D.Points = FindSegmentIntersectionPoints(segmentOrigin,
                                                  segmentEnd,
                                                  segmentTangent,
                                                  mesh3D,
                                                  mesh3D_geometricData,
                                                  segment_starting_cell3D.StartingCell3DIndex);

    mesh1D.Segments = CreateIntersectionSegments(mesh3D, mesh1D.Points);

    return mesh1D;
}
// ***************************************************************************
} // namespace Gedim
