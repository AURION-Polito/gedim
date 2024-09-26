#include "MeshUtilities.hpp"

#include "Eigen_Utilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  Gedim::MeshUtilities::Intersect_mesh_polyhedron_result
  Gedim::MeshUtilities::Intersect_mesh_polyhedron(const Gedim::GeometryUtilities& geometry_utilities,
                                                  const Eigen::MatrixXd& polyhedron_vertices,
                                                  const Eigen::MatrixXi& polyhedron_edges,
                                                  const std::vector<Eigen::MatrixXd>& polyhedron_edges_vertices,
                                                  const Eigen::MatrixXd& polyhedron_edges_tangent,
                                                  const std::vector<Eigen::MatrixXd>& polyhedron_edges_boudingBox,
                                                  const std::vector<Eigen::MatrixXi>& polyhedron_faces,
                                                  const std::vector<Eigen::MatrixXd>& polyhedron_faces_vertices,
                                                  const std::vector<Eigen::MatrixXd>& polyhedron_faces_rotated_vertices,
                                                  const std::vector<Eigen::Vector3d>& polyhedron_faces_normals,
                                                  const std::vector<bool>& polyhedron_faces_normal_direction,
                                                  const std::vector<Eigen::Vector3d>& polyhedron_faces_translation,
                                                  const std::vector<Eigen::Matrix3d>& polyhedron_faces_rotation_matrix,
                                                  const std::vector<Eigen::MatrixXd>& polyhedron_faces_boudingBox,
                                                  const Eigen::MatrixXd& polyhedron_boudingBox,
                                                  const IMeshDAO& mesh,
                                                  const std::vector<Eigen::MatrixXd>& mesh_cell1Ds_boudingBox,
                                                  const std::vector<Eigen::MatrixXd>& mesh_cell1Ds_vertices,
                                                  const std::vector<Eigen::Vector3d>& mesh_cell1Ds_tangent,
                                                  const std::vector<Eigen::MatrixXd>& mesh_cell2Ds_vertices,
                                                  const std::vector<Eigen::Vector3d>& mesh_cell2Ds_normal,
                                                  const std::vector<Eigen::MatrixXd>& mesh_cell2Ds_2D_vertices,
                                                  const std::vector<Eigen::Vector3d>& mesh_cell2Ds_translation,
                                                  const std::vector<Eigen::Matrix3d>& mesh_cell2Ds_rotation_matrix,
                                                  const std::vector<Eigen::MatrixXd>& mesh_cell2Ds_boudingBox,
                                                  const std::vector<Eigen::MatrixXd>& mesh_cell3Ds_boudingBox,
                                                  const std::vector<std::vector<Eigen::MatrixXi>>& mesh_cell3Ds_faces,
                                                  const std::vector<std::vector<Eigen::MatrixXd>>& mesh_cell3Ds_faces_vertices,
                                                  const std::vector<std::vector<Eigen::MatrixXd>>& mesh_cell3Ds_faces_2D_vertices,
                                                  const std::vector<std::vector<Eigen::Vector3d>>& mesh_cell3Ds_faces_normal,
                                                  const std::vector<std::vector<bool>>& mesh_cell3Ds_faces_normal_directions,
                                                  const std::vector<std::vector<Eigen::Vector3d>>& mesh_cell3Ds_faces_translation,
                                                  const std::vector<std::vector<Eigen::Matrix3d>>& mesh_cell3Ds_faces_rotation_matrix) const
  {
    struct Intersection final
    {
        unsigned int Intersection_index;
        Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types Type;
        unsigned int Geometry_index;
        std::list<unsigned int> Cell0Ds_index;
        std::set<unsigned int> Cell1Ds_index;
        std::set<unsigned int> Cell2Ds_index;
        std::list<unsigned int> Cell3Ds_index;
        Eigen::Vector3d Intersection_coordinates;
    };

    std::map<unsigned int, std::list<unsigned int>> cell1Ds_intersections;
    std::map<unsigned int, std::list<unsigned int>> cell2Ds_intersections;
    std::map<unsigned int, std::list<unsigned int>> cell3Ds_intersections;

    auto find_intersection = [&geometry_utilities](
                             std::list<Intersection>& intersections,
        Eigen::Vector3d new_intersection_coordinates)
    {
      std::list<Intersection>::iterator it_found = intersections.end();
      for (std::list<Intersection>::iterator it = intersections.begin();
           it != intersections.end();
           it++)
      {
        const Intersection& intersection = *it;

        if (!geometry_utilities.PointsAreCoincident(intersection.Intersection_coordinates,
                                                    new_intersection_coordinates))
          continue;

        it_found = it;
        break;
      }

      return it_found;
    };

    std::list<Intersection> intersections;
    Intersect_mesh_polyhedron_result::Mesh_Intersections mesh_intersections;

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); ++p)
    {
      if (!mesh.Cell0DIsActive(p))
        continue;

      const auto cell0D_coordinates = mesh.Cell0DCoordinates(p);

      if (!geometry_utilities.IsPointInBoundingBox(cell0D_coordinates,
                                                   polyhedron_boudingBox))
        continue;


      const auto cell0D_polyhedron_position = geometry_utilities.PointPolyhedronPosition(cell0D_coordinates,
                                                                                         polyhedron_faces,
                                                                                         polyhedron_faces_vertices,
                                                                                         polyhedron_faces_rotated_vertices,
                                                                                         polyhedron_faces_normals,
                                                                                         polyhedron_faces_normal_direction,
                                                                                         polyhedron_faces_translation,
                                                                                         polyhedron_faces_rotation_matrix);

      switch (cell0D_polyhedron_position.Type)
      {
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Outside:
          continue;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Inside:
        {
          unsigned int new_intersection_index = intersections.size();
          intersections.push_back({});

          Intersection& new_intersection = intersections.back();
          new_intersection.Intersection_index = new_intersection_index;
          new_intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Polyhedron;
          new_intersection.Geometry_index = 0;
          new_intersection.Intersection_coordinates = cell0D_coordinates;
          new_intersection.Cell0Ds_index = { p };

          mesh_intersections.Cell0Ds_intersections.insert(std::make_pair(p, new_intersection_index));
        }
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace:
        {
          unsigned int new_intersection_index = intersections.size();
          intersections.push_back({});

          Intersection& new_intersection = intersections.back();
          new_intersection.Intersection_index = new_intersection_index;
          new_intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Face;
          new_intersection.Geometry_index = cell0D_polyhedron_position.BorderIndex;
          new_intersection.Intersection_coordinates = cell0D_coordinates;
          new_intersection.Cell0Ds_index = { p };

          mesh_intersections.Cell0Ds_intersections.insert(std::make_pair(p, new_intersection_index));
        }
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
        {
          unsigned int new_intersection_index = intersections.size();
          intersections.push_back({});

          Intersection& new_intersection = intersections.back();
          new_intersection.Intersection_index = new_intersection_index;
          new_intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Edge;
          new_intersection.Geometry_index = cell0D_polyhedron_position.BorderIndex;
          new_intersection.Intersection_coordinates = cell0D_coordinates;
          new_intersection.Cell0Ds_index = { p };

          mesh_intersections.Cell0Ds_intersections.insert(std::make_pair(p, new_intersection_index));
        }
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex:
        {
          unsigned int new_intersection_index = intersections.size();
          intersections.push_back({});

          Intersection& new_intersection = intersections.back();
          new_intersection.Intersection_index = new_intersection_index;
          new_intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Vertex;
          new_intersection.Geometry_index = cell0D_polyhedron_position.BorderIndex;
          new_intersection.Intersection_coordinates = cell0D_coordinates;
          new_intersection.Cell0Ds_index = { p };

          mesh_intersections.Cell0Ds_intersections.insert(std::make_pair(p, new_intersection_index));
        }
          break;
        default:
          throw std::runtime_error("Unexpected cell0D intersection");
      }
    }

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
    {
      if (!mesh.Cell1DIsActive(e))
        continue;

      if (!geometry_utilities.BoundingBoxesIntersects(mesh_cell1Ds_boudingBox.at(e),
                                                      polyhedron_boudingBox))
        continue;

      const auto& cell1D_vertices = mesh_cell1Ds_vertices.at(e);
      const auto& cell1D_tangent = mesh_cell1Ds_tangent.at(e);

      for (unsigned int p_f = 0; p_f < polyhedron_faces.size(); ++p_f)
      {
        if (!geometry_utilities.BoundingBoxesIntersects(mesh_cell1Ds_boudingBox.at(e),
                                                        polyhedron_faces_boudingBox.at(p_f)))
          continue;

        const auto& face = polyhedron_faces.at(p_f);
        const auto& face_vertices = polyhedron_faces_vertices.at(p_f);
        const auto& face_2D_vertices = polyhedron_faces_rotated_vertices.at(p_f);
        const auto& face_normal = polyhedron_faces_normals.at(p_f);
        const auto& face_translation = polyhedron_faces_translation.at(p_f);
        const auto& face_rotation_matrix = polyhedron_faces_rotation_matrix.at(p_f);

        const auto intersection_face_plane_cell1D = geometry_utilities.IntersectionSegmentPlane(cell1D_vertices.col(0),
                                                                                                cell1D_vertices.col(1),
                                                                                                face_normal,
                                                                                                face_vertices.col(0));

        switch (intersection_face_plane_cell1D.Type)
        {
          case Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection:
          case Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections:
            continue;
          case Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection:
            break;
          default:
            throw std::runtime_error("Unexpected cell2D intersection");
        }

        switch (intersection_face_plane_cell1D.SingleIntersection.Type)
        {
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::LeftTheSegment:
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::RightTheSegment:
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineAfterEnd:
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin:
            continue;
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin:
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd:
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment:
            break;
          default:
            throw std::runtime_error("Unexpected cell2D edge intersection");
        }

        const Eigen::Vector3d intersection_coordinate = cell1D_vertices.col(0) +
                                                        intersection_face_plane_cell1D.SingleIntersection.CurvilinearCoordinate *
                                                        cell1D_tangent;

        const Eigen::Vector3d intersection_2D = geometry_utilities.RotatePointsFrom3DTo2D(intersection_coordinate,
                                                                                          face_rotation_matrix.transpose(),
                                                                                          face_translation);

        const auto intersection_face_position = geometry_utilities.PointPolygonPosition(intersection_2D,
                                                                                        face_2D_vertices);

        if (intersection_face_position.Type ==
            Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside)
          continue;

        auto intersection_found = find_intersection(intersections,
                                                    intersection_coordinate);

        const bool new_intersection = (intersection_found == intersections.end());

        if (new_intersection)
        {
          unsigned int new_intersection_index = intersections.size();
          intersections.push_back({});

          intersection_found = std::next(intersections.end(), -1);

          Intersection& new_intersection = *intersection_found;
          new_intersection.Intersection_index = new_intersection_index;
          new_intersection.Intersection_coordinates = intersection_coordinate;
        }

        Intersection& intersection = *intersection_found;
        if (intersection.Cell1Ds_index.find(e) == intersection.Cell1Ds_index.end())
          intersection.Cell1Ds_index.insert(e);
        if (cell1Ds_intersections.find(e) == cell1Ds_intersections.end())
          cell1Ds_intersections.insert(std::make_pair(e, std::list<unsigned int>()));
        cell1Ds_intersections[e].push_back(intersection.Intersection_index);

        if (!new_intersection)
          continue;

        switch (intersection_face_position.Type)
        {
          case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex:
          {
            intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Vertex;
            intersection.Geometry_index = face(0, intersection_face_position.BorderIndex);
          }
            break;
          case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge:
          {
            intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Edge;
            intersection.Geometry_index = face(1, intersection_face_position.BorderIndex);
          }
            break;
          case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside:
          {
            intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Face;
            intersection.Geometry_index = p_f;
          }
            break;
          default:
            throw std::runtime_error("Unexpected face intersection type");
        }
      }
    }

    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); ++f)
    {
      if (!mesh.Cell2DIsActive(f))
        continue;

      if (!geometry_utilities.BoundingBoxesIntersects(mesh_cell2Ds_boudingBox.at(f),
                                                      polyhedron_boudingBox))
        continue;

      const auto& cell2D_vertices = mesh_cell2Ds_vertices.at(f);
      const auto& cell2D_2D_vertices = mesh_cell2Ds_2D_vertices.at(f);
      const auto& cell2D_normal = mesh_cell2Ds_normal.at(f);
      const auto& cell2D_translation = mesh_cell2Ds_translation.at(f);
      const auto& cell2D_rotation_matrix = mesh_cell2Ds_rotation_matrix.at(f);

      for (unsigned int p_e = 0; p_e < polyhedron_edges.cols(); ++p_e)
      {
        if (!geometry_utilities.BoundingBoxesIntersects(mesh_cell2Ds_boudingBox.at(f),
                                                        polyhedron_edges_boudingBox.at(p_e)))
          continue;

        const auto& polyhedron_edge_vertices = polyhedron_edges_vertices.at(p_e);
        const auto& polyhedron_edge_tangent = polyhedron_edges_tangent.col(p_e);

        const auto intersection_cell2D_plane_edge = geometry_utilities.IntersectionSegmentPlane(polyhedron_edge_vertices.col(0),
                                                                                                polyhedron_edge_vertices.col(1),
                                                                                                cell2D_normal,
                                                                                                cell2D_vertices.col(0));

        switch (intersection_cell2D_plane_edge.Type)
        {
          case Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection:
          case Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections:
            continue;
          case Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection:
            break;
          default:
            throw std::runtime_error("Unexpected cell2D intersection");
        }

        switch (intersection_cell2D_plane_edge.SingleIntersection.Type)
        {
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::LeftTheSegment:
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::RightTheSegment:
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineAfterEnd:
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin:
            continue;
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin:
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd:
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment:
            break;
          default:
            throw std::runtime_error("Unexpected cell2D edge intersection");
        }

        const Eigen::Vector3d intersection_coordinate = polyhedron_edge_vertices.col(0) +
                                                        intersection_cell2D_plane_edge.SingleIntersection.CurvilinearCoordinate *
                                                        polyhedron_edge_tangent;

        const Eigen::Vector3d intersection_2D = geometry_utilities.RotatePointsFrom3DTo2D(intersection_coordinate,
                                                                                          cell2D_rotation_matrix.transpose(),
                                                                                          cell2D_translation);

        if (!geometry_utilities.IsPointInsidePolygon(intersection_2D,
                                                     cell2D_2D_vertices))
          continue;

        auto intersection_found = find_intersection(intersections,
                                                    intersection_coordinate);

        const bool new_intersection = (intersection_found == intersections.end());

        if (new_intersection)
        {
          unsigned int new_intersection_index = intersections.size();
          intersections.push_back({});

          intersection_found = std::next(intersections.end(), -1);

          Intersection& new_intersection = *intersection_found;
          new_intersection.Intersection_index = new_intersection_index;
          new_intersection.Intersection_coordinates = intersection_coordinate;
        }

        Intersection& intersection = *intersection_found;
        if (intersection.Cell2Ds_index.find(f) == intersection.Cell2Ds_index.end())
          intersection.Cell2Ds_index.insert(f);
        if (cell2Ds_intersections.find(f) == cell2Ds_intersections.end())
          cell2Ds_intersections.insert(std::make_pair(f, std::list<unsigned int>()));
        cell2Ds_intersections[f].push_back(intersection.Intersection_index);

        if (!new_intersection)
          continue;

        switch (intersection_cell2D_plane_edge.SingleIntersection.Type)
        {
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin:
          {
            intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Vertex;
            intersection.Geometry_index = polyhedron_edges(0, p_e);
          }
            break;
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd:
          {
            intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Vertex;
            intersection.Geometry_index = polyhedron_edges(1, p_e);
          }
            break;
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment:
          {
            intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Edge;
            intersection.Geometry_index = p_e;
          }
            break;
          default:
            throw std::runtime_error("Unexpected cell1D intersection type");
        }
      }
    }

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    {
      if (!mesh.Cell3DIsActive(c))
        continue;

      if (!geometry_utilities.BoundingBoxesIntersects(mesh_cell3Ds_boudingBox.at(c),
                                                      polyhedron_boudingBox))
        continue;

      const auto& cell3D_faces = mesh_cell3Ds_faces.at(c);
      const auto& cell3D_faces_vertices = mesh_cell3Ds_faces_vertices.at(c);
      const auto& cell3D_faces_2D_vertices = mesh_cell3Ds_faces_2D_vertices.at(c);
      const auto& cell3D_faces_normal = mesh_cell3Ds_faces_normal.at(c);
      const auto& cell3D_faces_normal_directions = mesh_cell3Ds_faces_normal_directions.at(c);
      const auto& cell3D_faces_translation = mesh_cell3Ds_faces_translation.at(c);
      const auto& cell3D_faces_rotation_matrix = mesh_cell3Ds_faces_rotation_matrix.at(c);

      for (unsigned int p_v = 0; p_v < polyhedron_vertices.cols(); ++p_v)
      {
        const Eigen::Vector3d vertex = polyhedron_vertices.col(p_v);

        if (!geometry_utilities.IsPointInBoundingBox(vertex,
                                                     mesh_cell3Ds_boudingBox.at(c)))
          continue;

        const auto point_cell3D_position = geometry_utilities.PointPolyhedronPosition(vertex,
                                                                                      cell3D_faces,
                                                                                      cell3D_faces_vertices,
                                                                                      cell3D_faces_2D_vertices,
                                                                                      cell3D_faces_normal,
                                                                                      cell3D_faces_normal_directions,
                                                                                      cell3D_faces_translation,
                                                                                      cell3D_faces_rotation_matrix);
        switch (point_cell3D_position.Type)
        {
          case GeometryUtilities::PointPolyhedronPositionResult::Types::Outside:
            continue;
          case GeometryUtilities::PointPolyhedronPositionResult::Types::Inside:
          case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace:
          case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
          case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex:
          {
            auto intersection_found = find_intersection(intersections,
                                                        vertex);

            const bool new_intersection = (intersection_found == intersections.end());

            if (new_intersection)
            {
              unsigned int new_intersection_index = intersections.size();
              intersections.push_back({});

              intersection_found = std::next(intersections.end(), -1);

              Intersection& new_intersection = *intersection_found;
              new_intersection.Intersection_index = new_intersection_index;
              new_intersection.Intersection_coordinates = vertex;
            }

            Intersection& intersection = *intersection_found;
            intersection.Cell3Ds_index.push_back(c);

            if (cell3Ds_intersections.find(c) == cell3Ds_intersections.end())
              cell3Ds_intersections.insert(std::make_pair(c, std::list<unsigned int>()));
            cell3Ds_intersections[c].push_back(intersection.Intersection_index);

            if (!new_intersection)
              continue;

            intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Polyhedron;
            intersection.Geometry_index = 0;
          }
            break;
          default:
            throw std::runtime_error("Unexpected cell3D intersection");
        }
      }
    }

    if (intersections.empty())
    {
      return
      {
        Intersect_mesh_polyhedron_result::Types::None,
        {},
        {},
        {}
      };
    }

    Intersect_mesh_polyhedron_result result;
    result.Type = Intersect_mesh_polyhedron_result::Types::Vertices;
    result.Mesh_intersections = mesh_intersections;

    result.Intersections_Coordinates.resize(3, intersections.size());
    result.Polyhedron_intersections.resize(intersections.size());

    unsigned int i = 0;
    for (const auto& intersection : intersections)
    {
      result.Intersections_Coordinates.col(i)<< intersection.Intersection_coordinates;

      result.Polyhedron_intersections[i].Type = intersection.Type;
      result.Polyhedron_intersections[i].Cell0Ds_index = std::vector<unsigned int>(intersection.Cell0Ds_index.begin(),
                                                                                   intersection.Cell0Ds_index.end());
      result.Polyhedron_intersections[i].Cell1Ds_index = std::vector<unsigned int>(intersection.Cell1Ds_index.begin(),
                                                                                   intersection.Cell1Ds_index.end());
      result.Polyhedron_intersections[i].Cell2Ds_index = std::vector<unsigned int>(intersection.Cell2Ds_index.begin(),
                                                                                   intersection.Cell2Ds_index.end());
      result.Polyhedron_intersections[i].Cell3Ds_index = std::vector<unsigned int>(intersection.Cell3Ds_index.begin(),
                                                                                   intersection.Cell3Ds_index.end());

      i++;
    }

    for (const auto& cell1D_intersections : cell1Ds_intersections)
    {
      result.Mesh_intersections.Cell1Ds_intersections.insert(std::make_pair(cell1D_intersections.first,
                                                                            std::vector<unsigned int>(cell1D_intersections.second.begin(),
                                                                                                      cell1D_intersections.second.end())));
    }

    for (const auto& cell2D_intersections : cell2Ds_intersections)
    {
      result.Mesh_intersections.Cell2Ds_intersections.insert(std::make_pair(cell2D_intersections.first,
                                                                            std::vector<unsigned int>(cell2D_intersections.second.begin(),
                                                                                                      cell2D_intersections.second.end())));
    }

    for (const auto& cell3D_intersections : cell3Ds_intersections)
    {
      result.Mesh_intersections.Cell3Ds_intersections.insert(std::make_pair(cell3D_intersections.first,
                                                                            std::vector<unsigned int>(cell3D_intersections.second.begin(),
                                                                                                      cell3D_intersections.second.end())));
    }

    return result;
  }
  // ***************************************************************************
}
