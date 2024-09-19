#include "MeshUtilities.hpp"

#include "Eigen_Utilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  Gedim::MeshUtilities::Intersect_mesh_polyhedron_result
  Gedim::MeshUtilities::Intersect_mesh_polyhedron(const Gedim::GeometryUtilities& geometry_utilities,
                                                  const Eigen::MatrixXd& polyhedron_vertices,
                                                  const Eigen::MatrixXi& polyhedron_edges,
                                                  const std::vector<Eigen::MatrixXi>& polyhedron_faces,
                                                  const std::vector<Eigen::MatrixXd>& polyhedron_faces_vertices,
                                                  const std::vector<Eigen::MatrixXd>& polyhedron_faces_rotated_vertices,
                                                  const std::vector<Eigen::Vector3d>& polyhedron_faces_normals,
                                                  const std::vector<bool>& polyhedron_faces_normal_direction,
                                                  const std::vector<Eigen::Vector3d>& polyhedron_face_translations,
                                                  const std::vector<Eigen::Matrix3d>& polyhedron_face_rotation_matrix,
                                                  const Eigen::MatrixXd& polyhedron_boudingBox,
                                                  const IMeshDAO& mesh,
                                                  const std::vector<Eigen::MatrixXd>& mesh_cell1Ds_boudingBox,
                                                  const std::vector<Eigen::MatrixXd>& mesh_cell1Ds_vertices,
                                                  const std::vector<Eigen::Vector3d>& mesh_cell1Ds_tangent,
                                                  const std::vector<Eigen::MatrixXd>& mesh_cell2Ds_boudingBox,
                                                  const std::vector<Eigen::MatrixXd>& mesh_cell3Ds_boudingBox) const
  {
    struct Intersection final
    {
        Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types Type;
        unsigned int Geometry_index;
        std::list<unsigned int> Cell0Ds_index;
        std::list<unsigned int> Cell1Ds_index;
        std::list<unsigned int> Cell2Ds_index;
        std::list<unsigned int> Cell3Ds_index;
        Eigen::Vector3d Intersection_coordinates;
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
                                                                                         polyhedron_face_translations,
                                                                                         polyhedron_face_rotation_matrix);

      switch (cell0D_polyhedron_position.Type)
      {
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Outside:
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Inside:
        {
          unsigned int new_intersection_index = intersections.size();
          Intersection new_intersection;
          new_intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Polyhedron;
          new_intersection.Geometry_index = 0;
          new_intersection.Intersection_coordinates = cell0D_coordinates;
          new_intersection.Cell0Ds_index = { p };
          intersections.push_back(new_intersection);
          mesh_intersections.Cell0Ds_intersections.insert(std::make_pair(p, new_intersection_index));
        }
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace:
        {
          unsigned int new_intersection_index = intersections.size();
          Intersection new_intersection;
          new_intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Face;
          new_intersection.Geometry_index = cell0D_polyhedron_position.BorderIndex;
          new_intersection.Intersection_coordinates = cell0D_coordinates;
          new_intersection.Cell0Ds_index = { p };
          intersections.push_back(new_intersection);
          mesh_intersections.Cell0Ds_intersections.insert(std::make_pair(p, new_intersection_index));
        }
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
        {
          unsigned int new_intersection_index = intersections.size();
          Intersection new_intersection;
          new_intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Edge;
          new_intersection.Geometry_index = cell0D_polyhedron_position.BorderIndex;
          new_intersection.Intersection_coordinates = cell0D_coordinates;
          new_intersection.Cell0Ds_index = { p };
          intersections.push_back(new_intersection);
          mesh_intersections.Cell0Ds_intersections.insert(std::make_pair(p, new_intersection_index));
        }
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex:
        {
          unsigned int new_intersection_index = intersections.size();
          Intersection new_intersection;
          new_intersection.Type = Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types::Vertex;
          new_intersection.Geometry_index = cell0D_polyhedron_position.BorderIndex;
          new_intersection.Intersection_coordinates = cell0D_coordinates;
          new_intersection.Cell0Ds_index = { p };
          intersections.push_back(new_intersection);
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

      const auto intersection_polyhedron_line = geometry_utilities.IntersectionPolyhedronLine(polyhedron_vertices,
                                                                                              polyhedron_edges,
                                                                                              polyhedron_faces,
                                                                                              polyhedron_faces_normals,
                                                                                              polyhedron_faces_normal_direction,
                                                                                              mesh_cell1Ds_vertices.at(e).col(0),
                                                                                              mesh_cell1Ds_tangent.at(e));
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

    return result;
  }
  // ***************************************************************************
}
