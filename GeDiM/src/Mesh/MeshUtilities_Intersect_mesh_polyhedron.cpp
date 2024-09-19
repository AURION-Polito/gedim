#include "MeshUtilities.hpp"

#include "Eigen_Utilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  Gedim::MeshUtilities::Intersect_mesh_polyhedron_result
  Gedim::MeshUtilities::Intersect_mesh_polyhedron(const Gedim::GeometryUtilities& geometry_utilities,
                                                  const std::vector<Eigen::MatrixXi>& polyhedron_faces,
                                                  const std::vector<Eigen::MatrixXd>& polyhedron_faces_vertices,
                                                  const std::vector<Eigen::MatrixXd>& polyhedron_faces_rotated_vertices,
                                                  const std::vector<Eigen::Vector3d>& polyhedron_faces_normals,
                                                  const std::vector<bool>& polyhedron_faces_normal_direction,
                                                  const std::vector<Eigen::Vector3d>& polyhedron_face_translations,
                                                  const std::vector<Eigen::Matrix3d>& polyhedron_face_rotation_matrix,
                                                  const Eigen::MatrixXd& polyhedron_boudingBox,
                                                  const IMeshDAO& mesh) const
  {
    std::list<Eigen::Vector3d> intersections_coordinates;
    Intersect_mesh_polyhedron_result::Mesh_Intersections mesh_intersections;
    std::list<Intersect_mesh_polyhedron_result::Polyhedron_Intersection> polyhedron_intersections;

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
          break;

        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
          break;
        default:
          throw std::runtime_error("Unexpected cell0D intersection");
      }
    }

    if (intersections_coordinates.empty())
    {
      return
      {
        Intersect_mesh_polyhedron_result::Types::None,
        {},
        {},
        {}
      };
    }

    return
    {
      Intersect_mesh_polyhedron_result::Types::Vertices,
          Eigen_Utilities::CoordinatesToMatrix(intersections_coordinates),
          std::vector<Intersect_mesh_polyhedron_result::Polyhedron_Intersection>(polyhedron_intersections.begin(),
                                                                                 polyhedron_intersections.end()),
          mesh_intersections
    };
  }
  // ***************************************************************************
}
