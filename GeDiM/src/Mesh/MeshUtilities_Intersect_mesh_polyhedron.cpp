#include "MeshUtilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  Gedim::MeshUtilities::Intersect_mesh_polyhedron_result
  Gedim::MeshUtilities::Intersect_mesh_polyhedron(const Gedim::GeometryUtilities& geometry_utilities,
                                                  const Gedim::GeometryUtilities::Polyhedron& polyhedron,
                                                  const Eigen::MatrixXd& polyhedron_boudingBox,
                                                  const IMeshDAO& mesh) const
  {
    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); ++p)
    {
      if (!mesh.Cell0DIsActive(p))
        continue;

      const auto cell0D_coordinates = mesh.Cell0DCoordinates(p);

      if (!geometry_utilities.IsPointInBoundingBox(cell0D_coordinates,
                                                   polyhedron_boudingBox))
        continue;
    }

    return {};
  }
  // ***************************************************************************
}
