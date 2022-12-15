#include "test_StringsUtilities.hpp"
#include "test_Configurations.hpp"
#include "test_geometry.hpp"
#include "test_geometry_intersection.hpp"
#include "test_geometry_point.hpp"
#include "test_geometry_segment.hpp"
#include "test_geometry_polygon.hpp"
#include "test_geometry_polyhedron.hpp"
#include "test_geometry_split.hpp"
#include "test_geometry_splitPolyhedron.hpp"
#include "test_LAPACK_utilities.hpp"
#include "test_Eigen.hpp"
#include "test_Pardiso.hpp"
#include "test_map.hpp"
#include "test_mesh.hpp"
#include "test_meshUtilities.hpp"
#include "test_intersectorMesh2DSegment.hpp"
#include "test_unionMeshSegment.hpp"
#include "test_conformMesh.hpp"
#include "test_export_conformMesh.hpp"
#include "test_VTKUtilities.hpp"
#include "test_Quadrature1D.hpp"
#include "test_Quadrature2D.hpp"
#include "test_Quadrature3D.hpp"
#include "test_conformMeshUtilities.hpp"

#include <gtest/gtest.h>

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
