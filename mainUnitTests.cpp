#include "test_geometry.hpp"
#include "test_mesh.hpp"
#include "test_intersectorMesh2DSegment.hpp"
#include "test_unionMeshSegment.hpp"
#include "test_conformMesh.hpp"
#include "test_export_conformMesh.hpp"

#include <gtest/gtest.h>

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
