#ifndef __TEST_GEOMETRY_POLYHEDRON_H
#define __TEST_GEOMETRY_POLYHEDRON_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronBarycenter)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube barycenter
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateCubeWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                               Eigen::Vector3d(0.0,1.0,0.0));
        ASSERT_TRUE(geometryUtility.PointsAreCoincident(geometryUtility.PolyhedronBarycenter(cube.Vertices),
                                                        Eigen::Vector3d(0.5, 0.5, 0.5)));
      }

      // check tetrahedron barycenter
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));

        ASSERT_TRUE(geometryUtility.PointsAreCoincident(geometryUtility.PolyhedronBarycenter(tetrahedron.Vertices),
                                                        Eigen::Vector3d(0.25, 0.25, 0.25)));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_POLYHEDRON_H
