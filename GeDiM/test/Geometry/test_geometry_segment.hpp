#ifndef __TEST_GEOMETRY_SEGMENT_H
#define __TEST_GEOMETRY_SEGMENT_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

  TEST(TestGeometryUtilities, TestSegmentLine)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check line passing in the origin positive
      {
        Eigen::Vector3d segmentOrigin(1.0, 1.0, 0.0);
        Eigen::Vector3d segmentEnd(2.0, 2.0, 0.0);

        const double expectedSlope = 1.0;
        const double expectedIntercept = 0.0;

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentSlope(segmentOrigin,
                                                                                    segmentEnd),
                                                     expectedSlope, geometryUtilities.Tolerance1D()));

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentIntercept(segmentOrigin,
                                                                                        segmentEnd),
                                                     expectedIntercept, geometryUtilities.Tolerance1D()));
      }

      // check line passing in the origin negative
      {
        Eigen::Vector3d segmentOrigin(-1.0, 1.0, 0.0);
        Eigen::Vector3d segmentEnd(2.0, -2.0, 0.0);

        const double expectedSlope = -1.0;
        const double expectedIntercept = 0.0;

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentSlope(segmentOrigin,
                                                                                    segmentEnd),
                                                     expectedSlope, geometryUtilities.Tolerance1D()));

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentIntercept(segmentOrigin,
                                                                                        segmentEnd),
                                                     expectedIntercept, geometryUtilities.Tolerance1D()));
      }

      // check line not passing in the origin positive
      {
        Eigen::Vector3d segmentOrigin(0.0, 1.0, 0.0);
        Eigen::Vector3d segmentEnd(-1.0, 0.0, 0.0);

        const double expectedSlope = 1.0;
        const double expectedIntercept = 1.0;

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentSlope(segmentOrigin,
                                                                                    segmentEnd),
                                                     expectedSlope, geometryUtilities.Tolerance1D()));

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentIntercept(segmentOrigin,
                                                                                        segmentEnd),
                                                     expectedIntercept, geometryUtilities.Tolerance1D()));
      }

      // check line not passing in the origin negtive
      {
        Eigen::Vector3d segmentOrigin(0.0, -1.0, 0.0);
        Eigen::Vector3d segmentEnd(-1.0, 0.0, 0.0);

        const double expectedSlope = -1.0;
        const double expectedIntercept = -1.0;

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentSlope(segmentOrigin,
                                                                                    segmentEnd),
                                                     expectedSlope, geometryUtilities.Tolerance1D()));

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentIntercept(segmentOrigin,
                                                                                        segmentEnd),
                                                     expectedIntercept, geometryUtilities.Tolerance1D()));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_SEGMENT_H
