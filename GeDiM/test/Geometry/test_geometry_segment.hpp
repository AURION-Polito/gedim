#ifndef __TEST_GEOMETRY_SEGMENT_H
#define __TEST_GEOMETRY_SEGMENT_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

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

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentSlope(segmentOrigin, segmentEnd),
                                                         expectedSlope,
                                                         geometryUtilities.Tolerance1D()));

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentIntercept(segmentOrigin, segmentEnd),
                                                         expectedIntercept,
                                                         geometryUtilities.Tolerance1D()));
        }

        // check line passing in the origin negative
        {
            Eigen::Vector3d segmentOrigin(-1.0, 1.0, 0.0);
            Eigen::Vector3d segmentEnd(2.0, -2.0, 0.0);

            const double expectedSlope = -1.0;
            const double expectedIntercept = 0.0;

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentSlope(segmentOrigin, segmentEnd),
                                                         expectedSlope,
                                                         geometryUtilities.Tolerance1D()));

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentIntercept(segmentOrigin, segmentEnd),
                                                         expectedIntercept,
                                                         geometryUtilities.Tolerance1D()));
        }

        // check line not passing in the origin positive
        {
            Eigen::Vector3d segmentOrigin(0.0, 1.0, 0.0);
            Eigen::Vector3d segmentEnd(-1.0, 0.0, 0.0);

            const double expectedSlope = 1.0;
            const double expectedIntercept = 1.0;

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentSlope(segmentOrigin, segmentEnd),
                                                         expectedSlope,
                                                         geometryUtilities.Tolerance1D()));

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentIntercept(segmentOrigin, segmentEnd),
                                                         expectedIntercept,
                                                         geometryUtilities.Tolerance1D()));
        }

        // check line not passing in the origin negtive
        {
            Eigen::Vector3d segmentOrigin(0.0, -1.0, 0.0);
            Eigen::Vector3d segmentEnd(-1.0, 0.0, 0.0);

            const double expectedSlope = -1.0;
            const double expectedIntercept = -1.0;

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentSlope(segmentOrigin, segmentEnd),
                                                         expectedSlope,
                                                         geometryUtilities.Tolerance1D()));

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(geometryUtilities.SegmentIntercept(segmentOrigin, segmentEnd),
                                                         expectedIntercept,
                                                         geometryUtilities.Tolerance1D()));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, MakeConcatenation)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    {
        Eigen::MatrixXi segments;

        ASSERT_EQ(Eigen::MatrixXi(), geometryUtilities.MakeConcatenation(segments, 1));
    }

    {
        Eigen::MatrixXi segments(2, 1);
        segments.col(0) << 5, 4;

        Eigen::MatrixXi expected_result(2, 1);
        expected_result.col(0) << 4, 0;

        ASSERT_EQ(expected_result, geometryUtilities.MakeConcatenation(segments, 4));
    }

    {
        Eigen::MatrixXi segments(2, 5);
        segments.col(0) << 5, 4;
        segments.col(1) << 5, 3;
        segments.col(2) << 2, 4;
        segments.col(3) << 0, 1;
        segments.col(4) << 0, 3;

        ASSERT_EQ(Eigen::MatrixXi(), geometryUtilities.MakeConcatenation(segments, 12));
    }

    {
        Eigen::MatrixXi segments(2, 5);
        segments.col(0) << 5, 4;
        segments.col(1) << 5, 3;
        segments.col(2) << 2, 4;
        segments.col(3) << 0, 1;
        segments.col(4) << 0, 3;

        Eigen::MatrixXi expected_result(2, 5);
        expected_result.col(0) << 1, 3;
        expected_result.col(1) << 0, 4;
        expected_result.col(2) << 3, 1;
        expected_result.col(3) << 5, 0;
        expected_result.col(4) << 4, 2;

        ASSERT_EQ(expected_result, geometryUtilities.MakeConcatenation(segments, 1));
    }

    {
        Eigen::MatrixXi segments(2, 5);
        segments.col(0) << 6, 12;
        segments.col(1) << 6, 3;
        segments.col(2) << 2, 12;
        segments.col(3) << 0, 1;
        segments.col(4) << 0, 3;

        Eigen::MatrixXi expected_result(2, 5);
        expected_result.col(0) << 1, 3;
        expected_result.col(1) << 0, 4;
        expected_result.col(2) << 3, 1;
        expected_result.col(3) << 6, 0;
        expected_result.col(4) << 12, 2;

        ASSERT_EQ(expected_result, geometryUtilities.MakeConcatenation(segments, 1));
    }

    {
        Eigen::MatrixXi segments(2, 5);
        segments.col(0) << 6, 12;
        segments.col(1) << 6, 3;
        segments.col(2) << 2, 12;
        segments.col(3) << 0, 2;
        segments.col(4) << 0, 3;

        Eigen::MatrixXi expected_result(2, 5);
        expected_result.col(0) << 2, 2;
        expected_result.col(1) << 12, 0;
        expected_result.col(2) << 6, 1;
        expected_result.col(3) << 3, 4;
        expected_result.col(4) << 0, 3;
        ASSERT_EQ(expected_result, geometryUtilities.MakeConcatenation(segments, 2));
    }
}
} // namespace GedimUnitTesting

#endif // __TEST_GEOMETRY_SEGMENT_H
