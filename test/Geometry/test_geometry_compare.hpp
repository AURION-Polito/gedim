#ifndef __TEST_GEOMETRY_COMPARE_H
#define __TEST_GEOMETRY_COMPARE_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

  TEST(TestGeometryUtilities, TestCompareValues)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check Compare1DValues
      {
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, 1.0), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::Coincident);
        ASSERT_EQ(geometryUtility.Compare1DValues(1.0, -1.0), Gedim::GeometryUtilities::CompareTypes::SecondBeforeFirst);
      }

      // check Compare2DValues
      {
        ASSERT_EQ(geometryUtility.Compare2DValues(0.0, 1.0), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare2DValues(0.0, geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::Coincident);
        ASSERT_EQ(geometryUtility.Compare2DValues(1.0, -1.0), Gedim::GeometryUtilities::CompareTypes::SecondBeforeFirst);
      }

      // check Compare3DValues
      {
        ASSERT_EQ(geometryUtility.Compare3DValues(0.0, 1.0), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare3DValues(0.0, geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare3DValues(0.0, geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::Coincident);
        ASSERT_EQ(geometryUtility.Compare3DValues(1.0, -1.0), Gedim::GeometryUtilities::CompareTypes::SecondBeforeFirst);
      }

      // check IsLenghtPositive
      {
        ASSERT_FALSE(geometryUtility.IsValue1DPositive(0.0));
        ASSERT_FALSE(geometryUtility.IsValue1DPositive(geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue1DPositive(-1.0));
        ASSERT_TRUE(geometryUtility.IsValue1DPositive(2 * geometryUtilityConfig.Tolerance));
        ASSERT_TRUE(geometryUtility.IsValue1DPositive(10.0));
      }

      // check IsAreaPositive
      {
        ASSERT_FALSE(geometryUtility.IsValue2DPositive(0.0));
        ASSERT_TRUE(geometryUtility.IsValue2DPositive(geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue2DPositive(geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue2DPositive(-1.0));
        ASSERT_TRUE(geometryUtility.IsValue2DPositive(10.0));
      }

      // check IsLenghtPositive
      {
        ASSERT_FALSE(geometryUtility.IsValue3DPositive(0.0));
        ASSERT_TRUE(geometryUtility.IsValue3DPositive(geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue3DPositive(geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue3DPositive(-1.0));
        ASSERT_TRUE(geometryUtility.IsValue3DPositive(10.0));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_COMPARE_H
