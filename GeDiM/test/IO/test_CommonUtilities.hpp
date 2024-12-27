#ifndef __TEST_COMMONUTILITIES_H
#define __TEST_COMMONUTILITIES_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "CommonUtilities.hpp"

namespace GedimUnitTesting
{
// ***************************************************************************
TEST(TestCommonUtilities, TestSortArrayIndices)
{
    const std::vector<unsigned int> sortIndices =
        Gedim::Utilities::SortArrayIndices(std::vector<unsigned int>{5, 6, 1, 4, 8, 0, 3, 2, 7, 9});

    EXPECT_EQ(std::vector<unsigned int>({5, 2, 7, 6, 3, 0, 1, 8, 4, 9}), sortIndices);
}
// ***************************************************************************
TEST(TestCommonUtilities, TestSlope)
{
    std::vector<double> x{6, 5, 11, 7, 5, 4, 4};
    std::vector<double> y{2, 3, 9, 1, 8, 7, 5};

    ASSERT_DOUBLE_EQ(0.3055555555555555, Gedim::Utilities::Slope(x.begin(), y.begin(), x.size()));

    ASSERT_DOUBLE_EQ(0.5833333333333333, Gedim::Utilities::Slope(x.begin() + 1, y.begin() + 1, 4));
}
// ***************************************************************************
} // namespace GedimUnitTesting

#endif // __TEST_COMMONUTILITIES_H
