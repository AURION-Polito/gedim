#ifndef __TEST_COMMONUTILITIES_H
#define __TEST_COMMONUTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "CommonUtilities.hpp"

namespace GedimUnitTesting
{
  // ***************************************************************************
  TEST(TestCommonUtilities, TestSortArrayIndices)
  {
    const std::vector<unsigned int> sortIndices = Gedim::Utilities::SortArrayIndices(std::vector<unsigned int> { 5,6,1,4,8,0,3,2,7,9 });


    EXPECT_EQ(std::vector<unsigned int>({ 5,2,7,6,3,0,1,8,4,9 }), sortIndices);
  }
  // ***************************************************************************
}

#endif // __TEST_COMMONUTILITIES_H
