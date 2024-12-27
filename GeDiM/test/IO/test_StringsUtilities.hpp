#ifndef __TEST_STRINGSUTILITIES_H
#define __TEST_STRINGSUTILITIES_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "StringsUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
TEST(TestStringsUtilities, TestStringSplit_Test)
{
    vector<string> splitString = Gedim::StringsUtilities::Split("test_pippo_due", '_');

    ASSERT_EQ(3, splitString.size());
    ASSERT_EQ("test", splitString[0]);
    ASSERT_EQ("pippo", splitString[1]);
    ASSERT_EQ("due", splitString[2]);

    vector<string> splitStringCharacters = Gedim::StringsUtilities::Split("test:pippo_due:tre", {'>', '_', ':', '<'});

    ASSERT_EQ(4, splitStringCharacters.size());
    ASSERT_EQ("test", splitStringCharacters[0]);
    ASSERT_EQ("pippo", splitStringCharacters[1]);
    ASSERT_EQ("due", splitStringCharacters[2]);
    ASSERT_EQ("tre", splitStringCharacters[3]);
}
// ***************************************************************************
TEST(TestStringsUtilities, TestFindSeparator_Test)
{
    ASSERT_EQ(':', Gedim::StringsUtilities::FindSeparator("id:type", "id", "type"));
    ASSERT_EQ('-', Gedim::StringsUtilities::FindSeparator("type-id", "id", "type"));
}
// ***************************************************************************
TEST(TestStringsUtilities, TestToUpper_Test)
{
    ASSERT_EQ("ASAWCF", Gedim::StringsUtilities::ToUpper("ASaWcF"));
}
// ***************************************************************************
TEST(TestStringsUtilities, TestToLower_Test)
{
    ASSERT_EQ("asawcf", Gedim::StringsUtilities::ToLower("ASaWcF"));
}
// ***************************************************************************
TEST(TestStringsUtilities, TestParserInt_Test)
{
    string stringValue = "1";
    int correctValue = 1;

    ASSERT_EQ(correctValue, Gedim::StringsUtilities::Parse<int>(stringValue));
}
// ***************************************************************************
TEST(TestStringsUtilities, TestParserChar_Test)
{
    string stringValue = "b";
    char correctValue = 'b';

    ASSERT_EQ(correctValue, Gedim::StringsUtilities::Parse<char>(stringValue));
}
// ***************************************************************************
TEST(TestStringsUtilities, TestParserString_Test)
{
    string stringValue = "b";
    string correctValue = "b";

    ASSERT_EQ(correctValue, Gedim::StringsUtilities::Parse<string>(stringValue));
}
// ***************************************************************************
TEST(TestStringsUtilities, TestParserDouble_Test)
{
    string stringValue = "1.5";
    double correctValue = 1.5;

    ASSERT_EQ(correctValue, Gedim::StringsUtilities::Parse<double>(stringValue));
}
// ***************************************************************************
TEST(TestStringsUtilities, TestParserVectorInt_Test)
{
    string stringValue = "[1 2 4]";
    vector<int> correctValue{1, 2, 4};

    ASSERT_EQ(correctValue, Gedim::StringsUtilities::Parse<vector<int>>(stringValue));
}
// ***************************************************************************
TEST(TestStringsUtilities, TestParserVectorDouble_Test)
{
    string stringValue = "[1.2 2.7 4.4]";
    vector<double> correctValue{1.2, 2.7, 4.4};

    ASSERT_EQ(correctValue, Gedim::StringsUtilities::Parse<vector<double>>(stringValue));
}
// ***************************************************************************
TEST(TestStringsUtilities, TestParserBool_Test)
{
    string stringValue = "false";
    bool correctValue = false;

    ASSERT_EQ(correctValue, Gedim::StringsUtilities::Parse<bool>(stringValue));
}
// ***************************************************************************
} // namespace GedimUnitTesting

#endif // __TEST_STRINGSUTILITIES_H
