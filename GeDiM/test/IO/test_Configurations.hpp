#ifndef __TEST_CONFIGURATIONS_H
#define __TEST_CONFIGURATIONS_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "Configurations.hpp"
#include "IOUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
  // ***************************************************************************
  TEST(TestConfigurations, AddGenericProperty_Test)
  {
    string propertyId = "Test";
    string propertyIdTwo = "Test1";

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_EQ(0, Gedim::Configurations::NumberProperties());

    ASSERT_NO_THROW(Gedim::Configurations::AddProperty<double>(propertyId));
    ASSERT_EQ(1, Gedim::Configurations::NumberProperties());
    ASSERT_TRUE(Gedim::Configurations::ExistsProperty(propertyId));
    ASSERT_TRUE(!Gedim::Configurations::ExistsProperty(propertyIdTwo));
    ASSERT_TRUE(Gedim::Configurations::CheckPropertyType<double>(propertyId));
    ASSERT_TRUE(!Gedim::Configurations::CheckPropertyType<vector<double>>(propertyId));

    ASSERT_NO_THROW(Gedim::Configurations::AddProperty<double>(propertyIdTwo));
    ASSERT_EQ(2, Gedim::Configurations::NumberProperties());

    ASSERT_NO_THROW(Gedim::Configurations::RemoveProperty(propertyIdTwo));
    ASSERT_EQ(1, Gedim::Configurations::NumberProperties());

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_EQ(0, Gedim::Configurations::NumberProperties());
  }
  // ***************************************************************************
  TEST(TestConfigurations, AddIntProperty_Test)
  {
    string propertyId = "Test";
    int value = 2;

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty<int>(propertyId));
    ASSERT_NO_THROW(Gedim::Configurations::SetPropertyValue(propertyId, value));
    ASSERT_EQ(value, Gedim::Configurations::GetPropertyValue<int>(propertyId));
  }
  // ***************************************************************************
  TEST(TestConfigurations, AddUnsignedIntProperty_Test)
  {
    string propertyId = "Test";
    unsigned int value = 2;

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty<unsigned int>(propertyId));
    ASSERT_NO_THROW(Gedim::Configurations::SetPropertyValue(propertyId, value));
    ASSERT_EQ(value, Gedim::Configurations::GetPropertyValue<unsigned int>(propertyId));
  }
  // ***************************************************************************
  TEST(TestConfigurations, AddCharProperty_Test)
  {
    string propertyId = "Test";
    char value = 'b';

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty<char>(propertyId));
    ASSERT_NO_THROW(Gedim::Configurations::SetPropertyValue(propertyId, value));
    ASSERT_EQ(value, Gedim::Configurations::GetPropertyValue<char>(propertyId));
  }
  // ***************************************************************************
  TEST(TestConfigurations, AddStringProperty_Test)
  {
    string propertyId = "Test";
    string value = "b";

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty<string>(propertyId));
    ASSERT_NO_THROW(Gedim::Configurations::SetPropertyValue(propertyId, value));
    ASSERT_EQ(value, Gedim::Configurations::GetPropertyValue<string>(propertyId));
  }
  // ***************************************************************************
  TEST(TestConfigurations, AddConstCharProperty_Test)
  {
    string propertyId = "Test";
    const char* value = "b";

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty(propertyId, value));
    ASSERT_NO_THROW(Gedim::Configurations::SetPropertyValue(propertyId, value));
    ASSERT_EQ(value, Gedim::Configurations::GetPropertyValue<string>(propertyId));
  }
  // ***************************************************************************
  TEST(TestConfigurations, AddDoubleProperty_Test)
  {
    string propertyId = "Test";
    double value = 3.4;

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty<double>(propertyId));
    ASSERT_NO_THROW(Gedim::Configurations::SetPropertyValue(propertyId, value));
    ASSERT_EQ(value, Gedim::Configurations::GetPropertyValue<double>(propertyId));
  }
  // ***************************************************************************
  TEST(TestConfigurations, AddVectorIntProperty_Test)
  {
    string propertyId = "Test";
    vector<int> value{ 11, 22, 31 };

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty<vector<int>>(propertyId));
    ASSERT_NO_THROW(Gedim::Configurations::SetPropertyValue(propertyId, value));
    ASSERT_EQ(value, Gedim::Configurations::GetPropertyValue<vector<int>>(propertyId));
  }
  // ***************************************************************************
  TEST(TestConfigurations, AddVectorUnsignedIntProperty_Test)
  {
    string propertyId = "Test";
    vector<unsigned int> value{ 11, 22, 31 };

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty<vector<unsigned int>>(propertyId));
    ASSERT_NO_THROW(Gedim::Configurations::SetPropertyValue(propertyId, value));
    ASSERT_EQ(value, Gedim::Configurations::GetPropertyValue<vector<unsigned int>>(propertyId));
  }
  // ***************************************************************************
  TEST(TestConfigurations, AddVectorDoubleProperty_Test)
  {
    string propertyId = "Test";
    vector<double> value{ 11.4, 22.7, 31.9 };

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty<vector<double>>(propertyId));
    ASSERT_NO_THROW(Gedim::Configurations::SetPropertyValue(propertyId, value));
    ASSERT_EQ(value, Gedim::Configurations::GetPropertyValue<vector<double>>(propertyId));
  }
  // ***************************************************************************
  TEST(TestConfigurations, AddBoolProperty_Test)
  {
    string propertyId = "Test";
    bool value = false;

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty<bool>(propertyId));
    ASSERT_NO_THROW(Gedim::Configurations::SetPropertyValue(propertyId, value));
    ASSERT_EQ(value, Gedim::Configurations::GetPropertyValue<bool>(propertyId));
  }
  // ***************************************************************************
  TEST(TestConfigurations, InitializeArgcArgv_Test)
  {
    vector<char*> argv = {
      (char*)"FakeTest",
      (char*)"TestString=202",
      (char*)"FakeTest",
      (char*)"TestInt:int=202",
      (char*)"TestUnsignedInt:uint=103",
      (char*)"TestChar:char=b",
      (char*)"TestDouble:double=202.79",
      (char*)"TestBool:bool=true",
      (char*)"TestVectorInt:vector_int=[1,2,3]",
      (char*)"TestVectorUInt:vector_uint=[7,4,5]",
      (char*)"TestVectorDouble:vector_double=[1.5,1.6]",
    };

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::Initialize(argv.size(), argv.data()));
    ASSERT_EQ(9, Gedim::Configurations::NumberProperties());
    ASSERT_TRUE(Gedim::Configurations::ExistsProperty("TestInt"));
    ASSERT_TRUE(Gedim::Configurations::ExistsProperty("TestUnsignedInt"));
    ASSERT_TRUE(Gedim::Configurations::ExistsProperty("TestChar"));
    ASSERT_TRUE(Gedim::Configurations::ExistsProperty("TestDouble"));
    ASSERT_TRUE(Gedim::Configurations::ExistsProperty("TestString"));
    ASSERT_TRUE(Gedim::Configurations::ExistsProperty("TestBool"));
    ASSERT_TRUE(Gedim::Configurations::ExistsProperty("TestVectorInt"));
    ASSERT_TRUE(Gedim::Configurations::ExistsProperty("TestVectorUInt"));
    ASSERT_TRUE(Gedim::Configurations::ExistsProperty("TestVectorDouble"));

    ASSERT_EQ(202, Gedim::Configurations::GetPropertyValue<int>("TestInt"));
    ASSERT_EQ(103, Gedim::Configurations::GetPropertyValue<unsigned int>("TestUnsignedInt"));
    ASSERT_EQ('b', Gedim::Configurations::GetPropertyValue<char>("TestChar"));
    ASSERT_EQ(202.79, Gedim::Configurations::GetPropertyValue<double>("TestDouble"));
    ASSERT_EQ("202", Gedim::Configurations::GetPropertyValue<string>("TestString"));
    ASSERT_EQ(true, Gedim::Configurations::GetPropertyValue<bool>("TestBool"));
    ASSERT_EQ(vector<int>({1, 2, 3}), Gedim::Configurations::GetPropertyValue<vector<int>>("TestVectorInt"));
    ASSERT_EQ(vector<unsigned int>({7, 4, 5}), Gedim::Configurations::GetPropertyValue<vector<unsigned int>>("TestVectorUInt"));
    ASSERT_EQ(vector<double>({1.5, 1.6}), Gedim::Configurations::GetPropertyValue<vector<double>>("TestVectorDouble"));
  }
  // ***************************************************************************
  TEST(TestConfigurations, ExportToCsv_Test)
  {
    ASSERT_NO_THROW(Gedim::Configurations::Reset());

    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestInt", 10, "int test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestChar", 'b', "char test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestDouble", 12.6, "double test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestString", "pippo", "string test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestVectorInt", vector<int>{1, 7, 8}, "vector<int> test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestVectorDouble", vector<double>{1.6, 8.7, 10.8}, "vector<double> test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestBool", true, "bool test"));

    string exportFolder = "./Export/TestConfigurations/ExportToCsv_Test";
    Gedim::Output::CreateFolder(exportFolder);
    ASSERT_NO_THROW(Gedim::Configurations::ExportToCsv(exportFolder + "/Parameters.csv", false, ';'));

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::InitializeFromCsv(exportFolder + "/Parameters.csv", true, ';'));
    ASSERT_EQ(7, Gedim::Configurations::NumberProperties());

    ASSERT_EQ(10, Gedim::Configurations::GetPropertyValue<int>("TestInt"));
    ASSERT_EQ('b', Gedim::Configurations::GetPropertyValue<char>("TestChar"));
    ASSERT_EQ(12.6, Gedim::Configurations::GetPropertyValue<double>("TestDouble"));
    ASSERT_EQ("pippo", Gedim::Configurations::GetPropertyValue<string>("TestString"));
    ASSERT_EQ(true, Gedim::Configurations::GetPropertyValue<bool>("TestBool"));
    ASSERT_EQ(vector<int>({1, 7, 8}), Gedim::Configurations::GetPropertyValue<vector<int>>("TestVectorInt"));
    ASSERT_EQ(vector<double>({1.6, 8.7, 10.8}), Gedim::Configurations::GetPropertyValue<vector<double>>("TestVectorDouble"));
  }
  // ***************************************************************************
  TEST(TestConfigurations, ExportToIni_Test)
  {
    ASSERT_NO_THROW(Gedim::Configurations::Reset());

    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestInt", 10, "int test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestChar", 'b', "char test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestDouble", 12.6, "double test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestString", "pippo", "string test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestVectorInt", vector<int>{1, 7, 8}, "vector<int> test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestVectorDouble", vector<double>{1.6, 8.7, 10.8}, "vector<double> test"));
    ASSERT_NO_THROW(Gedim::Configurations::AddProperty("TestBool", true, "bool test"));

    string exportFolder = "./Export/TestConfigurations/ExportToIni_Test";
    Gedim::Output::CreateFolder(exportFolder);
    ASSERT_NO_THROW(Gedim::Configurations::ExportToIni(exportFolder + "/Parameters.ini", false, "Test"));

    ASSERT_NO_THROW(Gedim::Configurations::Reset());
    ASSERT_NO_THROW(Gedim::Configurations::InitializeFromIni(exportFolder + "/Parameters.ini"));
    ASSERT_EQ(7, Gedim::Configurations::NumberProperties());

    ASSERT_EQ(10, Gedim::Configurations::GetPropertyValue<int>("TestInt"));
    ASSERT_EQ('b', Gedim::Configurations::GetPropertyValue<char>("TestChar"));
    ASSERT_EQ(12.6, Gedim::Configurations::GetPropertyValue<double>("TestDouble"));
    ASSERT_EQ("pippo", Gedim::Configurations::GetPropertyValue<string>("TestString"));
    ASSERT_EQ(true, Gedim::Configurations::GetPropertyValue<bool>("TestBool"));
    ASSERT_EQ(vector<int>({1, 7, 8}), Gedim::Configurations::GetPropertyValue<vector<int>>("TestVectorInt"));
    ASSERT_EQ(vector<double>({1.6, 8.7, 10.8}), Gedim::Configurations::GetPropertyValue<vector<double>>("TestVectorDouble"));
  }
  // ***************************************************************************
}

#endif // __TEST_CONFIGURATIONS_H
