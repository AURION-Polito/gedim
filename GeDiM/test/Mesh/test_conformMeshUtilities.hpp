#ifndef __TEST_ConformMeshUtilities_H
#define __TEST_ConformMeshUtilities_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MeshMatrices.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "MeshMatrices_2D_4Cells_Mock.hpp"
#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "MeshMatricesDAO.hpp"

#include "MeshDAOExporterToCsv.hpp"
#include "MeshFromCsvUtilities.hpp"
#include "MeshDAOImporterFromCsv.hpp"
#include "FileTextReader.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

  TEST(TestConformMeshUtilities, TestConformMeshSingleSegment)
  {
    
  }
}

#endif // __TEST_ConformMeshUtilities_H
