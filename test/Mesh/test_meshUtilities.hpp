#ifndef __TEST_MESH_UTILITIES_H
#define __TEST_MESH_UTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

  TEST(TestMeshUtilities, TestMeshMatricesDAO)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    Eigen::MatrixXd polygon(3, 4);
    polygon.col(0)<< 0.0, 0.0, 0.0;
    polygon.col(1)<< 1.0, 0.0, 0.0;
    polygon.col(2)<< 1.0, 1.0, 0.0;
    polygon.col(3)<< 0.0, 1.0, 0.0;
    vector<unsigned int> vertexMarkers = { 1, 2, 3, 4 };
    vector<unsigned int> edgeMarkers = { 5, 6, 7, 8 };

    meshUtilities.Mesh2DFromPolygon(polygon,
                                    vertexMarkers,
                                    edgeMarkers,
                                    meshDao);

    EXPECT_EQ(meshDao.Dimension(), 2);
    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 4);
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 4);
    EXPECT_EQ(meshDao.Cell2DTotalNumber(), 1);
  }
}

#endif // __TEST_MESH_UTILITIES_H
