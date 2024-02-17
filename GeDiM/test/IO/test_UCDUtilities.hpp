#ifndef __TEST_UCDUtilities_H
#define __TEST_UCDUtilities_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "MeshMatrices_3D_329Cells_Mock.hpp"

#include "GeometryUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "IOUtilities.hpp"
#include "VTKUtilities.hpp"
#include "UCDUtilities.hpp"

namespace GedimUnitTesting
{
  // ***************************************************************************
  TEST(TestUCDUtilities, UCDUtilities_Test0Ds)
  {
    std::string exportFolder = "./Export/TestUCDUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    const unsigned int numGeometries = 4;

    Gedim::UCDUtilities exporter;

    // Export to UCD
    Eigen::MatrixXd points(3, numGeometries);
    vector<double> id(numGeometries);
    vector<double> data(2 * numGeometries);
    Eigen::VectorXi material(numGeometries);

    for (unsigned int g = 0; g < numGeometries; g++)
    {
      points.col(g)<< 1.0 + g,  0.0 + g, 0.0 + g;

      id[g] = g + 1;
      data[2 * g] = g + 1;
      data[2 * g + 1] = g + 1;
      material[g] = g + 1;
    }

    exporter.ExportPoints(exportFolder + "/Geometry0Ds.inp",
                          points,
                          {
                            {
                              "Id",
                              "kg",
                              static_cast<unsigned int>(id.size()),
                              1,
                              id.data()
                            },
                            {
                              "Data",
                              "m",
                              static_cast<unsigned int>(data.size()),
                              2,
                              data.data()
                            }
                          },
                          material);

  }
  // ***************************************************************************
  TEST(TestUCDUtilities, UCDUtilities_Test1Ds)
  {
    std::string exportFolder = "./Export/TestUCDUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    const unsigned int numGeometries = 5;

    Gedim::UCDUtilities exporter;

    // Export to UCD
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 4)<< 0.0, 1.0, 1.0, 0.0,
                                    0.0, 0.0, 1.0, 1.0,
                                    2.0, 2.0, 2.0, 2.0).finished();
    const Eigen::MatrixXi edges = (Eigen::MatrixXi(2, 5)<< 0, 1, 2, 3, 0,
                                   1, 2, 3, 0, 2).finished();

    vector<double> id_points = { 1, 2, 3, 4 };
    vector<double> id(numGeometries);
    vector<double> data(2 * numGeometries);
    Eigen::VectorXi material(numGeometries);

    for (unsigned int g = 0; g < numGeometries; g++)
    {
      id[g] = g + 1;
      data[2 * g] = g + 1;
      data[2 * g + 1] = g + 1;
      material[g] = g + 1;
    }

    exporter.ExportSegments(exportFolder + "/Geometry1Ds.inp",
                            points,
                            edges,
                            {
                              {
                                "id_points",
                                "kg",
                                static_cast<unsigned int>(id_points.size()),
                                1,
                                id_points.data()
                              }
                            },
                            {
                              {
                                "Id",
                                "kg",
                                static_cast<unsigned int>(id.size()),
                                1,
                                id.data()
                              },
                              {
                                "Data",
                                "m",
                                static_cast<unsigned int>(data.size()),
                                2,
                                data.data()
                              }
                            },
                            material);
  }
  // ***************************************************************************
  TEST(TestUCDUtilities, UCDUtilities_Test2Ds)
  {
    std::string exportFolder = "./Export/TestUCDUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    const unsigned int numGeometries = 3;

    Gedim::UCDUtilities exporter;

    // Export to UCD
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 8)<< 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
                                    0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
                                    2.0, 2.0, 2.0, 2.0, 4.0, 4.0, 4.0, 4.0).finished();
    const vector<vector<unsigned int>> polygons =
    {
      { 0, 1, 2 },
      { 0, 2, 3 },
      { 4, 5, 6, 7 }
    };

    vector<double> id_points = { 1, 2, 3, 4, 5, 6, 7, 8 };
    vector<double> id(numGeometries);
    vector<double> data(2 * numGeometries);
    Eigen::VectorXi material(numGeometries);

    for (unsigned int g = 0; g < numGeometries; g++)
    {
      id[g] = g + 1;
      data[2 * g] = g + 1;
      data[2 * g + 1] = g + 1;
      material[g] = g + 1;
    }

    exporter.ExportPolygons(exportFolder + "/Geometry2Ds.inp",
                            points,
                            polygons,
                            {
                              {
                                "id_points",
                                "kg",
                                static_cast<unsigned int>(id_points.size()),
                                1,
                                id_points.data()
                              }
                            },
                            {
                              {
                                "Id",
                                "kg",
                                static_cast<unsigned int>(id.size()),
                                1,
                                id.data()
                              },
                              {
                                "Data",
                                "m",
                                static_cast<unsigned int>(data.size()),
                                2,
                                data.data()
                              }
                            },
                            material);
  }
  // ***************************************************************************
  TEST(TestUCDUtilities, UCDUtilities_Test3D)
  {
    std::string exportFolder = "./Export/TestUCDUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    const unsigned int numGeometries = 2;

    Gedim::UCDUtilities exporter;

    // Export to UCD
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 9)<<
                                    0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
                                    0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0).finished();
    const vector<vector<unsigned int>> polyhedra =
    {
      { 0, 2, 1, 5 },
      { 5, 4, 7, 8 }
    };

    vector<double> id_points = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    vector<double> id(numGeometries);
    vector<double> data(2 * numGeometries);
    Eigen::VectorXi material(numGeometries);

    for (unsigned int g = 0; g < numGeometries; g++)
    {
      id[g] = g + 1;
      data[2 * g] = g + 1;
      data[2 * g + 1] = g + 1;
      material[g] = g + 1;
    }

    exporter.ExportPolyhedra(exportFolder + "/Geometry3Ds.inp",
                             points,
                             polyhedra,
                             {
                               {
                                 "id_points",
                                 "kg",
                                 static_cast<unsigned int>(id_points.size()),
                                 1,
                                 id_points.data()
                               }
                             },
                             {
                               {
                                 "Id",
                                 "kg",
                                 static_cast<unsigned int>(id.size()),
                                 1,
                                 id.data()
                               },
                               {
                                 "Data",
                                 "m",
                                 static_cast<unsigned int>(data.size()),
                                 2,
                                 data.data()
                               }
                             },
                             material);
  }
  // ***************************************************************************
  TEST(TestUCDUtilities, UCDUtilities_TestMesh3D)
  {
    std::string exportFolder = "./Export/TestUCDUtilities/TestMesh3D";
    Gedim::Output::CreateFolder(exportFolder);

    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);

    Gedim::MeshUtilities meshUtilities;
    meshUtilities.ExportMeshToUCD(mesh,
                                  exportFolder,
                                  "Mesh3D",
                                  false);
  }
  // ***************************************************************************

}

#endif // __TEST_UCDUtilities_H
