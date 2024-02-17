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
    const unsigned int numGeometries = 1;


    double angle = M_PI / 4.0; // 45 degrees
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix.row(0) << cos(angle), 0.0, sin(angle);
    rotationMatrix.row(1) << 0.0, 1.0, 0.0;
    rotationMatrix.row(2) << -sin(angle), 0.0, cos(angle);

    Eigen::Vector3d translation(0.0, 2.0, 0.0);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      Gedim::GeometryUtilities::Polyhedron geometry = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0 + g,
                                                                                                             0.0 + g,
                                                                                                             0.0 + g),
                                                                                             1.0);

      geometry.Vertices = geometryUtilities.RotatePoints(geometry.Vertices,
                                                         rotationMatrix,
                                                         translation);

      vector<double> id(1, g);
      vector<double> data(geometry.Vertices.cols());

      for (unsigned int v = 0; v < geometry.Vertices.cols(); v++)
        data[v] = 10.8 + g + v;

      vtpUtilities.AddPolyhedron(geometry.Vertices,
                                 geometry.Edges,
                                 geometry.Faces,
                                 {
                                   {
                                     "Id",
                                     Gedim::VTPProperty::Formats::Cells,
                                     static_cast<unsigned int>(id.size()),
                                     id.data()
                                   },
                                   {
                                     "Data",
                                     Gedim::VTPProperty::Formats::Points,
                                     static_cast<unsigned int>(data.size()),
                                     data.data()
                                   }
                                 });
    }

    std::string exportFolder = "./Export/TestUCDUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry3D.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestUCDUtilities, UCDUtilities_Test3Ds)
  {
    const unsigned int numGeometries = 1;


    double angle = M_PI / 4.0; // 45 degrees
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix.row(0) << cos(angle), 0.0, sin(angle);
    rotationMatrix.row(1) << 0.0, 1.0, 0.0;
    rotationMatrix.row(2) << -sin(angle), 0.0, cos(angle);

    Eigen::Vector3d translation(0.0, 2.0, 0.0);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      Gedim::GeometryUtilities::Polyhedron geometry = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0 + g,
                                                                                                             0.0 + g,
                                                                                                             0.0 + g),
                                                                                             1.0);

      geometry.Vertices = geometryUtilities.RotatePoints(geometry.Vertices,
                                                         rotationMatrix,
                                                         translation);
      const vector<vector<vector<unsigned int>>> polyhedronsFaces =
      {
        {
          { 0, 1, 2 },
          { 0, 1, 5 },
          { 0, 2, 5 },
          { 1, 2, 5 }
        },
        {
          { 0, 2, 3 },
          { 0, 2, 4 },
          { 0, 3, 4 },
          { 2, 3, 4 }
        }
      };


      vector<double> id(polyhedronsFaces.size());
      id[0] = 10 + g;
      id[1] = 20 + g;
      vector<double> data(geometry.Vertices.cols());

      for (unsigned int v = 0; v < geometry.Vertices.cols(); v++)
        data[v] = 10.8 + g + v;

      vtpUtilities.AddPolyhedrons(geometry.Vertices,
                                  polyhedronsFaces,
                                  {
                                    {
                                      "Id",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(id.size()),
                                      id.data()
                                    },
                                    {
                                      "Data",
                                      Gedim::VTPProperty::Formats::Points,
                                      static_cast<unsigned int>(data.size()),
                                      data.data()
                                    }
                                  });
    }

    std::string exportFolder = "./Export/TestUCDUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry3Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
}

#endif // __TEST_UCDUtilities_H
