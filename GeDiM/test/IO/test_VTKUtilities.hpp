#ifndef __TEST_VTPUtilities_H
#define __TEST_VTPUtilities_H

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

namespace GedimUnitTesting
{
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test0D)
  {
    const unsigned int numGeometries = 4;

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      Eigen::Vector3d geometry(1.0 + g,
                               0.0 + g,
                               0.0 + g);

      vector<double> id(1, g + 1);
      vector<double> data(1, 10.8 + g);

      vtpUtilities.AddPoint(geometry,
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

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry0D.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  TEST(TestVTPUtilities, VTPUtilities_Test0Ds)
  {
    const unsigned int numGeometries = 3;

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      Eigen::MatrixXd geometry = (Eigen::MatrixXd(3, 4)<< 0.0, 1.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0, 1.0,
                                  2.0 + g, 2.0 + g, 2.0 + g, 2.0 + g).finished();

      vector<double> id(4);
      id[0] = 1 + g;
      id[1] = 2 + g;
      id[2] = 3 + g;
      id[3] = 4 + g;
      vector<double> data(4);
      data[0] = 10 + g;
      data[1] = 20 + g;
      data[2] = 30 + g;
      data[3] = 40 + g;

      vtpUtilities.AddPoints(geometry,
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

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry0Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test1D)
  {
    const unsigned int numGeometries = 4;

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      Eigen::MatrixXd geometry(3, 2);

      geometry.col(0)<< 1.0 + g, 0.0 + g, 0.0 + g;
      geometry.col(1)<< 0.0 + g, 1.0 + g, 1.0 + g;

      vector<double> id(1, g);
      vector<double> data(geometry.cols());

      for (unsigned int v = 0; v < geometry.cols(); v++)
        data[v] = 10.8 + g + v;

      vtpUtilities.AddSegment(geometry,
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

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry1D.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test1Ds)
  {
    const unsigned int numGeometries = 3;

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      const Eigen::MatrixXd vertices = (Eigen::MatrixXd(3, 4)<< 0.0, 1.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0, 1.0,
                                        2.0 + g, 2.0 + g, 2.0 + g, 2.0 + g).finished();
      const Eigen::MatrixXi edges = (Eigen::MatrixXi(2, 5)<< 0, 1, 2, 3, 0,
                                     1, 2, 3, 0, 2).finished();

      vector<double> id(5);
      vector<double> data(4);

      for (unsigned int e = 0; e < 4; e++)
        id[e] = 10.8 + g + e;

      for (unsigned int v = 0; v < 4; v++)
        data[v] = 10.8 + g + v;

      vtpUtilities.AddSegments(vertices,
                               edges,
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

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry1Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test2D)
  {
    const unsigned int numGeometries = 4;

    const double angle = M_PI / 4.0; // 45 degrees
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix.row(0) << cos(angle), 0.0, sin(angle);
    rotationMatrix.row(1) << 0.0, 1.0, 0.0;
    rotationMatrix.row(2) << -sin(angle), 0.0, cos(angle);

    const Eigen::Vector3d translation(0.0, 2.0, 0.0);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      Eigen::MatrixXd geometry(3, 4);
      geometry.col(0)<< 0.0 + g, 0.0 + g, 0.0;
      geometry.col(1)<< 1.0 + g, 0.0 + g, 0.0;
      geometry.col(2)<< 1.0 + g, 1.0 + g, 0.0;
      geometry.col(3)<< 0.0 + g, 1.0 + g, 0.0;

      geometry = geometryUtilities.RotatePointsFrom2DTo3D(geometry,
                                                          rotationMatrix,
                                                          translation);

      vector<double> id(1, g);
      vector<double> data(geometry.cols());

      for (unsigned int v = 0; v < geometry.cols(); v++)
        data[v] = 10.8 + g + v;

      vtpUtilities.AddPolygon(geometry,
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

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry2D.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test2Ds)
  {
    const unsigned int numGeometries = 4;

    const double angle = M_PI / 4.0; // 45 degrees
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix.row(0) << cos(angle), 0.0, sin(angle);
    rotationMatrix.row(1) << 0.0, 1.0, 0.0;
    rotationMatrix.row(2) << -sin(angle), 0.0, cos(angle);

    const Eigen::Vector3d translation(0.0, 2.0, 0.0);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < numGeometries; g++)
    {
      Eigen::MatrixXd vertices(3, 4);
      vertices.col(0)<< 0.0 + g, 0.0 + g, 0.0;
      vertices.col(1)<< 1.0 + g, 0.0 + g, 0.0;
      vertices.col(2)<< 1.0 + g, 1.0 + g, 0.0;
      vertices.col(3)<< 0.0 + g, 1.0 + g, 0.0;

      vertices = geometryUtilities.RotatePointsFrom2DTo3D(vertices,
                                                          rotationMatrix,
                                                          translation);

      vector<vector<unsigned int>> polygons =
      {
        { 0, 1, 2 },
        { 0, 2, 3 }
      };

      vector<double> id(polygons.size());
      vector<double> data(vertices.cols());

      for (unsigned int f = 0; f < id.size(); f++)
        id[f] = 1 + g + f;

      for (unsigned int v = 0; v < vertices.cols(); v++)
        data[v] = 10.8 + g + v;

      vtpUtilities.AddPolygons(vertices,
                               polygons,
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

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry2Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test3D)
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

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry3D.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_Test3Ds)
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

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Geometry3Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh2D_Cell0Ds)
  {
    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
    {
      const vector<double> id(1, g);
      const vector<double> marker(1, mesh.Cell0DMarker(g));

      vtpUtilities.AddPoint(mesh.Cell0DCoordinates(g),
                            {
                              {
                                "Id",
                                Gedim::VTPProperty::Formats::Cells,
                                static_cast<unsigned int>(id.size()),
                                id.data()
                              },
                              {
                                "Marker",
                                Gedim::VTPProperty::Formats::Points,
                                static_cast<unsigned int>(marker.size()),
                                marker.data()
                              }
                            });
    }

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Mesh2D_Cell0Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh2D_Cell1Ds)
  {
    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell1DTotalNumber(); g++)
    {
      const vector<double> id(1, g);
      const vector<double> marker(2, mesh.Cell1DMarker(g));

      vtpUtilities.AddSegment(mesh.Cell1DCoordinates(g),
                              {
                                {
                                  "Id",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(id.size()),
                                  id.data()
                                },
                                {
                                  "Marker",
                                  Gedim::VTPProperty::Formats::Points,
                                  static_cast<unsigned int>(marker.size()),
                                  marker.data()
                                }
                              });
    }

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Mesh2D_Cell1Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh2D_Cell2Ds)
  {
    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell2DTotalNumber(); g++)
    {
      const vector<double> id(1, g);
      const vector<double> marker(mesh.Cell2DNumberVertices(g), mesh.Cell2DMarker(g));

      vtpUtilities.AddPolygon(mesh.Cell2DVerticesCoordinates(g),
                              {
                                {
                                  "Id",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(id.size()),
                                  id.data()
                                },
                                {
                                  "Marker",
                                  Gedim::VTPProperty::Formats::Points,
                                  static_cast<unsigned int>(marker.size()),
                                  marker.data()
                                }
                              });
    }

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Mesh2D_Cell2Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_SinglePoints_Cell0Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
    {
      const vector<double> id(1, g);
      const vector<double> marker(1, mesh.Cell0DMarker(g));

      vtpUtilities.AddPoint(mesh.Cell0DCoordinates(g),
                            {
                              {
                                "Id",
                                Gedim::VTPProperty::Formats::Cells,
                                static_cast<unsigned int>(id.size()),
                                id.data()
                              },
                              {
                                "Marker",
                                Gedim::VTPProperty::Formats::Cells,
                                static_cast<unsigned int>(marker.size()),
                                marker.data()
                              }
                            });
    }

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Mesh3D_SinglePoints_Cell0Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_GlobalPoints_Cell0Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);

    // Export to VTK
    vector<double> id(mesh.Cell0DTotalNumber());
    vector<double> marker(mesh.Cell0DTotalNumber());

    for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
    {
      id[g] = g;
      marker[g] = mesh.Cell0DMarker(g);
    }

    Gedim::VTKUtilities vtpUtilities;
    vtpUtilities.AddPoints(mesh.Cell0DsCoordinates(),
                           {
                             {
                               "Id",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(id.size()),
                               id.data()
                             },
                             {
                               "Marker",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(marker.size()),
                               marker.data()
                             }
                           });

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Mesh3D_GlobalPoints_Cell0Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_SingleSegments_Cell1Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell1DTotalNumber(); g++)
    {
      const vector<double> id(1, g);
      const vector<double> marker(1, mesh.Cell1DMarker(g));

      vtpUtilities.AddSegment(mesh.Cell1DCoordinates(g),
                              {
                                {
                                  "Id",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(id.size()),
                                  id.data()
                                },
                                {
                                  "Marker",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(marker.size()),
                                  marker.data()
                                }
                              });
    }

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Mesh3D_SingleSegments_Cell1Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_GlobalSegments_Cell1Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    vector<double> id(mesh.Cell1DTotalNumber());
    vector<double> marker(mesh.Cell1DTotalNumber());

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell1DTotalNumber(); g++)
    {
      id[g] = g;
      marker[g] = mesh.Cell1DMarker(g);
    }

    vtpUtilities.AddSegments(mesh.Cell0DsCoordinates(),
                             mesh.Cell1DsExtremes(),
                             {
                               {
                                 "Id",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(id.size()),
                                 id.data()
                               },
                               {
                                 "Marker",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(marker.size()),
                                 marker.data()
                               }
                             });
    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Mesh3D_GlobalSegments_Cell1Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_SinglePolygons_Cell2Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell2DTotalNumber(); g++)
    {
      const vector<double> id(1, g);
      const vector<double> marker(1, mesh.Cell2DMarker(g));

      vtpUtilities.AddPolygon(mesh.Cell2DVerticesCoordinates(g),
                              {
                                {
                                  "Id",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(id.size()),
                                  id.data()
                                },
                                {
                                  "Marker",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(marker.size()),
                                  marker.data()
                                }
                              });
    }

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Mesh3D_SinglePolygons_Cell2Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_GlobalPolygons_Cell2Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    vector<double> id(mesh.Cell2DTotalNumber());
    vector<double> marker(mesh.Cell2DTotalNumber());
    for (unsigned int g = 0; g < mesh.Cell2DTotalNumber(); g++)
    {
      id[g] = g;
      marker[g] = mesh.Cell2DMarker(g);
    }

    vtpUtilities.AddPolygons(mesh.Cell0DsCoordinates(),
                             mesh.Cell2DsVertices(),
                             {
                               {
                                 "Id",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(id.size()),
                                 id.data()
                               },
                               {
                                 "Marker",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(marker.size()),
                                 marker.data()
                               }
                             });

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Mesh3D_GlobalPolygons_Cell2Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_SinglePolyhedrons_Cell3Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;
    Gedim::MeshUtilities meshUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell3DTotalNumber(); g++)
    {
      const Gedim::MeshUtilities::VTPPolyhedron polyhedron = meshUtilities.MeshCell3DToVTPPolyhedron(mesh,
                                                                                                     g);

      const vector<double> id(1, g);
      const vector<double> marker(polyhedron.Vertices.cols(),
                                  mesh.Cell3DMarker(g));

      vtpUtilities.AddPolyhedron(polyhedron.Vertices,
                                 polyhedron.PolyhedronFaces,
                                 {
                                   {
                                     "Id",
                                     Gedim::VTPProperty::Formats::Cells,
                                     static_cast<unsigned int>(id.size()),
                                     id.data()
                                   },
                                   {
                                     "Marker",
                                     Gedim::VTPProperty::Formats::Points,
                                     static_cast<unsigned int>(marker.size()),
                                     marker.data()
                                   }
                                 });
    }

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Mesh3D_SinglePolyhedrons_Cell3Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_GlobalPolyhedrons_Cell3Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;
    Gedim::MeshUtilities meshUtilities;

    // Export to VTK
    vector<double> id(mesh.Cell3DTotalNumber());
    vector<double> marker(mesh.Cell3DTotalNumber());
    for (unsigned int g = 0; g < mesh.Cell3DTotalNumber(); g++)
    {
      id[g] = g;
      marker[g] = mesh.Cell3DMarker(g);
    }

    vtpUtilities.AddPolyhedrons(mesh.Cell0DsCoordinates(),
                                mesh.Cell3DsFacesVertices(),
                                {
                                  {
                                    "Id",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(id.size()),
                                    id.data()
                                  },
                                  {
                                    "Marker",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(marker.size()),
                                    marker.data()
                                  }
                                });

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    vtpUtilities.Export(exportFolder + "/Mesh3D_GlobalPolyhedrons_Cell3Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_ExportMesh)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;
    Gedim::MeshUtilities meshUtilities;

    std::string exportFolder = "./Export/TestVTPUtilities";
    Gedim::Output::CreateFolder(exportFolder);

    meshUtilities.ExportMeshToVTU(mesh,
                                  exportFolder,
                                  "ExportMesh");
  }
  // ***************************************************************************
}

#endif // __TEST_VTPUtilities_H
