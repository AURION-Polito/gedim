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
  TEST(TestVTPUtilities, VTPUtilities_Test2D)
  {
    const unsigned int numGeometries = 4;


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

      vector<double> id(geometry.Faces.size(), g);
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
  TEST(TestVTPUtilities, VTPUtilities_TestMesh2D_Cell0Ds)
  {
    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
    {
      vector<double> id(1, mesh.Cell0DId(g));
      vector<double> marker(1, mesh.Cell0DMarker(g));

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
      vector<double> id(1, mesh.Cell1DId(g));
      vector<double> marker(2, mesh.Cell1DMarker(g));

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
      vector<double> id(1, mesh.Cell2DId(g));
      vector<double> marker(mesh.Cell2DNumberVertices(g), mesh.Cell2DMarker(g));

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
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_Cell0Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
    {
      vector<double> id(1, mesh.Cell0DId(g));
      vector<double> marker(1, mesh.Cell0DMarker(g));

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

    vtpUtilities.Export(exportFolder + "/Mesh3D_Cell0Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_Cell1Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell1DTotalNumber(); g++)
    {
      vector<double> id(1, mesh.Cell1DId(g));
      vector<double> marker(1, mesh.Cell1DMarker(g));

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

    vtpUtilities.Export(exportFolder + "/Mesh3D_Cell1Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_Cell2Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell2DTotalNumber(); g++)
    {
      vector<double> id(1, mesh.Cell2DId(g));
      vector<double> marker(1, mesh.Cell2DMarker(g));

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

    vtpUtilities.Export(exportFolder + "/Mesh3D_Cell2Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
  TEST(TestVTPUtilities, VTPUtilities_TestMesh3D_Cell3Ds)
  {
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);
    Gedim::VTKUtilities vtpUtilities;
    Gedim::MeshUtilities meshUtilities;

    // Export to VTK
    for (unsigned int g = 0; g < mesh.Cell3DTotalNumber(); g++)
    {
      const Gedim::GeometryUtilities::Polyhedron polyhedron = meshUtilities.MeshCell3DToPolyhedron(mesh,
                                                                                                   g);

      vector<double> id(polyhedron.Vertices.cols(),
                        mesh.Cell3DId(g));
      vector<double> marker(polyhedron.Vertices.cols(),
                            mesh.Cell3DMarker(g));

      vtpUtilities.AddPolyhedron(polyhedron.Vertices,
                                 polyhedron.Edges,
                                 polyhedron.Faces,
                                 {
                                   {
                                     "Id",
                                     Gedim::VTPProperty::Formats::Points,
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

    vtpUtilities.Export(exportFolder + "/Mesh3D_Cell3Ds.vtu",
                        Gedim::VTKUtilities::Ascii);
  }
  // ***************************************************************************
}

#endif // __TEST_VTPUtilities_H
