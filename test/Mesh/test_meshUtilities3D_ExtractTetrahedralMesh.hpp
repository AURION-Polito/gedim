// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __TEST_meshUtilities3D_ExtractTetrahedralMesh_H
#define __TEST_meshUtilities3D_ExtractTetrahedralMesh_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshUtilities.hpp"
#include "TetgenInterface.hpp"
#include "VTKUtilities.hpp"
#include "VtkMeshExporter.hpp"
#include "VoroInterface.hpp"

namespace GedimUnitTesting
{
  TEST(TestMeshUtilities, TestMeshUtilities_2D_ExtractTriangularMesh)
  {
    std::string exportFolder = "./Export/TestMeshUtilities/TestMeshUtilities_2D_ExtractTriangularMesh";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-6;
    geometry_utilities_config.Tolerance2D = 1e-12;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
      const auto square = geometry_utilities.CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);

      Gedim::External::MeshUtilities::CreatePolygonalMesh(geometry_utilities,
                                                          square,
                                                          10,
                                                          1000,
                                                          mesh,
                                                          10);
    }

    {
      const std::string mesh_folder = exportFolder + "/mesh";
      Gedim::Output::CreateFolder(mesh_folder);
      Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_folder, "mesh");
    }

    const auto cell2Ds_triangles = Gedim::MeshUtilities::ComputeMeshCell2DsTriangles(geometry_utilities,
                                                                                     Gedim::MeshUtilities::MeshCell2DsTrianglesData::TrianglesTypes::ByBarycenter,
                                                                                     mesh);

    const auto extract_tetra = Gedim::MeshUtilities::ExtractTriangularMesh(cell2Ds_triangles,
                                                                           mesh);


    {
      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < extract_tetra.vertices.size(); ++c)
        {
          const auto& vertex = extract_tetra.vertices.at(c);
          const double original_cell0D = vertex.original_cell0D_index == Gedim::MeshUtilities::None ?
                                           -1.0 :
                                           vertex.original_cell0D_index;

          exporter.AddPoint(
                vertex.coordinate,
                {
                  {"original_cell0D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &original_cell0D}
                });
        }

        exporter.Export(exportFolder + "/extract_vertices.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int e = 0; e < extract_tetra.edges.size(); ++e)
        {
          const auto& edge = extract_tetra.edges.at(e);

          const double original_cell1D = edge.original_cell1D_index == Gedim::MeshUtilities::None ?
                                           -1.0 :
                                           edge.original_cell1D_index;

          exporter.AddSegment(
                extract_tetra.vertices.at(edge.extremes.at(0)).coordinate,
                extract_tetra.vertices.at(edge.extremes.at(1)).coordinate,
                {
                  {"original_cell1D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &original_cell1D},
                });
        }

        exporter.Export(exportFolder + "/extract_edges.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int f = 0; f < extract_tetra.triangles.size(); ++f)
        {
          const auto& triangle = extract_tetra.triangles.at(f);

          const double original_cell2D = triangle.original_cell2D_index == Gedim::MeshUtilities::None ?
                                           -1.0 :
                                           triangle.original_cell2D_index;

          Eigen::Matrix3d vertices;
          for (unsigned int t_v = 0; t_v < 3; ++t_v)
            vertices.col(t_v) = extract_tetra.vertices.at(triangle.vertices.at(t_v)).coordinate;

          exporter.AddPolygon(vertices,
                              {
                                {"original_cell2D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &original_cell2D},
                              });
        }

        exporter.Export(exportFolder + "/extract_triangles.vtu");
      }
    }

    Gedim::MeshMatrices triangles_mesh_data;
    Gedim::MeshMatricesDAO triangles_mesh(triangles_mesh_data);
    Gedim::MeshUtilities::TetrahedralMeshDataToMeshDAO(extract_tetra,
                                                       mesh,
                                                       triangles_mesh);

    {
      const std::string mesh_folder = exportFolder + "/triangles_mesh";
      Gedim::Output::CreateFolder(mesh_folder);
      Gedim::External::MeshUtilities::ExportMeshToVTU(triangles_mesh, mesh_folder, "mesh");
    }
  }

  TEST(TestMeshUtilities, TestMeshUtilities_3D_ExtractTetrahedralMesh)
  {
    std::string exportFolder = "./Export/TestMeshUtilities/TestMeshUtilities_3D_ExtractTetrahedralMesh";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-6;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
      const auto cube = geometry_utilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);

      Gedim::External::MeshUtilities::CreateTetrahedralMesh(cube.Vertices,
                                                           cube.Edges,
                                                           cube.Faces,
                                                           0.1,
                                                           mesh);

//      Gedim::External::MeshUtilities::CreatePolyhedralMesh(geometry_utilities,
//                                                           cube.Vertices,
//                                                           cube.Edges,
//                                                           cube.Faces,
//                                                           10,
//                                                           1000,
//                                                           mesh,
//                                                           10);
    }

    {
      const std::string mesh_folder = exportFolder + "/mesh";
      Gedim::Output::CreateFolder(mesh_folder);
      Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_folder, "mesh");
    }

    const auto cell2Ds_triangles = Gedim::MeshUtilities::ComputeMeshCell2DsTriangles(geometry_utilities,
                                                                                     Gedim::MeshUtilities::MeshCell2DsTrianglesData::TrianglesTypes::ByBarycenter,
                                                                                     mesh);

    const auto cell3Ds_tetra = Gedim::MeshUtilities::ComputeMeshCell3DsTetrahedra(geometry_utilities,
                                                                                  cell2Ds_triangles,
                                                                                  Gedim::MeshUtilities::MeshCell3DsTetrahedraData::TetrahedraTypes::ByBarycenter,
                                                                                  mesh);

    const auto extract_tetra = Gedim::MeshUtilities::ExtractTetrahedralMesh(cell2Ds_triangles,
                                                                            cell3Ds_tetra,
                                                                            mesh);


    {
      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < extract_tetra.vertices.size(); ++c)
        {
          const auto& vertex = extract_tetra.vertices.at(c);
          const double original_cell0D = vertex.original_cell0D_index == Gedim::MeshUtilities::None ?
                                           -1.0 :
                                           vertex.original_cell0D_index;

          exporter.AddPoint(
                vertex.coordinate,
                {
                  {"original_cell0D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &original_cell0D}
                });
        }

        exporter.Export(exportFolder + "/extract_vertices.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int e = 0; e < extract_tetra.edges.size(); ++e)
        {
          const auto& edge = extract_tetra.edges.at(e);

          const double original_cell1D = edge.original_cell1D_index == Gedim::MeshUtilities::None ?
                                           -1.0 :
                                           edge.original_cell1D_index;

          exporter.AddSegment(
                extract_tetra.vertices.at(edge.extremes.at(0)).coordinate,
                extract_tetra.vertices.at(edge.extremes.at(1)).coordinate,
                {
                  {"original_cell1D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &original_cell1D},
                });
        }

        exporter.Export(exportFolder + "/extract_edges.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int f = 0; f < extract_tetra.triangles.size(); ++f)
        {
          const auto& triangle = extract_tetra.triangles.at(f);

          const double original_cell2D = triangle.original_cell2D_index == Gedim::MeshUtilities::None ?
                                           -1.0 :
                                           triangle.original_cell2D_index;

          Eigen::Matrix3d vertices;
          for (unsigned int t_v = 0; t_v < 3; ++t_v)
            vertices.col(t_v) = extract_tetra.vertices.at(triangle.vertices.at(t_v)).coordinate;

          exporter.AddPolygon(vertices,
                              {
                                {"original_cell2D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &original_cell2D},
                              });
        }

        exporter.Export(exportFolder + "/extract_triangles.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int t = 0; t < extract_tetra.tetrahedra.size(); ++t)
        {
          const auto& tetrahedron = extract_tetra.tetrahedra.at(t);

          const double original_cell3D = tetrahedron.original_cell3D_index == Gedim::MeshUtilities::None ?
                                           -1.0 :
                                           tetrahedron.original_cell3D_index;

          Eigen::MatrixXd vertices(3, 4);
          std::map<unsigned int, unsigned int> vertex_position;
          for (unsigned int t_v = 0; t_v < 4; ++t_v)
          {
            vertices.col(t_v) = extract_tetra.vertices.at(tetrahedron.vertices.at(t_v)).coordinate;
            vertex_position.insert(std::make_pair(tetrahedron.vertices.at(t_v),
                                                  t_v));
          }

          std::vector<std::vector<unsigned int>> facets(tetrahedron.faces.size());
          for (unsigned int t_f = 0; t_f < tetrahedron.faces.size(); ++t_f)
          {
            const auto& triangle = extract_tetra.triangles.at(tetrahedron.faces.at(t_f));
            facets.at(t_f).resize(3);
            for (unsigned int t_f_v = 0; t_f_v < 3; ++t_f_v)
              facets.at(t_f).at(t_f_v) = vertex_position.at(triangle.vertices.at(t_f_v));
          }

          exporter.AddPolyhedron(vertices,
                                 facets,
                                 {
                                   {"original_cell3D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &original_cell3D},
                                 });
        }

        exporter.Export(exportFolder + "/extract_tetrahedra.vtu");
      }

    }

    Gedim::MeshMatrices tetrahedra_mesh_data;
    Gedim::MeshMatricesDAO tetrahedra_mesh(tetrahedra_mesh_data);
    Gedim::MeshUtilities::TetrahedralMeshDataToMeshDAO(extract_tetra,
                                                       mesh,
                                                       tetrahedra_mesh);

    {
      const std::string mesh_folder = exportFolder + "/tetrahedra_mesh";
      Gedim::Output::CreateFolder(mesh_folder);
      Gedim::External::MeshUtilities::ExportMeshToVTU(tetrahedra_mesh, mesh_folder, "mesh");
    }

  }
} // namespace GedimUnitTesting

#endif // __TEST_meshUtilities3D_ExtractTetrahedralMesh_H
