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
#include "VTKUtilities.hpp"
#include "VtkMeshExporter.hpp"
#include "VoroInterface.hpp"

namespace GedimUnitTesting
{
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

      Gedim::External::MeshUtilities::CreatePolyhedralMesh(geometry_utilities,
                                                           cube.Vertices,
                                                           cube.Edges,
                                                           cube.Faces,
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

//      {
//        Gedim::VTKUtilities exporter;

//        exporter.AddPolygons(mesh.Cell0DsCoordinates(),
//                             mesh.Cell2DsVertices(),
//                             {{"time", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(result.vertices_time.size()), result.vertices_time.data()}});

//        exporter.Export(exportFolder + "/cell2Ds_time.vtu");
//      }
//      {
//        Gedim::VTKUtilities exporter;

//        exporter.AddPolyhedrons(mesh.Cell0DsCoordinates(),
//                             mesh.Cell3DsFacesVertices(),
//                             {{"time", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(result.vertices_time.size()), result.vertices_time.data()}});

//        exporter.Export(exportFolder + "/cell3Ds_time.vtu");
//      }

    }

  }
} // namespace GedimUnitTesting

#endif // __TEST_meshUtilities3D_ExtractTetrahedralMesh_H
