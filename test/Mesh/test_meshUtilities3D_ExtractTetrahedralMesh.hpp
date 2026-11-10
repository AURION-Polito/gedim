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
                                                           100,
                                                           mesh,
                                                           10);
    }

    {
      const std::string mesh_folder = exportFolder + "/mesh";
      Gedim::Output::CreateFolder(mesh_folder);
      Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_folder, "mesh");
    }

    const auto cell2Ds_triangulation = Gedim::MeshUtilities::ComputeMeshCell2DsTriangles(geometry_utilities,
                                                                                         Gedim::MeshUtilities::MeshCell2DsTrianglesData::TrianglesTypes::ByBarycenter,
                                                                                         mesh);

    const auto cell3Ds_tetra = Gedim::MeshUtilities::ComputeMeshCell3DsTetrahedra(geometry_utilities,
                                                                                  cell2Ds_triangulation,
                                                                                  mesh);

  }
} // namespace GedimUnitTesting

#endif // __TEST_meshUtilities3D_ExtractTetrahedralMesh_H
