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

#ifndef __TEST_meshUtilities3D_ImportRegnFaceMesh_H
#define __TEST_meshUtilities3D_ImportRegnFaceMesh_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshUtilities.hpp"

namespace GedimUnitTesting
{
TEST(TestMeshUtilities, TestMeshUtilities_3D_ImportRegnFaceMesh)
{
    std::string exportFolder = "./Export/TestMeshUtilities/TestMeshUtilities_3D_ImportRegnFaceMesh";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-12;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshUtilities mesh_utilities;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    const std::string import_folder = "./Mesh/regn_face_test_mesh";

    mesh_utilities.ImportRegnFaceMesh(geometry_utilities,
                                      import_folder + "/mesh.node",
                                      import_folder + "/mesh.ele",
                                      mesh);
   
    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");

    Gedim::Output::CreateFolder(exportFolder + "/Mesh");
    mesh_utilities.ExportMeshToCsv(mesh, ';', exportFolder + "/Mesh");

    {
        Gedim::MeshUtilities::CheckMesh3DConfiguration config;
        mesh_utilities.CheckMesh3D(config, geometry_utilities, mesh);
    }

    {
        const auto mesh_geometric_data = mesh_utilities.FillMesh3DGeometricData(geometry_utilities, mesh);

        Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration config;
        mesh_utilities.CheckMeshGeometricData3D(config, geometry_utilities, mesh, mesh_geometric_data);
    }
}
} // namespace GedimUnitTesting

#endif // __TEST_MESH_H
