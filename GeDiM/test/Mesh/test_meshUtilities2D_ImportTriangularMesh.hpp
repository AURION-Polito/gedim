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

#ifndef __TEST_meshUtilities2D_ImportTriangularMesh_H
#define __TEST_meshUtilities2D_ImportTriangularMesh_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshUtilities.hpp"

namespace GedimUnitTesting
{
TEST(TestMeshUtilities, TestMeshUtilities_2D_ImportTriangularMesh)
{
    std::string exportFolder = "./Export/TestMeshUtilities/TestMeshUtilities_2D_ImportTriangularMesh";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-12;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshUtilities mesh_utilities;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    const std::string import_folder = "./Mesh/Triangular_Test_Mesh";

    mesh_utilities.ImportTriangularMesh(geometry_utilities,
                                        import_folder + "/Cell0Ds.csv",
                                        import_folder + "/Cell2Ds.csv",
                                        import_folder + "/Cell2DsMarker.csv",
                                        ',',
                                        mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");

    Gedim::Output::CreateFolder(exportFolder + "/Mesh");
    mesh_utilities.ExportMeshToCsv(mesh, ';', exportFolder + "/Mesh");

    {
        Gedim::MeshUtilities::CheckMesh2DConfiguration config;
        config.Cell1D_CheckNeighbours = false;
        mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
    }

    {
        const auto mesh_geometric_data = mesh_utilities.FillMesh2DGeometricData(geometry_utilities, mesh);

        Gedim::MeshUtilities::CheckMeshGeometricData2DConfiguration config;
        mesh_utilities.CheckMeshGeometricData2D(config, geometry_utilities, mesh, mesh_geometric_data);
    }
}
} // namespace GedimUnitTesting

#endif // __TEST_MESH_H
