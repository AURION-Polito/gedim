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

#ifndef __TEST_meshUtilities3D_Mesh3DFromPolyhedra_H
#define __TEST_meshUtilities3D_Mesh3DFromPolyhedra_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"

namespace GedimUnitTesting
{

TEST(TestMeshUtilities, TestMeshUtilities_3D_Mesh3DFromPolyhedra)
{
    std::string exportFolder = "./Export/TestMeshUtilities/TestMeshUtilities_3D_Mesh3DFromPolyhedra";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshUtilities mesh_utilities;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    std::vector<Gedim::GeometryUtilities::Polyhedron> polyhedra(4);
    polyhedra.at(0) = geometry_utilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
    polyhedra.at(1) = geometry_utilities.CreateCubeWithOrigin(Eigen::Vector3d(1.0, 0.0, 0.0), 1.0);
    polyhedra.at(2) = geometry_utilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 1.0, 0.0), 1.0);
    polyhedra.at(3) = geometry_utilities.CreateCubeWithOrigin(Eigen::Vector3d(1.0, 1.0, 0.0), 1.0);

    {
        Gedim::VTKUtilities exporter;

        for (const auto &polyhedron : polyhedra)
        {
            exporter.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
        }

        exporter.Export(exportFolder + "/polyhedra.vtu");
    }

    mesh_utilities.Mesh3DFromPolyhedra(geometry_utilities, polyhedra, mesh);

    ASSERT_EQ(18, mesh.Cell0DTotalNumber());
    ASSERT_EQ(33, mesh.Cell1DTotalNumber());
    ASSERT_EQ(20, mesh.Cell2DTotalNumber());
    ASSERT_EQ(4, mesh.Cell3DTotalNumber());

    mesh_utilities.ComputeCell2DCell3DNeighbours(mesh);

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

#endif // __TEST_meshUtilities3D_Mesh3DFromPolyhedra_H
