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

#ifndef __TEST_VORO_INTERFACE_H
#define __TEST_VORO_INTERFACE_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "VoroInterface.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

TEST(TestVoroInterface, TestVoroInterface2D)
{

#if ENABLE_VORO == 0
    GTEST_SKIP_("Voro module not activated.");
#endif

    std::string exportFolder = "./Export/TestVoro";
    Gedim::Output::CreateFolder(exportFolder);

    std::string exportFolderTest = exportFolder + "/Voro2D";
    Gedim::Output::CreateFolder(exportFolderTest);

    Gedim::MeshUtilities mesh_utilities;
    Gedim::MeshUtilities::MeshGeometricData3D mesh_geometric_data;
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1.0e-10;
    geometry_utilities_config.Tolerance2D = 1.0e-12;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Eigen::MatrixXd vertices = Eigen::MatrixXd::Zero(3, 4);
    vertices.row(0) << 0.0, 2.0, 2.0, 0.0; // x
    vertices.row(1) << 0.0, 0.0, 1.0, 1.0; // y
    vertices.row(2) << 0.0, 0.0, 0.0, 0.0; // z

    {
        Gedim::VTKUtilities vtk_utilities;
        vtk_utilities.AddPolygon(vertices);
        vtk_utilities.Export(exportFolderTest + "/Domain.vtu");
    }

    Gedim::VoroInterface voro_interface(geometry_utilities);

    const unsigned int numIterations = 5;
    const unsigned int random_seed = 5;

    for (unsigned int num_cells = 8000; num_cells < 40000; num_cells++)
    {

        std::cout << num_cells << std::endl;

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);

        Eigen::MatrixXd VoronoiPoints = voro_interface.GenerateRandomPoints(vertices, num_cells, random_seed);
        voro_interface.GenerateVoronoiTassellations2D(vertices, numIterations, VoronoiPoints, mesh);

        mesh_utilities.ExportMeshToVTU(mesh, exportFolderTest, "Mesh");

        {
            Gedim::VTKUtilities vtk_utilities;
            vtk_utilities.AddPoints(VoronoiPoints);
            vtk_utilities.Export(exportFolderTest + "/Points.vtu");
        }

        Gedim::MeshUtilities::CheckMesh2DConfiguration config;

        mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);
        mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
        {
            Gedim::Output::Assert(mesh.Cell1DNumberNeighbourCell2D(e) > 0);

            unsigned int count_n = 0;
            for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell2D(e); n++)
            {
                if (!mesh.Cell1DHasNeighbourCell2D(e, n))
                    continue;

                const unsigned int cell2DIndex = mesh.Cell1DNeighbourCell2D(e, n);
                const unsigned int cell2DNumEdges = mesh.Cell2DNumberEdges(cell2DIndex);

                // check edge orientation
                const unsigned int cell2DEdgeIndex = mesh.Cell2DFindEdge(cell2DIndex, e);
                const unsigned int edgeOrigin = mesh.Cell2DVertex(cell2DIndex, (cell2DEdgeIndex + 1) % cell2DNumEdges);
                const unsigned int edgeEnd = mesh.Cell2DVertex(cell2DIndex, cell2DEdgeIndex);

                Gedim::Output::Assert((mesh.Cell2DFindEdgeByExtremes(cell2DIndex, edgeOrigin, edgeEnd) == cell2DEdgeIndex &&
                                       mesh.Cell2DFindEdgeByExtremes(cell2DIndex, edgeEnd, edgeOrigin) == cell2DNumEdges) ||
                                      (mesh.Cell2DFindEdgeByExtremes(cell2DIndex, edgeOrigin, edgeEnd) == cell2DNumEdges &&
                                       mesh.Cell2DFindEdgeByExtremes(cell2DIndex, edgeEnd, edgeOrigin) == cell2DEdgeIndex));

                count_n++;
            }

            Gedim::Output::Assert((mesh.Cell1DMarker(e) != 0 && count_n != 2) || (mesh.Cell1DMarker(e) == 0 && count_n != 1));
        }
    }
}

TEST(TestVoroInterface, TestVoroInterface3D)
{

#if ENABLE_VORO == 0
    GTEST_SKIP_("Voro module not activated.");
#endif

    // Gedim::MeshUtilities mesh_utilities;
    // Gedim::MeshUtilities::MeshGeometricData3D mesh_geometric_data;
    // Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    // geometry_utilities_config.Tolerance1D = 1.0e-12;
    // geometry_utilities_config.Tolerance3D = 1.0e-14;
    // Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    // Eigen::MatrixXd vertices = Eigen::MatrixXd::Zero(3, 8);
    // vertices.row(0) << 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0; // x
    // vertices.row(1) << 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0; // y
    // vertices.row(2) << 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0; // z up
    // Eigen::MatrixXi edges = Eigen::MatrixXi::Zero(2, 12);
    // edges.row(0) << 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3;
    // edges.row(1) << 1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7;
    // std::vector<Eigen::MatrixXi> faces;
    // faces.assign(6, Eigen::MatrixXi::Zero(2, 4));
    // faces[0] << 3, 2, 1, 0, 2, 1, 0, 3;
    // faces[1] << 4, 5, 6, 7, 4, 5, 6, 7;
    // faces[2] << 0, 1, 5, 4, 0, 9, 4, 8;
    // faces[3] << 1, 2, 6, 5, 1, 10, 5, 9;
    // faces[4] << 2, 3, 7, 6, 2, 11, 6, 10;
    // faces[5] << 3, 0, 4, 7, 3, 8, 7, 11;

    // Gedim::VoroInterface voro_interface(geometry_utilities);

    // std::string exportFolder = "./Export/TestVoro";
    // Gedim::Output::CreateFolder(exportFolder);

    // std::string exportFolderTest = exportFolder + "/Voro3D";
    // Gedim::Output::CreateFolder(exportFolderTest);

    // const unsigned int numIterations = 5;
    // const unsigned int voro_seed = 6;
    // const unsigned int num_cells = 1000;

    // Gedim::MeshMatrices mesh_data;
    // Gedim::MeshMatricesDAO mesh(mesh_data);

    // Eigen::MatrixXd VoronoiPoints = voro_interface.GenerateRandomPoints(vertices, num_cells, voro_seed);
    // voro_interface.GenerateVoronoiTassellations3D(vertices, edges, faces, numIterations, VoronoiPoints, mesh);

    // mesh_utilities.ExportMeshToVTU(mesh, exportFolderTest, "Mesh");

    // {
    //     Gedim::VTKUtilities vtk_utilities;
    //     vtk_utilities.AddPoints(VoronoiPoints);
    //     vtk_utilities.Export(exportFolderTest + "/Points.vtu");
    // }

    // Gedim::MeshUtilities::CheckMesh3DConfiguration config;
    // mesh_utilities.CheckMesh3D(config, geometry_utilities, mesh);
}

} // namespace GedimUnitTesting

#endif
