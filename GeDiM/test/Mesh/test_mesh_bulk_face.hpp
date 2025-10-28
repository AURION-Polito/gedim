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

#ifndef __test_mesh_bulk_face_H
#define __test_mesh_bulk_face_H

#include "GeometryUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "SphereMeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numeric>

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

TEST(TestBulkFaceMesh, TestCreateCircle)
{

    std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateCircle";
    Gedim::Output::CreateFolder(exportFolder);

    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities mesh_utilities;
    const Gedim::SphereMeshUtilities sphere_mesh_utilities(geometryUtilities, mesh_utilities);

    {
        const Eigen::Vector3d center = Eigen::Vector3d::Zero();
        const double diameter = 1.0;
        const unsigned int num_points = 4;
        Eigen::MatrixXd vertices = sphere_mesh_utilities.circle(center, diameter, num_points);

        vector<unsigned int> vertexMarkers(vertices.cols());
        vector<unsigned int> edgeMarkers(vertices.cols());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), vertices.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        mesh_utilities.Mesh2DFromPolygon(vertices, vertexMarkers, edgeMarkers, mesh);

        Gedim::MeshUtilities::CheckMesh2DConfiguration config;
        mesh_utilities.CheckMesh2D(config, geometryUtilities, mesh);

        // Export to VTK
        {
            Gedim::VTKUtilities vtkUtilities;
            vtkUtilities.AddPolygon(vertices);
            vtkUtilities.Export(exportFolder + "/unit_circle.vtu", Gedim::VTKUtilities::Ascii);
        }
    }

    {
        Eigen::Vector3d center = Eigen::Vector3d::Zero();
        center << 1.5, 2.4, 0.0;
        const double diameter = 2.0;
        const unsigned int num_points = 30;
        Eigen::MatrixXd vertices = sphere_mesh_utilities.circle(center, diameter, num_points);

        vector<unsigned int> vertexMarkers(vertices.cols());
        vector<unsigned int> edgeMarkers(vertices.cols());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), vertices.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        mesh_utilities.Mesh2DFromPolygon(vertices, vertexMarkers, edgeMarkers, mesh);

        Gedim::MeshUtilities::CheckMesh2DConfiguration config;
        mesh_utilities.CheckMesh2D(config, geometryUtilities, mesh);

        // Export to VTK
        {
            Gedim::VTKUtilities vtkUtilities;
            vtkUtilities.AddPolygon(vertices);
            vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
        }
    }
}

TEST(TestBulkFaceMesh, TestCreateTriangularMeshCircle)
{

    std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateTriangularMeshCircle";
    Gedim::Output::CreateFolder(exportFolder);

    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities mesh_utilities;
    const Gedim::SphereMeshUtilities sphere_mesh_utilities(geometryUtilities, mesh_utilities);

    const Eigen::Vector3d center = Eigen::Vector3d::Zero();
    const double diameter = 1.0;
    const unsigned int num_points = 100;
    Eigen::MatrixXd vertices = sphere_mesh_utilities.circle(center, diameter, num_points);

#if ENABLE_TRIANGLE == 0
    throw std::runtime_error("Triangle library not active");
#endif
    const double max_relative_area = 0.1;
    const double max_cell_area = 2.0 * std::numbers::pi * max_relative_area;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
    mesh_utilities.CreateTriangularMesh(vertices, max_cell_area, mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    mesh_utilities.CheckMesh2D(config, geometryUtilities, mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "CircleDomain");
}

} // namespace GedimUnitTesting
#endif
