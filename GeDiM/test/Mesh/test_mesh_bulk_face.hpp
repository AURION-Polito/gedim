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
    const double max_relative_area = 0.001;
    const double max_cell_area = 2.0 * std::numbers::pi * max_relative_area;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
    mesh_utilities.CreateTriangularMesh(vertices, max_cell_area, mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    mesh_utilities.CheckMesh2D(config, geometryUtilities, mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "CircleDomain");
}

TEST(TestBulkFaceMesh, TestCreateIntersectionMesh)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateIntersectionMesh";
    Gedim::Output::CreateFolder(exportFolder);

    const auto domain_2D = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(), 1.0);
    const double domain_2D_area = 1.0;

    {
        std::vector<double> domain_2D_id(1, 1.0);
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(domain_2D,
                            {{"domain_id",
                              Gedim::VTPProperty::Formats::Cells,
                              static_cast<unsigned int>(domain_2D_id.size()),
                              domain_2D_id.data()}});
        exporter.Export(exportFolder + "/domain_2D.vtu");
    }

    const unsigned int num_domains_1D = 5;
    Eigen::MatrixXd vertices(3, num_domains_1D);
    vertices.col(0) << 1.0 / 2.0, 1.0 / 4.0, 0.0;
    vertices.col(1) << 5.0 / 8.0, 3.0 / 8.0, 0.0;
    vertices.col(2) << 3.0 / 4.0, 1.0 / 2.0, 0.0;
    vertices.col(3) << 1.0 / 2.0, 3.0 / 4.0, 0.0;
    vertices.col(4) << 1.0 / 4.0, 1.0 / 2.0, 0.0;

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(vertices);
        vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
    }

    Gedim::MeshUtilities mesh_utilities;
    // working meshes
    // 0.05, 0.01, 0.005, 0.001 (a = 0.25 OK, a = 0.1 NO)
    const double domain_2D_max_area = 0.005 * domain_2D_area; // 0.005

    Gedim::MeshMatrices domain_2D_mesh_data;
    Gedim::MeshMatricesDAO domain_2D_mesh(domain_2D_mesh_data);

    const vector<double> baseMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0, true);
    const vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0, true);

    Gedim::MeshMatrices original_mesh_data;
    Gedim::MeshMatricesDAO original_mesh(original_mesh_data);
    mesh_utilities.CreateRectangleMesh(domain_2D.col(0),
                                       domain_2D.col(1) - domain_2D.col(0),
                                       domain_2D.col(3) - domain_2D.col(0),
                                       baseMeshCurvilinearCoordinates,
                                       heightMeshCurvilinearCoordinates,
                                       domain_2D_mesh);

    // mesh_utilities.CreateTriangularMesh(domain_2D, domain_2D_max_area, domain_2D_mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(domain_2D_mesh);

    {
        mesh_utilities.ExportMeshToVTU(domain_2D_mesh, exportFolder, "original_mesh");
    }

    mesh_utilities.CreatePolygonIntersectionMesh(geometry_utilities, vertices, domain_2D_mesh);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    const auto filter_data = mesh_utilities.FilterActiveMesh(domain_2D_mesh);
    mesh_utilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, domain_2D_mesh, mesh);

    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    {
        mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");
    }

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
}

TEST(TestBulkFaceMesh, TestCreateIntersectionMesh_1)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateIntersectionMesh_1";
    Gedim::Output::CreateFolder(exportFolder);

    const auto domain_2D = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(), 1.0);
    const double domain_2D_area = 1.0;

    {
        std::vector<double> domain_2D_id(1, 1.0);
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(domain_2D,
                            {{"domain_id",
                              Gedim::VTPProperty::Formats::Cells,
                              static_cast<unsigned int>(domain_2D_id.size()),
                              domain_2D_id.data()}});
        exporter.Export(exportFolder + "/domain_2D.vtu");
    }

    const unsigned int num_domains_1D = 5;
    Eigen::MatrixXd vertices(3, num_domains_1D);
    vertices.col(0) << 5.0 / 8.0, 3.0 / 8.0, 0.0;
    vertices.col(1) << 3.0 / 4.0, 1.0 / 2.0, 0.0;
    vertices.col(2) << 1.0 / 2.0, 3.0 / 4.0, 0.0;
    vertices.col(3) << 1.0 / 4.0, 1.0 / 2.0, 0.0;
    vertices.col(4) << 1.0 / 2.0, 1.0 / 4.0, 0.0;

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(vertices);
        vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
    }

    Gedim::MeshUtilities mesh_utilities;
    // working meshes
    // 0.05, 0.01, 0.005, 0.001 (a = 0.25 OK, a = 0.1 NO)
    const double domain_2D_max_area = 0.005 * domain_2D_area; // 0.005

    Gedim::MeshMatrices domain_2D_mesh_data;
    Gedim::MeshMatricesDAO domain_2D_mesh(domain_2D_mesh_data);

    const vector<double> baseMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0, true);
    const vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0, true);

    Gedim::MeshMatrices original_mesh_data;
    Gedim::MeshMatricesDAO original_mesh(original_mesh_data);
    mesh_utilities.CreateRectangleMesh(domain_2D.col(0),
                                       domain_2D.col(1) - domain_2D.col(0),
                                       domain_2D.col(3) - domain_2D.col(0),
                                       baseMeshCurvilinearCoordinates,
                                       heightMeshCurvilinearCoordinates,
                                       domain_2D_mesh);

    // mesh_utilities.CreateTriangularMesh(domain_2D, domain_2D_max_area, domain_2D_mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(domain_2D_mesh);

    {
        mesh_utilities.ExportMeshToVTU(domain_2D_mesh, exportFolder, "original_mesh");
    }

    mesh_utilities.CreatePolygonIntersectionMesh(geometry_utilities, vertices, domain_2D_mesh);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    const auto filter_data = mesh_utilities.FilterActiveMesh(domain_2D_mesh);
    mesh_utilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, domain_2D_mesh, mesh);

    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    {
        mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");
    }

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
}

TEST(TestBulkFaceMesh, TestCreateIntersectionMesh_2)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateIntersectionMesh_2";
    Gedim::Output::CreateFolder(exportFolder);

    const auto domain_2D = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(), 1.0);
    const double domain_2D_area = 1.0;

    {
        std::vector<double> domain_2D_id(1, 1.0);
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(domain_2D,
                            {{"domain_id",
                              Gedim::VTPProperty::Formats::Cells,
                              static_cast<unsigned int>(domain_2D_id.size()),
                              domain_2D_id.data()}});
        exporter.Export(exportFolder + "/domain_2D.vtu");
    }

    const unsigned int num_domains_1D = 5;
    Eigen::MatrixXd vertices(3, num_domains_1D);
    vertices.col(0) << 1.0 / 4.0, 1.0 / 4.0, 0.0;
    vertices.col(1) << 3.0 / 8.0, 1.0 / 4.0, 0.0;
    vertices.col(2) << 3.0 / 4.0, 1.0 / 4.0, 0.0;
    vertices.col(3) << 3.0 / 4.0, 3.0 / 4.0, 0.0;
    vertices.col(4) << 1.0 / 4.0, 3.0 / 4.0, 0.0;

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(vertices);
        vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
    }

    Gedim::MeshUtilities mesh_utilities;
    // working meshes
    // 0.05, 0.01, 0.005, 0.001 (a = 0.25 OK, a = 0.1 NO)
    const double domain_2D_max_area = 0.005 * domain_2D_area; // 0.005

    Gedim::MeshMatrices domain_2D_mesh_data;
    Gedim::MeshMatricesDAO domain_2D_mesh(domain_2D_mesh_data);

    const vector<double> baseMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0, true);
    const vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0, true);

    Gedim::MeshMatrices original_mesh_data;
    Gedim::MeshMatricesDAO original_mesh(original_mesh_data);
    mesh_utilities.CreateRectangleMesh(domain_2D.col(0),
                                       domain_2D.col(1) - domain_2D.col(0),
                                       domain_2D.col(3) - domain_2D.col(0),
                                       baseMeshCurvilinearCoordinates,
                                       heightMeshCurvilinearCoordinates,
                                       domain_2D_mesh);

    // mesh_utilities.CreateTriangularMesh(domain_2D, domain_2D_max_area, domain_2D_mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(domain_2D_mesh);

    {
        mesh_utilities.ExportMeshToVTU(domain_2D_mesh, exportFolder, "original_mesh");
    }

    mesh_utilities.CreatePolygonIntersectionMesh(geometry_utilities, vertices, domain_2D_mesh);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    const auto filter_data = mesh_utilities.FilterActiveMesh(domain_2D_mesh);
    mesh_utilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, domain_2D_mesh, mesh);

    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    {
        mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");
    }

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
}

TEST(TestBulkFaceMesh, TestCreateIntersectionMesh_3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateIntersectionMesh_3";
    Gedim::Output::CreateFolder(exportFolder);

    const auto domain_2D = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(), 1.0);
    const double domain_2D_area = 1.0;

    {
        std::vector<double> domain_2D_id(1, 1.0);
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(domain_2D,
                            {{"domain_id",
                              Gedim::VTPProperty::Formats::Cells,
                              static_cast<unsigned int>(domain_2D_id.size()),
                              domain_2D_id.data()}});
        exporter.Export(exportFolder + "/domain_2D.vtu");
    }

    const unsigned int num_domains_1D = 6;
    Eigen::MatrixXd vertices(3, num_domains_1D);
    vertices.col(0) << 1.0 / 4.0, 3.0 / 8.0, 0.0;
    vertices.col(1) << 1.0 / 4.0, 1.0 / 4.0, 0.0;
    vertices.col(2) << 3.0 / 8.0, 1.0 / 4.0, 0.0;
    vertices.col(3) << 3.0 / 4.0, 1.0 / 4.0, 0.0;
    vertices.col(4) << 3.0 / 4.0, 3.0 / 4.0, 0.0;
    vertices.col(5) << 1.0 / 4.0, 3.0 / 4.0, 0.0;

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(vertices);
        vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
    }

    Gedim::MeshUtilities mesh_utilities;
    // working meshes
    // 0.05, 0.01, 0.005, 0.001 (a = 0.25 OK, a = 0.1 NO)
    const double domain_2D_max_area = 0.005 * domain_2D_area; // 0.005

    Gedim::MeshMatrices domain_2D_mesh_data;
    Gedim::MeshMatricesDAO domain_2D_mesh(domain_2D_mesh_data);

    const vector<double> baseMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0, true);
    const vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0, true);

    Gedim::MeshMatrices original_mesh_data;
    Gedim::MeshMatricesDAO original_mesh(original_mesh_data);
    mesh_utilities.CreateRectangleMesh(domain_2D.col(0),
                                       domain_2D.col(1) - domain_2D.col(0),
                                       domain_2D.col(3) - domain_2D.col(0),
                                       baseMeshCurvilinearCoordinates,
                                       heightMeshCurvilinearCoordinates,
                                       domain_2D_mesh);

    // mesh_utilities.CreateTriangularMesh(domain_2D, domain_2D_max_area, domain_2D_mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(domain_2D_mesh);

    {
        mesh_utilities.ExportMeshToVTU(domain_2D_mesh, exportFolder, "original_mesh");
    }

    mesh_utilities.CreatePolygonIntersectionMesh(geometry_utilities, vertices, domain_2D_mesh);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    const auto filter_data = mesh_utilities.FilterActiveMesh(domain_2D_mesh);
    mesh_utilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, domain_2D_mesh, mesh);

    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    {
        mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");
    }

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
}

TEST(TestBulkFaceMesh, TestCreateIntersectionMesh_5)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);
    Gedim::MeshUtilities mesh_utilities;

    const std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateIntersectionMesh_5";
    Gedim::Output::CreateFolder(exportFolder);

    const auto domain_2D = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(), 1.0);
    const double domain_2D_area = 1.0;

    {
        std::vector<double> domain_2D_id(1, 1.0);
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(domain_2D,
                            {{"domain_id",
                              Gedim::VTPProperty::Formats::Cells,
                              static_cast<unsigned int>(domain_2D_id.size()),
                              domain_2D_id.data()}});
        exporter.Export(exportFolder + "/domain_2D.vtu");
    }

    const unsigned int num_domains_1D = 3;
    Eigen::MatrixXd vertices(3, num_domains_1D);
    vertices.col(0) << 1.0 / 4.0, 6.0 / 8.0, 0.0;
    vertices.col(1) << 6.0 / 8.0, 1.0 / 4.0, 0.0;
    vertices.col(2) << 6.0 / 8.0, 6.0 / 8.0, 0.0;

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(vertices);
        vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
    }

    // working meshes
    // 0.05, 0.01, 0.005, 0.001 (a = 0.25 OK, a = 0.1 NO)
    const double domain_2D_max_area = 0.005 * domain_2D_area; // 0.005

    Gedim::MeshMatrices domain_2D_mesh_data;
    Gedim::MeshMatricesDAO domain_2D_mesh(domain_2D_mesh_data);

    const vector<double> baseMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(2 + 1, 0.0, 1.0, true);
    const vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(2 + 1, 0.0, 1.0, true);

    Gedim::MeshMatrices original_mesh_data;
    Gedim::MeshMatricesDAO original_mesh(original_mesh_data);
    mesh_utilities.CreateRectangleMesh(domain_2D.col(0),
                                       domain_2D.col(1) - domain_2D.col(0),
                                       domain_2D.col(3) - domain_2D.col(0),
                                       baseMeshCurvilinearCoordinates,
                                       heightMeshCurvilinearCoordinates,
                                       domain_2D_mesh);

    // mesh_utilities.CreateTriangularMesh(domain_2D, domain_2D_max_area, domain_2D_mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(domain_2D_mesh);

    {
        mesh_utilities.ExportMeshToVTU(domain_2D_mesh, exportFolder, "original_mesh");
    }

    mesh_utilities.CreatePolygonIntersectionMesh(geometry_utilities, vertices, domain_2D_mesh);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    const auto filter_data = mesh_utilities.FilterActiveMesh(domain_2D_mesh);
    mesh_utilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, domain_2D_mesh, mesh);

    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    {
        mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");
    }

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    config.Cell2D_CheckConvexity = false;
    mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
}

TEST(TestBulkFaceMesh, TestCreateIntersectionMesh_6)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);
    Gedim::MeshUtilities mesh_utilities;

    const std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateIntersectionMesh_6";
    Gedim::Output::CreateFolder(exportFolder);

    const auto domain_2D = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(), 1.0);
    const double domain_2D_area = 1.0;

    {
        std::vector<double> domain_2D_id(1, 1.0);
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(domain_2D,
                            {{"domain_id",
                              Gedim::VTPProperty::Formats::Cells,
                              static_cast<unsigned int>(domain_2D_id.size()),
                              domain_2D_id.data()}});
        exporter.Export(exportFolder + "/domain_2D.vtu");
    }

    const Gedim::SphereMeshUtilities sphere_mesh_utilities(geometry_utilities, mesh_utilities);
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    center << 0.5, 0.5, 0.0;
    const double diameter = 0.5;
    const unsigned int num_points = 50;
    Eigen::MatrixXd vertices = sphere_mesh_utilities.circle(center, diameter, num_points);

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(vertices);
        vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
    }

    // working meshes
    // 0.05, 0.01, 0.005, 0.001 (a = 0.25 OK, a = 0.1 NO)
    const double domain_2D_max_area = 0.005 * domain_2D_area; // 0.005

    Gedim::MeshMatrices domain_2D_mesh_data;
    Gedim::MeshMatricesDAO domain_2D_mesh(domain_2D_mesh_data);

    const vector<double> baseMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(2 + 1, 0.0, 1.0, true);
    const vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(2 + 1, 0.0, 1.0, true);

    Gedim::MeshMatrices original_mesh_data;
    Gedim::MeshMatricesDAO original_mesh(original_mesh_data);
    mesh_utilities.CreateRectangleMesh(domain_2D.col(0),
                                       domain_2D.col(1) - domain_2D.col(0),
                                       domain_2D.col(3) - domain_2D.col(0),
                                       baseMeshCurvilinearCoordinates,
                                       heightMeshCurvilinearCoordinates,
                                       domain_2D_mesh);

    // mesh_utilities.CreateTriangularMesh(domain_2D, domain_2D_max_area, domain_2D_mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(domain_2D_mesh);

    {
        mesh_utilities.ExportMeshToVTU(domain_2D_mesh, exportFolder, "original_mesh");
    }

    mesh_utilities.CreatePolygonIntersectionMesh(geometry_utilities, vertices, domain_2D_mesh);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    const auto filter_data = mesh_utilities.FilterActiveMesh(domain_2D_mesh);
    mesh_utilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, domain_2D_mesh, mesh);

    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    {
        mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");
    }

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    config.Cell2D_CheckConvexity = false;
    mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
}

TEST(TestBulkFaceMesh, TestCreateIntersectionMesh_4)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);
    Gedim::MeshUtilities mesh_utilities;

    const std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateIntersectionMesh_4";
    Gedim::Output::CreateFolder(exportFolder);

    const auto domain_2D = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(), 1.0);
    const double domain_2D_area = 1.0;

    {
        std::vector<double> domain_2D_id(1, 1.0);
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(domain_2D,
                            {{"domain_id",
                              Gedim::VTPProperty::Formats::Cells,
                              static_cast<unsigned int>(domain_2D_id.size()),
                              domain_2D_id.data()}});
        exporter.Export(exportFolder + "/domain_2D.vtu");
    }

    const Gedim::SphereMeshUtilities sphere_mesh_utilities(geometry_utilities, mesh_utilities);
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    center << 0.5, 0.5, 0.0;
    const double diameter = 0.5;
    const unsigned int num_points = 50;
    Eigen::MatrixXd shifted_vertices = sphere_mesh_utilities.circle(center, diameter, num_points);
    Eigen::MatrixXd vertices = Eigen::MatrixXd::Zero(3, num_points);
    vertices << shifted_vertices.rightCols(22), shifted_vertices.leftCols(28);

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(vertices);
        vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
    }

    // working meshes
    // 0.05, 0.01, 0.005, 0.001 (a = 0.25 OK, a = 0.1 NO)
    const double domain_2D_max_area = 0.005 * domain_2D_area; // 0.005

    Gedim::MeshMatrices domain_2D_mesh_data;
    Gedim::MeshMatricesDAO domain_2D_mesh(domain_2D_mesh_data);

    const vector<double> baseMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(2 + 1, 0.0, 1.0, true);
    const vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(2 + 1, 0.0, 1.0, true);

    Gedim::MeshMatrices original_mesh_data;
    Gedim::MeshMatricesDAO original_mesh(original_mesh_data);
    mesh_utilities.CreateRectangleMesh(domain_2D.col(0),
                                       domain_2D.col(1) - domain_2D.col(0),
                                       domain_2D.col(3) - domain_2D.col(0),
                                       baseMeshCurvilinearCoordinates,
                                       heightMeshCurvilinearCoordinates,
                                       domain_2D_mesh);

    // mesh_utilities.CreateTriangularMesh(domain_2D, domain_2D_max_area, domain_2D_mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(domain_2D_mesh);

    {
        mesh_utilities.ExportMeshToVTU(domain_2D_mesh, exportFolder, "original_mesh");
    }

    mesh_utilities.CreatePolygonIntersectionMesh(geometry_utilities, vertices, domain_2D_mesh);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    const auto filter_data = mesh_utilities.FilterActiveMesh(domain_2D_mesh);
    mesh_utilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, domain_2D_mesh, mesh);

    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    {
        mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");
    }

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    config.Cell2D_CheckConvexity = false;
    mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
}

TEST(TestBulkFaceMesh, TestCreateIntersectionMesh_7)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);
    Gedim::MeshUtilities mesh_utilities;

    const std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateIntersectionMesh_7";
    Gedim::Output::CreateFolder(exportFolder);

    const auto domain_2D = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(), 1.0);
    const double domain_2D_area = 1.0;

    {
        std::vector<double> domain_2D_id(1, 1.0);
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(domain_2D,
                            {{"domain_id",
                              Gedim::VTPProperty::Formats::Cells,
                              static_cast<unsigned int>(domain_2D_id.size()),
                              domain_2D_id.data()}});
        exporter.Export(exportFolder + "/domain_2D.vtu");
    }

    const Gedim::SphereMeshUtilities sphere_mesh_utilities(geometry_utilities, mesh_utilities);
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    center << 0.52, 0.64, 0.0;
    const double diameter = 0.3;
    const unsigned int num_points = 50;
    Eigen::MatrixXd vertices = sphere_mesh_utilities.circle(center, diameter, num_points);

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(vertices);
        vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
    }

    // working meshes
    std::vector<double> area_relativa = {0.25, 0.1, 0.05, 0.01, 0.005, 0.001};

    for (unsigned int a = 0; a < area_relativa.size(); a++)
    {
        const double domain_2D_max_area = area_relativa[a] * domain_2D_area; // 0.005

        Gedim::MeshMatrices domain_2D_mesh_data;
        Gedim::MeshMatricesDAO domain_2D_mesh(domain_2D_mesh_data);

        Gedim::MeshMatrices original_mesh_data;
        Gedim::MeshMatricesDAO original_mesh(original_mesh_data);

        mesh_utilities.CreateTriangularMesh(domain_2D, domain_2D_max_area, domain_2D_mesh);
        mesh_utilities.ComputeCell1DCell2DNeighbours(domain_2D_mesh);

        {
            mesh_utilities.ExportMeshToVTU(domain_2D_mesh, exportFolder, "original_mesh" + to_string(a));
        }

        mesh_utilities.CreatePolygonIntersectionMesh(geometry_utilities, vertices, domain_2D_mesh);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);

        const auto filter_data = mesh_utilities.FilterActiveMesh(domain_2D_mesh);
        mesh_utilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, domain_2D_mesh, mesh);

        mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

        {
            mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh" + to_string(a));
        }

        Gedim::MeshUtilities::CheckMesh2DConfiguration config;
        config.Cell2D_CheckConvexity = false;
        mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
    }
}

} // namespace GedimUnitTesting
#endif
