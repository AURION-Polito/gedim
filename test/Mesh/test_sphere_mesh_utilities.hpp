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

#ifndef __TEST_SPHERE_MESH_UTILITIES_H
#define __TEST_SPHERE_MESH_UTILITIES_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numeric>

#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "SphereMeshUtilities.hpp"
#include "VTKUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
TEST(TestSphereMeshUtilities, TestUVSphere)
{
    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;
    const Gedim::SphereMeshUtilities sphereMeshUtilities(geometryUtilities, meshUtilities);

    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron = sphereMeshUtilities.uv_sphere(4, 2);

        vector<unsigned int> vertexMarkers(polyhedron.Vertices.cols());
        vector<unsigned int> edgeMarkers(polyhedron.Edges.cols());
        vector<unsigned int> faceMarkers(polyhedron.Faces.size());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), polyhedron.Vertices.cols() + 1);
        std::iota(faceMarkers.begin(), faceMarkers.end(), polyhedron.Vertices.cols() + polyhedron.Edges.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        meshUtilities.Mesh3DFromPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);

        Gedim::MeshUtilities::CheckMesh3DConfiguration config;
        meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);

        // Export to VTK
        std::string exportFolder = "./Export/TestSphereMeshUtilities/TestUVSphere";
        Gedim::Output::CreateFolder(exportFolder);

        {
            Gedim::VTKUtilities vtpUtilities;

            //  original polyhedron
            vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

            vtpUtilities.Export(exportFolder + "/ref_0.vtu", Gedim::VTKUtilities::Ascii);
        }
    }

    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron = sphereMeshUtilities.uv_sphere(8, 7);

        vector<unsigned int> vertexMarkers(polyhedron.Vertices.cols());
        vector<unsigned int> edgeMarkers(polyhedron.Edges.cols());
        vector<unsigned int> faceMarkers(polyhedron.Faces.size());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), polyhedron.Vertices.cols() + 1);
        std::iota(faceMarkers.begin(), faceMarkers.end(), polyhedron.Vertices.cols() + polyhedron.Edges.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        meshUtilities.Mesh3DFromPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);

        Gedim::MeshUtilities::CheckMesh3DConfiguration config;
        meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);

        // Export to VTK
        std::string exportFolder = "./Export/TestSphereMeshUtilities/TestUVSphere";
        Gedim::Output::CreateFolder(exportFolder);

        {
            Gedim::VTKUtilities vtpUtilities;

            //  original polyhedron
            vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

            vtpUtilities.Export(exportFolder + "/ref_1.vtu", Gedim::VTKUtilities::Ascii);
        }
    }

    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron = sphereMeshUtilities.uv_sphere(23, 11);

        vector<unsigned int> vertexMarkers(polyhedron.Vertices.cols());
        vector<unsigned int> edgeMarkers(polyhedron.Edges.cols());
        vector<unsigned int> faceMarkers(polyhedron.Faces.size());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), polyhedron.Vertices.cols() + 1);
        std::iota(faceMarkers.begin(), faceMarkers.end(), polyhedron.Vertices.cols() + polyhedron.Edges.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        meshUtilities.Mesh3DFromPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);

        Gedim::MeshUtilities::CheckMesh3DConfiguration config;
        meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);

        // Export to VTK
        std::string exportFolder = "./Export/TestSphereMeshUtilities/TestUVSphere";
        Gedim::Output::CreateFolder(exportFolder);

        {
            Gedim::VTKUtilities vtpUtilities;

            //  original polyhedron
            vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

            vtpUtilities.Export(exportFolder + "/ref_2.vtu", Gedim::VTKUtilities::Ascii);
        }
    }
}

TEST(TestSphereMeshUtilities, TestCreateCircle)
{

    std::string exportFolder = "./Export/TestSphereMeshUtilities/TestCreateCircle";
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

TEST(TestSphereMeshUtilities, TestCreateEllipse)
{

    std::string exportFolder = "./Export/TestSphereMeshUtilities/TestCreateEllipse";
    Gedim::Output::CreateFolder(exportFolder);

    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities mesh_utilities;
    const Gedim::SphereMeshUtilities sphere_mesh_utilities(geometryUtilities, mesh_utilities);

    {
        const Eigen::Vector3d center = Eigen::Vector3d::Zero();
        const double radius_1 = 1.0;
        const double radius_2 = 1.0;
        const unsigned int num_points = 100;
        Eigen::MatrixXd vertices = sphere_mesh_utilities.ellipse(center, radius_1, radius_2, 0.0, num_points);

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

        const double radius_1 = 5.0;
        const double radius_2 = 2.3;
        const double rotation_angle = std::numbers::pi * 0.25;
        const unsigned int num_points = 100;
        Eigen::MatrixXd vertices = sphere_mesh_utilities.ellipse(center, radius_1, radius_2, rotation_angle, num_points);

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
            vtkUtilities.Export(exportFolder + "/ellipse.vtu", Gedim::VTKUtilities::Ascii);
        }
    }
}

} // namespace GedimUnitTesting

#endif // __TEST_SPHERE_MESH_UTILITIES_H
