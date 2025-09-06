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

#ifndef __TEST_PLATONIC_SOLID_H
#define __TEST_PLATONIC_SOLID_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numeric>

#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "PlatonicSolid.hpp"
#include "VTKUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

TEST(TestPlatonicSolid, TestTetrahedron)
{
    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;

    const Gedim::PlatonicSolid platonicSolid = Gedim::PlatonicSolid(geometryUtilities, meshUtilities);

    const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.tetrahedron();

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

    meshUtilities.ComputeCell2DCell3DNeighbours(mesh);
    const Gedim::MeshUtilities::MeshGeometricData3D geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
    const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
    meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestTetrahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
        vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
        meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh");
    }
}

TEST(TestPlatonicSolid, TestDualTetrahedron)
{
    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;

    const Gedim::PlatonicSolid platonicSolid = Gedim::PlatonicSolid(geometryUtilities, meshUtilities);

    const Gedim::GeometryUtilities::Polyhedron tetrahedron = platonicSolid.tetrahedron();

    const Gedim::GeometryUtilities::Polyhedron dual = platonicSolid.dual_polyhedron(tetrahedron);

    vector<unsigned int> vertexMarkers(dual.Vertices.cols());
    vector<unsigned int> edgeMarkers(dual.Edges.cols());
    vector<unsigned int> faceMarkers(dual.Faces.size());

    std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
    std::iota(edgeMarkers.begin(), edgeMarkers.end(), dual.Vertices.cols() + 1);
    std::iota(faceMarkers.begin(), faceMarkers.end(), dual.Vertices.cols() + dual.Edges.cols() + 1);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
    meshUtilities.Mesh3DFromPolyhedron(dual.Vertices, dual.Edges, dual.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);

    Gedim::MeshUtilities::CheckMesh3DConfiguration config;
    meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);

    meshUtilities.ComputeCell2DCell3DNeighbours(mesh);
    const Gedim::MeshUtilities::MeshGeometricData3D geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
    const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
    meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestDualTetrahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(dual.Vertices, dual.Edges, dual.Faces);
        vtpUtilities.Export(exportFolder + "/Dual.vtu", Gedim::VTKUtilities::Ascii);
        meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh");
    }
}

TEST(TestPlatonicSolid, TestHexahedron)
{
    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;

    const Gedim::PlatonicSolid platonicSolid = Gedim::PlatonicSolid(geometryUtilities, meshUtilities);

    const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.hexahedron();

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

    meshUtilities.ComputeCell2DCell3DNeighbours(mesh);
    const Gedim::MeshUtilities::MeshGeometricData3D geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
    const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
    meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestHexahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
        vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
        meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh");
    }
}

TEST(TestPlatonicSolid, TestOctahedron)
{
    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;

    const Gedim::PlatonicSolid platonicSolid = Gedim::PlatonicSolid(geometryUtilities, meshUtilities);

    const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.octahedron();

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

    meshUtilities.ComputeCell2DCell3DNeighbours(mesh);
    const Gedim::MeshUtilities::MeshGeometricData3D geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
    const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
    meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestOctahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
        vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
        meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh");
    }
}

TEST(TestPlatonicSolid, TestIcosahedron)
{
    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;

    const Gedim::PlatonicSolid platonicSolid = Gedim::PlatonicSolid(geometryUtilities, meshUtilities);

    const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.icosahedron();

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

    meshUtilities.ComputeCell2DCell3DNeighbours(mesh);
    const Gedim::MeshUtilities::MeshGeometricData3D geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
    const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
    meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestIcosahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
        vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
        meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh");
    }
}

TEST(TestPlatonicSolid, TestDodecahedron)
{
    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;

    const Gedim::PlatonicSolid platonicSolid = Gedim::PlatonicSolid(geometryUtilities, meshUtilities);

    const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.dodecahedron();

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

    meshUtilities.ComputeCell2DCell3DNeighbours(mesh);
    const Gedim::MeshUtilities::MeshGeometricData3D geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
    const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
    meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestDodecahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
        vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
        meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh");
    }
}

TEST(TestPlatonicSolid, TestTriangulateI)
{
    std::string exportFolder = "./Export/TestPlatonicSolid/TestTriangulateI";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-14;
    geometryUtilitiesConfig.Tolerance3D = 1.0e-10;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;
    const Gedim::PlatonicSolid platonicSolid(geometryUtilities, meshUtilities);
    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.tetrahedron();

        for (unsigned int i = 1; i < 4; i++)
        {
            Gedim::MeshMatrices mesh_data;
            Gedim::MeshMatricesDAO mesh(mesh_data);

            platonicSolid.first_class_geodesic_polyhedron(polyhedron, i, mesh);

            // Export to VTK
            {
                meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_tetrahedron_Triangle_" + to_string(i));
            }

            meshUtilities.ComputeCell1DCell2DNeighbours(mesh);
            meshUtilities.ComputeCell2DCell3DNeighbours(mesh);

            Gedim::MeshUtilities::CheckMesh3DConfiguration config;
            meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);
            const Gedim::MeshUtilities::MeshGeometricData3D geometricData =
                meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
            const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
            meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);
        }
    }

    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.octahedron();

        for (unsigned int i = 1; i < 4; i++)
        {
            Gedim::MeshMatrices mesh_data;
            Gedim::MeshMatricesDAO mesh(mesh_data);

            platonicSolid.first_class_geodesic_polyhedron(polyhedron, i, mesh);

            // Export to VTK
            {
                meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_octahedron_Triangle_" + to_string(i));
            }

            meshUtilities.ComputeCell1DCell2DNeighbours(mesh);
            meshUtilities.ComputeCell2DCell3DNeighbours(mesh);

            Gedim::MeshUtilities::CheckMesh3DConfiguration config;
            meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);
            const Gedim::MeshUtilities::MeshGeometricData3D geometricData =
                meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
            const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
            meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);
        }
    }

    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.icosahedron();

        for (unsigned int i = 1; i < 4; i++)
        {
            Gedim::MeshMatrices mesh_data;
            Gedim::MeshMatricesDAO mesh(mesh_data);

            platonicSolid.first_class_geodesic_polyhedron(polyhedron, i, mesh);

            // Export to VTK
            {
                meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_icosahedron_Triangle_" + to_string(i));
            }

            meshUtilities.ComputeCell1DCell2DNeighbours(mesh);
            meshUtilities.ComputeCell2DCell3DNeighbours(mesh);

            Gedim::MeshUtilities::CheckMesh3DConfiguration config;
            meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);
            const Gedim::MeshUtilities::MeshGeometricData3D geometricData =
                meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
            const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
            meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);
        }
    }
}

TEST(TestPlatonicSolid, TestTriangulateII)
{
    std::string exportFolder = "./Export/TestPlatonicSolid/TestTriangulateII";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-14;
    geometryUtilitiesConfig.Tolerance3D = 1.0e-10;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;

    const Gedim::PlatonicSolid platonicSolid(geometryUtilities, meshUtilities);

    {
        Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.tetrahedron();

        for (unsigned int i = 1; i < 4; i++)
        {
            Gedim::MeshMatrices mesh_data;
            Gedim::MeshMatricesDAO mesh(mesh_data);

            platonicSolid.second_class_geodesic_polyhedron(polyhedron, i, mesh);
            auto a = mesh.Cell1DsExtremes();

            const unsigned int num_vertices = polyhedron.Vertices.cols() + polyhedron.Edges.cols() * (2 * i - 1) +
                                              polyhedron.Faces.size() * ((i * i - i) * 1.5 + 1);
            ASSERT_TRUE(num_vertices == mesh.Cell0DTotalNumber());

            const unsigned int num_edges = polyhedron.Edges.cols() * 2 * i + 3 * polyhedron.Faces.size() * ((3 * i * i + i) * 0.5);
            ASSERT_TRUE(num_edges == mesh.Cell1DTotalNumber());

            const unsigned int num_faces = polyhedron.Faces.size() * 3 * (i * i + i);
            ASSERT_TRUE(num_faces == mesh.Cell2DTotalNumber());

            ASSERT_TRUE(num_faces + num_vertices == num_edges + 2);

            // Export to VTK
            {
                meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_tetrahedron_Triangle_" + to_string(i));
            }

            meshUtilities.ComputeCell1DCell2DNeighbours(mesh);
            meshUtilities.ComputeCell2DCell3DNeighbours(mesh);

            Gedim::MeshUtilities::CheckMesh3DConfiguration config;
            meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);
            const Gedim::MeshUtilities::MeshGeometricData3D geometricData =
                meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
            const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
            meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);
        }
    }

    {
        Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.octahedron();

        for (unsigned int i = 1; i < 4; i++)
        {
            Gedim::MeshMatrices mesh_data;
            Gedim::MeshMatricesDAO mesh(mesh_data);

            platonicSolid.second_class_geodesic_polyhedron(polyhedron, i, mesh);
            auto a = mesh.Cell1DsExtremes();

            const unsigned int num_vertices = polyhedron.Vertices.cols() + polyhedron.Edges.cols() * (2 * i - 1) +
                                              polyhedron.Faces.size() * ((i * i - i) * 1.5 + 1);
            ASSERT_TRUE(num_vertices == mesh.Cell0DTotalNumber());

            const unsigned int num_edges = polyhedron.Edges.cols() * 2 * i + 3 * polyhedron.Faces.size() * ((3 * i * i + i) * 0.5);
            ASSERT_TRUE(num_edges == mesh.Cell1DTotalNumber());

            const unsigned int num_faces = polyhedron.Faces.size() * 3 * (i * i + i);
            ASSERT_TRUE(num_faces == mesh.Cell2DTotalNumber());

            ASSERT_TRUE(num_faces + num_vertices == num_edges + 2);

            // Export to VTK
            {
                meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_octahedron_Triangle_" + to_string(i));
            }

            meshUtilities.ComputeCell1DCell2DNeighbours(mesh);
            meshUtilities.ComputeCell2DCell3DNeighbours(mesh);

            Gedim::MeshUtilities::CheckMesh3DConfiguration config;
            meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);
            const Gedim::MeshUtilities::MeshGeometricData3D geometricData =
                meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
            const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
            meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);
        }
    }

    {
        Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.icosahedron();

        for (unsigned int i = 1; i < 4; i++)
        {
            Gedim::MeshMatrices mesh_data;
            Gedim::MeshMatricesDAO mesh(mesh_data);

            platonicSolid.second_class_geodesic_polyhedron(polyhedron, i, mesh);
            auto a = mesh.Cell1DsExtremes();

            const unsigned int num_vertices = polyhedron.Vertices.cols() + polyhedron.Edges.cols() * (2 * i - 1) +
                                              polyhedron.Faces.size() * ((i * i - i) * 1.5 + 1);
            ASSERT_TRUE(num_vertices == mesh.Cell0DTotalNumber());

            const unsigned int num_edges = polyhedron.Edges.cols() * 2 * i + 3 * polyhedron.Faces.size() * ((3 * i * i + i) * 0.5);
            ASSERT_TRUE(num_edges == mesh.Cell1DTotalNumber());

            const unsigned int num_faces = polyhedron.Faces.size() * 3 * (i * i + i);
            ASSERT_TRUE(num_faces == mesh.Cell2DTotalNumber());

            ASSERT_TRUE(num_faces + num_vertices == num_edges + 2);

            // Export to VTK
            {
                meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_icosahedron_Triangle_" + to_string(i));
            }

            meshUtilities.ComputeCell1DCell2DNeighbours(mesh);
            meshUtilities.ComputeCell2DCell3DNeighbours(mesh);

            Gedim::MeshUtilities::CheckMesh3DConfiguration config;
            meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);
            const Gedim::MeshUtilities::MeshGeometricData3D geometricData =
                meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh);
            const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration configuration;
            meshUtilities.CheckMeshGeometricData3D(configuration, geometryUtilities, mesh, geometricData);
        }
    }
}

TEST(TestPlatonicSolid, TestGeodesicPolyhedron)
{
    std::string exportFolder = "./Export/TestPlatonicSolid/TestGeodesicPolyhedron";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-14;
    geometryUtilitiesConfig.Tolerance3D = 1.0e-10;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;
    const Gedim::PlatonicSolid platonicSolid(geometryUtilities, meshUtilities);
    for (unsigned int i = 1; i < 4; i++)
    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.geodesic_polyhedron(3, 3, i, 0);

        vector<unsigned int> vertexMarkers(polyhedron.Vertices.cols());
        vector<unsigned int> edgeMarkers(polyhedron.Edges.cols());
        vector<unsigned int> faceMarkers(polyhedron.Faces.size());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), polyhedron.Vertices.cols() + 1);
        std::iota(faceMarkers.begin(), faceMarkers.end(), polyhedron.Vertices.cols() + polyhedron.Edges.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        meshUtilities.Mesh3DFromPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);

        // Export to VTK
        {
            Gedim::VTKUtilities vtpUtilities;

            //  original polyhedron
            vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
            vtpUtilities.Export(exportFolder + "/Geodesic_Tetrahedron_" + to_string(i) + ".vtu", Gedim::VTKUtilities::Ascii);
            meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_Tetrahedron_" + to_string(i));
        }
    }

    for (unsigned int i = 1; i < 4; i++)
    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.geodesic_polyhedron(3, 4, i, 0);

        vector<unsigned int> vertexMarkers(polyhedron.Vertices.cols());
        vector<unsigned int> edgeMarkers(polyhedron.Edges.cols());
        vector<unsigned int> faceMarkers(polyhedron.Faces.size());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), polyhedron.Vertices.cols() + 1);
        std::iota(faceMarkers.begin(), faceMarkers.end(), polyhedron.Vertices.cols() + polyhedron.Edges.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        meshUtilities.Mesh3DFromPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);

        // Export to VTK
        {
            Gedim::VTKUtilities vtpUtilities;

            //  original polyhedron
            vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
            vtpUtilities.Export(exportFolder + "/Geodesic_Octahedron_" + to_string(i) + ".vtu", Gedim::VTKUtilities::Ascii);
            meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_Octahedron_" + to_string(i));
        }
    }

    for (unsigned int i = 1; i < 4; i++)
    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.geodesic_polyhedron(3, 5, i, 0);

        vector<unsigned int> vertexMarkers(polyhedron.Vertices.cols());
        vector<unsigned int> edgeMarkers(polyhedron.Edges.cols());
        vector<unsigned int> faceMarkers(polyhedron.Faces.size());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), polyhedron.Vertices.cols() + 1);
        std::iota(faceMarkers.begin(), faceMarkers.end(), polyhedron.Vertices.cols() + polyhedron.Edges.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        meshUtilities.Mesh3DFromPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);

        // Export to VTK
        {
            Gedim::VTKUtilities vtpUtilities;

            //  original polyhedron
            vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
            vtpUtilities.Export(exportFolder + "/Geodesic_Icosahedron_" + to_string(i) + ".vtu", Gedim::VTKUtilities::Ascii);
            meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_Icosahedron_" + to_string(i));
        }
    }
}

TEST(TestPlatonicSolid, TestGoldbergPolyhedron)
{
    std::string exportFolder = "./Export/TestPlatonicSolid/TestGoldbergPolyhedron";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-14;
    geometryUtilitiesConfig.Tolerance3D = 1.0e-10;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;
    const Gedim::PlatonicSolid platonicSolid(geometryUtilities, meshUtilities);
    for (unsigned int i = 1; i < 4; i++)
    {

        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.goldberg_polyhedron(3, 3, i, 0);

        vector<unsigned int> vertexMarkers(polyhedron.Vertices.cols());
        vector<unsigned int> edgeMarkers(polyhedron.Edges.cols());
        vector<unsigned int> faceMarkers(polyhedron.Faces.size());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), polyhedron.Vertices.cols() + 1);
        std::iota(faceMarkers.begin(), faceMarkers.end(), polyhedron.Vertices.cols() + polyhedron.Edges.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        meshUtilities.Mesh3DFromPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);

        // Export to VTK
        {
            Gedim::VTKUtilities vtpUtilities;

            //  original polyhedron
            vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
            vtpUtilities.Export(exportFolder + "/Goldberg_Tetrahedron_" + to_string(i) + ".vtu", Gedim::VTKUtilities::Ascii);
            meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_Tetrahedron_" + to_string(i));
        }
    }

    for (unsigned int i = 1; i < 4; i++)
    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.goldberg_polyhedron(4, 3, i, 0);

        vector<unsigned int> vertexMarkers(polyhedron.Vertices.cols());
        vector<unsigned int> edgeMarkers(polyhedron.Edges.cols());
        vector<unsigned int> faceMarkers(polyhedron.Faces.size());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), polyhedron.Vertices.cols() + 1);
        std::iota(faceMarkers.begin(), faceMarkers.end(), polyhedron.Vertices.cols() + polyhedron.Edges.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        meshUtilities.Mesh3DFromPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);

        // Export to VTK
        {
            Gedim::VTKUtilities vtpUtilities;

            //  original polyhedron
            vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
            vtpUtilities.Export(exportFolder + "/Goldberg_Octahedron_" + to_string(i) + ".vtu", Gedim::VTKUtilities::Ascii);
            meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_Octahedron_" + to_string(i));
        }
    }

    for (unsigned int i = 1; i < 4; i++)
    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.goldberg_polyhedron(5, 3, i, 0);

        vector<unsigned int> vertexMarkers(polyhedron.Vertices.cols());
        vector<unsigned int> edgeMarkers(polyhedron.Edges.cols());
        vector<unsigned int> faceMarkers(polyhedron.Faces.size());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), polyhedron.Vertices.cols() + 1);
        std::iota(faceMarkers.begin(), faceMarkers.end(), polyhedron.Vertices.cols() + polyhedron.Edges.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        meshUtilities.Mesh3DFromPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);

        // Export to VTK
        {
            Gedim::VTKUtilities vtpUtilities;

            //  original polyhedron
            vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
            vtpUtilities.Export(exportFolder + "/Goldberg_Icosahedron_" + to_string(i) + ".vtu", Gedim::VTKUtilities::Ascii);
            meshUtilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_Icosahedron_" + to_string(i));
        }
    }
}

} // namespace GedimUnitTesting

#endif // __TEST_PLATONIC_SOLID_H
