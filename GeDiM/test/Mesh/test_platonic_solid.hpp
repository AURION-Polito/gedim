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

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestTetrahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

        vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
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

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestDualTetrahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(dual.Vertices, dual.Edges, dual.Faces);

        vtpUtilities.Export(exportFolder + "/Dual.vtu", Gedim::VTKUtilities::Ascii);
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

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestHexahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

        vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
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

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestOctahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

        vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
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

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestIcosahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

        vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
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

    // Export to VTK
    std::string exportFolder = "./Export/TestPlatonicSolid/TestDodecahedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

        vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
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

    for (unsigned int i = 1; i < 3; i++)
    {
        const Gedim::PlatonicSolid platonicSolid(geometryUtilities, meshUtilities);

        const Gedim::GeometryUtilities::Polyhedron ico = platonicSolid.tetrahedron();
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.first_class_geodesic_polyhedron(ico, i);

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
        }

        Gedim::MeshUtilities::CheckMesh3DConfiguration config;
        config.Cell3D_CheckConvexity = false;
        meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);
    }

    for (unsigned int i = 1; i < 3; i++)
    {
        const Gedim::PlatonicSolid platonicSolid(geometryUtilities, meshUtilities);

        const Gedim::GeometryUtilities::Polyhedron ico = platonicSolid.icosahedron();
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.first_class_geodesic_polyhedron(ico, i);

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
        }

        Gedim::MeshUtilities::CheckMesh3DConfiguration config;
        config.Cell3D_CheckConvexity = true;
        meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);
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

    for (unsigned int i = 1; i < 3; i++)
    {
        const Gedim::PlatonicSolid platonicSolid(geometryUtilities, meshUtilities);

        const Gedim::GeometryUtilities::Polyhedron ico = platonicSolid.tetrahedron();
        const Gedim::GeometryUtilities::Polyhedron dual = platonicSolid.first_class_geodesic_polyhedron(ico, i);
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.goldberg_polyhedron(dual);

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
        }

        Gedim::MeshUtilities::CheckMesh3DConfiguration config;
        config.Cell3D_CheckConvexity = false;
        meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);
    }

    for (unsigned int i = 1; i < 3; i++)
    {
        const Gedim::PlatonicSolid platonicSolid(geometryUtilities, meshUtilities);

        const Gedim::GeometryUtilities::Polyhedron ico = platonicSolid.icosahedron();
        const Gedim::GeometryUtilities::Polyhedron dual = platonicSolid.first_class_geodesic_polyhedron(ico, i);
        const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.goldberg_polyhedron(dual);

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
        }

        Gedim::MeshUtilities::CheckMesh3DConfiguration config;
        config.Cell3D_CheckConvexity = true;
        meshUtilities.CheckMesh3D(config, geometryUtilities, mesh);
    }
}

} // namespace GedimUnitTesting

#endif // __TEST_PLATONIC_SOLID_H
