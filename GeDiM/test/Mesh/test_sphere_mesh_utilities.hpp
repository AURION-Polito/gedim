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

} // namespace GedimUnitTesting

#endif // __TEST_SPHERE_MESH_UTILITIES_H
