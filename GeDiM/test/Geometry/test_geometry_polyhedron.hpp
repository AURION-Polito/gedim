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

#ifndef __TEST_GEOMETRY_POLYHEDRON_H
#define __TEST_GEOMETRY_POLYHEDRON_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numeric>

#include "GeometryUtilities.hpp"
#include "GraphUtilities.hpp"
#include "MapTetrahedron.hpp"
#include "MeshMatricesDAO.hpp"
#include "PlatonicSolid.hpp"
#include "QuadratureData.hpp"
#include "Quadrature_Gauss3D_Tetrahedron_PositiveWeights.hpp"
#include "VTKUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronBarycenter)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check cube barycenter
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));
            ASSERT_TRUE(geometryUtilities.PointsAreCoincident(geometryUtilities.PolyhedronBarycenter(cube.Vertices),
                                                              Eigen::Vector3d(0.5, 0.5, 0.5)));
        }

        // check tetrahedron barycenter
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));

            ASSERT_TRUE(geometryUtilities.PointsAreCoincident(geometryUtilities.PolyhedronBarycenter(tetrahedron.Vertices),
                                                              Eigen::Vector3d(0.25, 0.25, 0.25)));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronEdgeTangents)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check cube edge tangents
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));
            Eigen::MatrixXd expectedEdgeTangents(3, 12);
            expectedEdgeTangents.col(0) << 1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
            expectedEdgeTangents.col(1) << 0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00;
            expectedEdgeTangents.col(2) << -1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
            expectedEdgeTangents.col(3) << 0.0000000000000000e+00, -1.0000000000000000e+00, 0.0000000000000000e+00;
            expectedEdgeTangents.col(4) << 1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
            expectedEdgeTangents.col(5) << 0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00;
            expectedEdgeTangents.col(6) << -1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
            expectedEdgeTangents.col(7) << 0.0000000000000000e+00, -1.0000000000000000e+00, 0.0000000000000000e+00;
            expectedEdgeTangents.col(8) << 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;
            expectedEdgeTangents.col(9) << 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;
            expectedEdgeTangents.col(10) << 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;
            expectedEdgeTangents.col(11) << 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;

            ASSERT_EQ(geometryUtilities.PolyhedronEdgeTangents(cube.Vertices, cube.Edges), expectedEdgeTangents);
        }

        // check tetrahedron face vertices
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));

            Eigen::MatrixXd expectedEdgeTangents(3, 6);
            expectedEdgeTangents.col(0) << 1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
            expectedEdgeTangents.col(1) << 0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00;
            expectedEdgeTangents.col(2) << -1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00;
            expectedEdgeTangents.col(3) << 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;
            expectedEdgeTangents.col(4) << -1.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;
            expectedEdgeTangents.col(5) << 0.0000000000000000e+00, -1.0000000000000000e+00, 1.0000000000000000e+00;

            ASSERT_EQ(geometryUtilities.PolyhedronEdgeTangents(tetrahedron.Vertices, tetrahedron.Edges), expectedEdgeTangents);
            ;
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceVertices)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check cube face vertices
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));
            vector<Eigen::MatrixXd> expectedFaceVertices(6, Eigen::MatrixXd(3, 4));
            expectedFaceVertices[0] << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00;
            expectedFaceVertices[1] << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00;
            expectedFaceVertices[2] << 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00;
            expectedFaceVertices[3] << 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00;
            expectedFaceVertices[4] << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00;
            expectedFaceVertices[5] << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00;

            ASSERT_EQ(geometryUtilities.PolyhedronFaceVertices(cube.Vertices, cube.Faces), expectedFaceVertices);
        }

        // check tetrahedron face vertices
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));

            vector<Eigen::MatrixXd> expectedFaceVertices(4, Eigen::MatrixXd(3, 3));

            expectedFaceVertices[0] << 0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00;
            expectedFaceVertices[1] << 0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 1.0000000000000000e+00;
            expectedFaceVertices[2] << 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 1.0000000000000000e+00;
            expectedFaceVertices[3] << 1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 1.0000000000000000e+00;

            ASSERT_EQ(geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices, tetrahedron.Faces), expectedFaceVertices);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceNormals)
{
    try
    {
        std::string exportFolder = "./Export/TestPolyhedron_TestPolyhedronFaceNormals";
        Gedim::Output::CreateFolder(exportFolder);

        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check cube face normals
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));

            const Eigen::Vector3d barycenter(0.5, 0.5, 0.5);
            const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices, cube.Faces);
            const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);
            const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
            const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);
            const vector<Eigen::Matrix3d> faceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices, faceNormals, faceTranslations);

            const vector<Eigen::MatrixXd> face2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(faceVertices, faceTranslations, faceRotationMatrices);

            vector<Eigen::Vector3d> expectedFaceNormals(6);
            expectedFaceNormals[0] << +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
            expectedFaceNormals[1] << +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
            expectedFaceNormals[2] << +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceNormals[3] << +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceNormals[4] << +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceNormals[5] << +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;

            ASSERT_EQ(faceNormals, expectedFaceNormals);
            ASSERT_EQ(geometryUtilities.PolyhedronFaceNormalDirections(cube.Vertices,
                                                                       cube.Edges,
                                                                       cube.Faces,
                                                                       faceVertices,
                                                                       faceBarycenters,
                                                                       face2DVertices,
                                                                       faceNormals,
                                                                       faceTranslations,
                                                                       faceRotationMatrices),
                      vector<bool>({false, true, false, true, true, false}));
            ASSERT_EQ(geometryUtilities.PolyhedronFaceNormalDirections(faceVertices, barycenter, faceNormals),
                      vector<bool>({false, true, false, true, true, false}));
        }

        // check tetrahedron face normals
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));
            const Eigen::Vector3d barycenter(0.25, 0.25, 0.25);
            const vector<Eigen::MatrixXd> faceVertices =
                geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices, tetrahedron.Faces);
            const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);
            const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
            const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);
            const vector<Eigen::Matrix3d> faceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices, faceNormals, faceTranslations);

            const vector<Eigen::MatrixXd> face2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(faceVertices, faceTranslations, faceRotationMatrices);

            vector<Eigen::Vector3d> expectedFaceNormals(4);
            expectedFaceNormals[0] << +0.0000000000000000e+00, +0.0000000000000000e+00, 1.0000000000000000e+00;
            expectedFaceNormals[1] << +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceNormals[2] << +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceNormals[3] << +5.7735026918962584e-01, +5.7735026918962584e-01, +5.7735026918962584e-01;

            ASSERT_EQ(faceNormals, expectedFaceNormals);
            ASSERT_EQ(geometryUtilities.PolyhedronFaceNormalDirections(tetrahedron.Vertices,
                                                                       tetrahedron.Edges,
                                                                       tetrahedron.Faces,
                                                                       faceVertices,
                                                                       faceBarycenters,
                                                                       face2DVertices,
                                                                       faceNormals,
                                                                       faceTranslations,
                                                                       faceRotationMatrices),
                      vector<bool>({false, true, false, true}));
            ASSERT_EQ(geometryUtilities.PolyhedronFaceNormalDirections(faceVertices, barycenter, faceNormals),
                      vector<bool>({false, true, false, true}));
        }

        // check tetrahedron 2 face normals
        {
            Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
            geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
            Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

            Eigen::MatrixXd tetraVertices(3, 4);
            tetraVertices.row(0) << 4.0386589600299999e-01, 4.4444444444400000e-01, 3.3333333333300003e-01, 4.0386589600299999e-01;
            tetraVertices.row(1) << 1.0135526475499999e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
            tetraVertices.row(2) << 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0135526475499999e-01;
            Eigen::MatrixXi tetraEdges(2, 6);
            tetraEdges.row(0) << 0, 0, 3, 1, 2, 3;
            tetraEdges.row(1) << 2, 1, 0, 3, 1, 2;
            std::vector<Eigen::MatrixXi> tetraFaces(4, Eigen::MatrixXi(2, 3));
            tetraFaces[0].row(0) << 3, 2, 1;
            tetraFaces[0].row(1) << 5, 4, 3;
            tetraFaces[1].row(0) << 3, 0, 1;
            tetraFaces[1].row(1) << 2, 1, 3;
            tetraFaces[2].row(0) << 0, 3, 2;
            tetraFaces[2].row(1) << 2, 5, 0;
            tetraFaces[3].row(0) << 1, 2, 0;
            tetraFaces[3].row(1) << 4, 0, 1;

            const Eigen::Vector3d barycenter = geometryUtilities.PolyhedronBarycenter(tetraVertices);
            const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(tetraVertices, tetraFaces);
            const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);
            const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
            const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);
            const vector<Eigen::Matrix3d> faceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices, faceNormals, faceTranslations);

            const vector<Eigen::MatrixXd> face2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(faceVertices, faceTranslations, faceRotationMatrices);
            const vector<Eigen::Vector3d> face2DCentroid = geometryUtilities.PolyhedronFaceBarycenter(face2DVertices);
            const std::vector<std::vector<unsigned int>> faceTriangulations =
                geometryUtilities.PolyhedronFaceTriangulationsByEarClipping(tetraFaces.size(), face2DVertices);
            const std::vector<std::vector<Eigen::Matrix3d>> faces3DTriangulationPoints =
                geometryUtilities.PolyhedronFaceExtractTriangulationPoints(faceVertices, faceTriangulations);
            std::vector<Eigen::Vector3d> faceInternalPoints(tetraFaces.size());
            for (unsigned int f = 0; f < tetraFaces.size(); f++)
                faceInternalPoints[f] = geometryUtilities.PolygonBarycenter(faces3DTriangulationPoints[f][0]);

            const std::vector<bool> face2DNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(tetraVertices,
                                                                 tetraEdges,
                                                                 tetraFaces,
                                                                 faceVertices,
                                                                 faceBarycenters,
                                                                 face2DVertices,
                                                                 faceNormals,
                                                                 faceTranslations,
                                                                 faceRotationMatrices);
            const string exportTetraFolder = exportFolder + "/Tetra2";
            Gedim::Output::CreateFolder(exportTetraFolder);
            geometryUtilities.ExportPolyhedronToVTU(0,
                                                    tetraVertices,
                                                    tetraEdges,
                                                    tetraFaces,
                                                    {tetraVertices},
                                                    0.0,
                                                    barycenter,
                                                    faceVertices,
                                                    {0.0, 0.0, 0.0, 0.0},
                                                    face2DCentroid,
                                                    faceTranslations,
                                                    faceRotationMatrices,
                                                    faces3DTriangulationPoints,
                                                    faceInternalPoints,
                                                    faceNormals,
                                                    face2DNormalDirections,
                                                    exportTetraFolder);

            ASSERT_EQ(vector<bool>({true, false, false, true}), face2DNormalDirections);
            ASSERT_EQ(geometryUtilities.PolyhedronFaceNormalDirections(faceVertices, barycenter, faceNormals),
                      vector<bool>({true, false, false, true}));
        }

        // check tetrahedron 3 face normals
        {
            Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
            geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
            Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

            Eigen::MatrixXd tetraVertices(3, 4);
            tetraVertices.row(0) << 1.2896219846500001e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
            tetraVertices.row(1) << 0.0000000000000000e+00, 1.2896219846500001e-01, 0.0000000000000000e+00, 0.0000000000000000e+00;
            tetraVertices.row(2) << 2.0216932951100000e-01, 2.0216932951100000e-01, 2.2222222222200000e-01, 1.1111111111100000e-01;
            Eigen::MatrixXi tetraEdges(2, 6);
            tetraEdges.row(0) << 3, 1, 1, 2, 0, 3;
            tetraEdges.row(1) << 0, 0, 3, 1, 2, 2;
            std::vector<Eigen::MatrixXi> tetraFaces(4, Eigen::MatrixXi(2, 3));
            tetraFaces[0].row(0) << 3, 2, 1;
            tetraFaces[0].row(1) << 5, 3, 2;
            tetraFaces[1].row(0) << 3, 1, 0;
            tetraFaces[1].row(1) << 2, 1, 0;
            tetraFaces[2].row(0) << 3, 0, 2;
            tetraFaces[2].row(1) << 0, 4, 5;
            tetraFaces[3].row(0) << 1, 0, 2;
            tetraFaces[3].row(1) << 1, 4, 3;

            const Eigen::Vector3d barycenter = geometryUtilities.PolyhedronBarycenter(tetraVertices);
            const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(tetraVertices, tetraFaces);
            const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);
            const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
            const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);
            const vector<Eigen::Matrix3d> faceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices, faceNormals, faceTranslations);

            const vector<Eigen::MatrixXd> face2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(faceVertices, faceTranslations, faceRotationMatrices);
            const vector<Eigen::Vector3d> face2DCentroid = geometryUtilities.PolyhedronFaceBarycenter(face2DVertices);
            const std::vector<std::vector<unsigned int>> faceTriangulations =
                geometryUtilities.PolyhedronFaceTriangulationsByEarClipping(tetraFaces.size(), face2DVertices);
            const std::vector<std::vector<Eigen::Matrix3d>> faces3DTriangulationPoints =
                geometryUtilities.PolyhedronFaceExtractTriangulationPoints(faceVertices, faceTriangulations);
            std::vector<Eigen::Vector3d> faceInternalPoints(tetraFaces.size());
            for (unsigned int f = 0; f < tetraFaces.size(); f++)
                faceInternalPoints[f] = geometryUtilities.PolygonBarycenter(faces3DTriangulationPoints[f][0]);

            const std::vector<bool> face2DNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(tetraVertices,
                                                                 tetraEdges,
                                                                 tetraFaces,
                                                                 faceVertices,
                                                                 faceBarycenters,
                                                                 face2DVertices,
                                                                 faceNormals,
                                                                 faceTranslations,
                                                                 faceRotationMatrices);
            const string exportTetraFolder = exportFolder + "/Tetra3";
            Gedim::Output::CreateFolder(exportTetraFolder);
            geometryUtilities.ExportPolyhedronToVTU(0,
                                                    tetraVertices,
                                                    tetraEdges,
                                                    tetraFaces,
                                                    {tetraVertices},
                                                    0.0,
                                                    barycenter,
                                                    faceVertices,
                                                    {0.0, 0.0, 0.0, 0.0},
                                                    face2DCentroid,
                                                    faceTranslations,
                                                    faceRotationMatrices,
                                                    faces3DTriangulationPoints,
                                                    faceInternalPoints,
                                                    faceNormals,
                                                    face2DNormalDirections,
                                                    exportTetraFolder);

            ASSERT_EQ(geometryUtilities.PolyhedronFaceNormalDirections(faceVertices, barycenter, faceNormals),
                      vector<bool>({true, true, true, false}));
            ASSERT_EQ(vector<bool>({true, true, true, false}), face2DNormalDirections);
        }

        // check tetrahedron 4 face normals
        {
            Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
            geometryUtilitiesConfig.Tolerance1D = 1.0e-5;
            Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

            Eigen::MatrixXd tetraVertices(3, 4);
            tetraVertices.row(0) << 2.4069810210299999e-01, 2.5752285704000000e-01, 1.0671088379000000e-01, 1.2153801988400000e-01;
            tetraVertices.row(1) << 4.8751808986700002e-01, 3.8626821047799997e-01, 3.7923669784300001e-01, 4.3349091758300001e-01;
            tetraVertices.row(2) << 4.4591160662099999e-01, 4.2369762890099999e-01, 3.8280917906299999e-01, 4.8221637936799999e-01;
            Eigen::MatrixXi tetraEdges(2, 6);
            tetraEdges.row(0) << 0, 0, 0, 1, 2, 2;
            tetraEdges.row(1) << 2, 3, 1, 3, 1, 3;
            std::vector<Eigen::MatrixXi> tetraFaces(4, Eigen::MatrixXi(2, 3));
            tetraFaces[0].row(0) << 3, 2, 1;
            tetraFaces[0].row(1) << 5, 4, 3;
            tetraFaces[1].row(0) << 3, 1, 0;
            tetraFaces[1].row(1) << 3, 2, 1;
            tetraFaces[2].row(0) << 3, 0, 2;
            tetraFaces[2].row(1) << 1, 0, 5;
            tetraFaces[3].row(0) << 0, 2, 1;
            tetraFaces[3].row(1) << 0, 4, 2;

            const Eigen::Vector3d barycenter = geometryUtilities.PolyhedronBarycenter(tetraVertices);
            const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(tetraVertices, tetraFaces);
            const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);
            const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
            const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);
            const vector<Eigen::Matrix3d> faceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices, faceNormals, faceTranslations);

            const vector<Eigen::MatrixXd> face2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(faceVertices, faceTranslations, faceRotationMatrices);
            const vector<Eigen::Vector3d> face2DCentroid = geometryUtilities.PolyhedronFaceBarycenter(face2DVertices);
            const std::vector<std::vector<unsigned int>> faceTriangulations =
                geometryUtilities.PolyhedronFaceTriangulationsByEarClipping(tetraFaces.size(), face2DVertices);
            const std::vector<std::vector<Eigen::Matrix3d>> faces3DTriangulationPoints =
                geometryUtilities.PolyhedronFaceExtractTriangulationPoints(faceVertices, faceTriangulations);
            std::vector<Eigen::Vector3d> faceInternalPoints(tetraFaces.size());
            for (unsigned int f = 0; f < tetraFaces.size(); f++)
                faceInternalPoints[f] = geometryUtilities.PolygonBarycenter(faces3DTriangulationPoints[f][0]);

            const std::vector<bool> face2DNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(tetraVertices,
                                                                 tetraEdges,
                                                                 tetraFaces,
                                                                 faceVertices,
                                                                 faceBarycenters,
                                                                 face2DVertices,
                                                                 faceNormals,
                                                                 faceTranslations,
                                                                 faceRotationMatrices);
            const string exportTetraFolder = exportFolder + "/Tetra4";
            Gedim::Output::CreateFolder(exportTetraFolder);
            geometryUtilities.ExportPolyhedronToVTU(0,
                                                    tetraVertices,
                                                    tetraEdges,
                                                    tetraFaces,
                                                    {tetraVertices},
                                                    0.0,
                                                    barycenter,
                                                    faceVertices,
                                                    {0.0, 0.0, 0.0, 0.0},
                                                    face2DCentroid,
                                                    faceTranslations,
                                                    faceRotationMatrices,
                                                    faces3DTriangulationPoints,
                                                    faceInternalPoints,
                                                    faceNormals,
                                                    face2DNormalDirections,
                                                    exportTetraFolder);

            ASSERT_EQ(geometryUtilities.PolyhedronFaceNormalDirections(faceVertices, barycenter, faceNormals),
                      vector<bool>({true, true, true, false}));
            ASSERT_EQ(vector<bool>({true, true, true, false}), face2DNormalDirections);
        }

        // check tetrahedron 4 face normals
        {
            Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
            geometryUtilitiesConfig.MinTolerance = 1.0e-14;
            geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
            geometryUtilitiesConfig.Tolerance2D = 1.0e-10;
            geometryUtilitiesConfig.Tolerance3D = 1.0e-10;
            Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

            Eigen::MatrixXd tetraVertices(3, 4);
            tetraVertices.row(0) << 2.6689961207925128e-01, 2.9727523153850433e-01, 2.3553366387491570e-01, 2.6218536281077792e-01;
            tetraVertices.row(1) << 2.2030798286073699e-01, 2.3832083723535250e-01, 2.8404384733338389e-01, 1.7441173133044119e-01;
            tetraVertices.row(2) << 2.5263209528783409e-01, 1.8159146436527179e-01, 1.8365509494692039e-01, 1.8786269256096871e-01;
            Eigen::MatrixXi tetraEdges(2, 6);
            tetraEdges.row(0) << 0, 3, 1, 3, 1, 2;
            tetraEdges.row(1) << 2, 0, 0, 1, 2, 3;
            std::vector<Eigen::MatrixXi> tetraFaces(4, Eigen::MatrixXi(2, 3));
            tetraFaces[0].row(0) << 3, 2, 1;
            tetraFaces[0].row(1) << 5, 4, 3;
            tetraFaces[1].row(0) << 3, 0, 1;
            tetraFaces[1].row(1) << 1, 2, 3;
            tetraFaces[2].row(0) << 0, 3, 2;
            tetraFaces[2].row(1) << 1, 5, 0;
            tetraFaces[3].row(0) << 0, 2, 1;
            tetraFaces[3].row(1) << 0, 4, 2;

            const Eigen::Vector3d barycenter = geometryUtilities.PolyhedronBarycenter(tetraVertices);
            const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(tetraVertices, tetraFaces);
            const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);
            const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
            const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);
            const vector<Eigen::Matrix3d> faceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices, faceNormals, faceTranslations);

            const vector<Eigen::MatrixXd> face2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(faceVertices, faceTranslations, faceRotationMatrices);
            const vector<Eigen::Vector3d> face2DCentroid = geometryUtilities.PolyhedronFaceBarycenter(face2DVertices);
            const std::vector<std::vector<unsigned int>> faceTriangulations =
                geometryUtilities.PolyhedronFaceTriangulationsByEarClipping(tetraFaces.size(), face2DVertices);
            const std::vector<std::vector<Eigen::Matrix3d>> faces3DTriangulationPoints =
                geometryUtilities.PolyhedronFaceExtractTriangulationPoints(faceVertices, faceTriangulations);
            std::vector<Eigen::Vector3d> faceInternalPoints(tetraFaces.size());
            for (unsigned int f = 0; f < tetraFaces.size(); f++)
                faceInternalPoints[f] = geometryUtilities.PolygonBarycenter(faces3DTriangulationPoints[f][0]);

            const std::vector<bool> face2DNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(tetraVertices,
                                                                 tetraEdges,
                                                                 tetraFaces,
                                                                 faceVertices,
                                                                 faceBarycenters,
                                                                 face2DVertices,
                                                                 faceNormals,
                                                                 faceTranslations,
                                                                 faceRotationMatrices);
            const string exportTetraFolder = exportFolder + "/Tetra5";
            Gedim::Output::CreateFolder(exportTetraFolder);
            geometryUtilities.ExportPolyhedronToVTU(0,
                                                    tetraVertices,
                                                    tetraEdges,
                                                    tetraFaces,
                                                    {tetraVertices},
                                                    0.0,
                                                    barycenter,
                                                    faceVertices,
                                                    {0.0, 0.0, 0.0, 0.0},
                                                    face2DCentroid,
                                                    faceTranslations,
                                                    faceRotationMatrices,
                                                    faces3DTriangulationPoints,
                                                    faceInternalPoints,
                                                    faceNormals,
                                                    face2DNormalDirections,
                                                    exportTetraFolder);

            ASSERT_EQ(geometryUtilities.PolyhedronFaceNormalDirections(faceVertices, barycenter, faceNormals),
                      vector<bool>({true, false, false, false}));
            ASSERT_EQ(vector<bool>({true, false, false, false}), face2DNormalDirections);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceNormals_Concave)
{
    try
    {
        std::string exportFolder = "./Export/TestPolyhedron_TestPolyhedronFaceNormals_Concave";
        Gedim::Output::CreateFolder(exportFolder);

        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check concave face normals
        {
            Eigen::MatrixXd polygon(3, 4);
            polygon.col(0) << 0.0, 0.0, 0.0;
            polygon.col(1) << -1.0, -1.0, 0.0;
            polygon.col(2) << 1.0, 0.0, 0.0;
            polygon.col(3) << -1.0, 1.0, 0.0;

            const Gedim::GeometryUtilities::Polyhedron concave =
                geometryUtilities.CreatePolyhedronWithExtrusion(polygon, Eigen::Vector3d(0.0, 0.0, 1.0));

            {
                geometryUtilities.ExportPolyhedronToVTU(concave.Vertices, concave.Edges, concave.Faces, exportFolder);
            }

            const std::vector<Eigen::MatrixXd> faceVertices =
                geometryUtilities.PolyhedronFaceVertices(concave.Vertices, concave.Faces);
            const std::vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);

            const unsigned int numFaces = concave.Faces.size();

            const std::vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);

            const std::vector<Eigen::Matrix3d> faceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices, faceNormals, faceTranslations);
            std::vector<Eigen::MatrixXd> face2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(faceVertices, faceTranslations, faceRotationMatrices);

            std::vector<Eigen::MatrixXd> face2DVerticesCCW(numFaces);
            std::vector<std::vector<unsigned int>> face2DsCCW(numFaces);

            for (unsigned int f = 0; f < numFaces; f++)
            {
                const vector<unsigned int> faceConvexHull = geometryUtilities.ConvexHull(face2DVertices[f], false);

                if (geometryUtilities.PolygonOrientation(faceConvexHull) == Gedim::GeometryUtilities::PolygonOrientations::Clockwise)
                {
                    face2DsCCW[f] = geometryUtilities.ChangePolygonOrientation(face2DVertices[f].cols());
                    face2DVerticesCCW[f] = geometryUtilities.ExtractPoints(face2DVertices[f], face2DsCCW[f]);
                }
                else
                {
                    face2DVerticesCCW[f] = face2DVertices[f];

                    face2DsCCW[f].resize(face2DVertices[f].cols());
                    std::iota(face2DsCCW[f].begin(), face2DsCCW[f].end(), 0);
                }
            }

            std::vector<std::vector<unsigned int>> faceTriangulations =
                geometryUtilities.PolyhedronFaceTriangulationsByEarClipping(numFaces, face2DVerticesCCW);
            for (unsigned int f = 0; f < numFaces; f++)
            {
                for (unsigned int nt = 0; nt < faceTriangulations[f].size(); nt++)
                    faceTriangulations[f][nt] = face2DsCCW[f][faceTriangulations[f][nt]];
            }
            const std::vector<std::vector<Eigen::Matrix3d>> facesTriangulationPoints =
                geometryUtilities.PolyhedronFaceExtractTriangulationPoints(faceVertices, faceTriangulations);

            {
                Gedim::VTKUtilities exporter;

                for (unsigned int f = 0; f < facesTriangulationPoints.size(); f++)
                {
                    std::vector<double> faceIndex(1, f);
                    for (unsigned int t = 0; t < facesTriangulationPoints[f].size(); t++)
                    {
                        exporter.AddPolygon(facesTriangulationPoints[f][t],
                                            {{"Face",
                                              Gedim::VTPProperty::Formats::Cells,
                                              static_cast<unsigned int>(faceIndex.size()),
                                              faceIndex.data()}});
                    }
                }

                exporter.Export(exportFolder + "/FacesTriangulationPoints.vtu");
            }

            std::vector<Eigen::Vector3d> faceInternalPoints(concave.Faces.size());
            for (unsigned int f = 0; f < concave.Faces.size(); f++)
                faceInternalPoints[f] = geometryUtilities.PolygonBarycenter(facesTriangulationPoints[f][0]);

            {
                Gedim::VTKUtilities exporter;

                for (unsigned int f = 0; f < faceInternalPoints.size(); f++)
                {
                    std::vector<double> faceIndex(1, f);
                    exporter.AddPoint(faceInternalPoints[f],
                                      {{"Face",
                                        Gedim::VTPProperty::Formats::Cells,
                                        static_cast<unsigned int>(faceIndex.size()),
                                        faceIndex.data()}});
                }

                exporter.Export(exportFolder + "/FacesInternalPoint.vtu");
            }

            {
                Gedim::VTKUtilities exporter;

                for (unsigned int f = 0; f < faceInternalPoints.size(); f++)
                {
                    std::vector<double> faceIndex(1, f);
                    exporter.AddSegment(faceInternalPoints[f],
                                        faceInternalPoints[f] + faceNormals[f],
                                        {{"Face",
                                          Gedim::VTPProperty::Formats::Cells,
                                          static_cast<unsigned int>(faceIndex.size()),
                                          faceIndex.data()}});
                }

                exporter.Export(exportFolder + "/FacesNormal.vtu");
            }

            vector<Eigen::Vector3d> expectedFaceNormals(6);
            expectedFaceNormals[0] << +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
            expectedFaceNormals[1] << +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
            expectedFaceNormals[2] << -7.0710678118654746e-01, +7.0710678118654746e-01, +0.0000000000000000e+00;
            expectedFaceNormals[3] << +4.4721359549995793e-01, -8.9442719099991586e-01, +0.0000000000000000e+00;
            expectedFaceNormals[4] << +4.4721359549995793e-01, +8.9442719099991586e-01, +0.0000000000000000e+00;
            expectedFaceNormals[5] << -7.0710678118654746e-01, -7.0710678118654746e-01, +0.0000000000000000e+00;

            ASSERT_EQ(faceNormals, expectedFaceNormals);
            ASSERT_EQ(geometryUtilities.PolyhedronFaceNormalDirections(concave.Vertices,
                                                                       concave.Edges,
                                                                       concave.Faces,
                                                                       faceVertices,
                                                                       faceInternalPoints,
                                                                       face2DVertices,
                                                                       faceNormals,
                                                                       faceTranslations,
                                                                       faceRotationMatrices),
                      vector<bool>({false, true, true, true, true, true}));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceEdgeDirections)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check cube face edge directions
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));
            vector<vector<bool>> expectedFaceEdgeDirections(6);
            expectedFaceEdgeDirections[0] = vector<bool>({true, true, true, true});
            expectedFaceEdgeDirections[1] = vector<bool>({true, true, true, true});
            expectedFaceEdgeDirections[2] = vector<bool>({false, true, true, false});
            expectedFaceEdgeDirections[3] = vector<bool>({true, true, false, false});
            expectedFaceEdgeDirections[4] = vector<bool>({true, true, false, false});
            expectedFaceEdgeDirections[5] = vector<bool>({false, true, true, false});

            ASSERT_EQ(expectedFaceEdgeDirections,
                      geometryUtilities.PolyhedronFaceEdgeDirections(cube.Vertices, cube.Edges, cube.Faces));
        }

        // check tetrahedron face edge directions
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));

            vector<vector<bool>> expectedFaceEdgeDirections(4);
            expectedFaceEdgeDirections[0] = vector<bool>({true, true, false});
            expectedFaceEdgeDirections[1] = vector<bool>({true, true, false});
            expectedFaceEdgeDirections[2] = vector<bool>({true, true, false});
            expectedFaceEdgeDirections[3] = vector<bool>({true, true, false});

            ASSERT_EQ(
                expectedFaceEdgeDirections,
                geometryUtilities.PolyhedronFaceEdgeDirections(tetrahedron.Vertices, tetrahedron.Edges, tetrahedron.Faces));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceEdgeTangents)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check cube face edge tangents
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));

            const Eigen::MatrixXd edgeTangents = geometryUtilities.PolyhedronEdgeTangents(cube.Vertices, cube.Edges);
            const vector<vector<bool>> faceEdgeDirections =
                geometryUtilities.PolyhedronFaceEdgeDirections(cube.Vertices, cube.Edges, cube.Faces);

            vector<Eigen::MatrixXd> expectedFaceEdgeTangents(6);
            expectedFaceEdgeTangents[0] = (Eigen::MatrixXd(3, 4) << 1, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 0).finished();
            expectedFaceEdgeTangents[1] = (Eigen::MatrixXd(3, 4) << 1, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 0).finished();
            expectedFaceEdgeTangents[2] = (Eigen::MatrixXd(3, 4) << -0, 0, 0, -0, 1, 0, -1, -0, -0, 1, 0, -1).finished();
            expectedFaceEdgeTangents[3] = (Eigen::MatrixXd(3, 4) << 0, 0, -0, -0, 1, 0, -1, -0, 0, 1, -0, -1).finished();
            expectedFaceEdgeTangents[4] = (Eigen::MatrixXd(3, 4) << 1, 0, -1, -0, 0, 0, -0, -0, 0, 1, -0, -1).finished();
            expectedFaceEdgeTangents[5] = (Eigen::MatrixXd(3, 4) << 1, 0, -1, -0, -0, 0, 0, -0, -0, 1, 0, -1).finished();

            ASSERT_EQ(geometryUtilities.PolyhedronFaceEdgeTangents(cube.Vertices, cube.Edges, cube.Faces, faceEdgeDirections, edgeTangents),
                      expectedFaceEdgeTangents);
        }

        // check tetrahedron face edge tangents
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));
            const Eigen::MatrixXd edgeTangents =
                geometryUtilities.PolyhedronEdgeTangents(tetrahedron.Vertices, tetrahedron.Edges);
            const vector<vector<bool>> faceEdgeDirections =
                geometryUtilities.PolyhedronFaceEdgeDirections(tetrahedron.Vertices, tetrahedron.Edges, tetrahedron.Faces);

            vector<Eigen::MatrixXd> expectedFaceEdgeTangents(4);
            expectedFaceEdgeTangents[0] = (Eigen::MatrixXd(3, 3) << 1, -1, -0, 0, 1, -1, 0, 0, -0).finished();
            expectedFaceEdgeTangents[1] = (Eigen::MatrixXd(3, 3) << 1, -1, -0, 0, 0, -0, 0, 1, -1).finished();
            expectedFaceEdgeTangents[2] = (Eigen::MatrixXd(3, 3) << 0, 0, -0, 1, -1, -0, 0, 1, -1).finished();
            expectedFaceEdgeTangents[3] = (Eigen::MatrixXd(3, 3) << -1, 0, 1, 1, -1, -0, 0, 1, -1).finished();

            ASSERT_EQ(geometryUtilities.PolyhedronFaceEdgeTangents(tetrahedron.Vertices, tetrahedron.Edges, tetrahedron.Faces, faceEdgeDirections, edgeTangents),
                      expectedFaceEdgeTangents);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceTranslations)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check cube face translations
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));

            const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices, cube.Faces);

            vector<Eigen::Vector3d> expectedFaceTranslations(6);
            expectedFaceTranslations[0] << +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceTranslations[1] << +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
            expectedFaceTranslations[2] << +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceTranslations[3] << +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceTranslations[4] << +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceTranslations[5] << +0.0000000000000000e+00, +1.0000000000000000e+00, +0.0000000000000000e+00;

            ASSERT_EQ(expectedFaceTranslations, geometryUtilities.PolyhedronFaceTranslations(faceVertices));
        }

        // check tetrahedron face translations
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));
            const vector<Eigen::MatrixXd> faceVertices =
                geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices, tetrahedron.Faces);

            vector<Eigen::Vector3d> expectedFaceTranslations(4);
            expectedFaceTranslations[0] << +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceTranslations[1] << +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceTranslations[2] << +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
            expectedFaceTranslations[3] << +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;

            ASSERT_EQ(expectedFaceTranslations, geometryUtilities.PolyhedronFaceTranslations(faceVertices));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceRotationMatrices)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check cube face rotation matrices
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));

            const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices, cube.Faces);
            const Eigen::Vector3d barycenter = geometryUtilities.PolyhedronBarycenter(cube.Vertices);
            const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
            const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);

            vector<Eigen::Matrix3d> expectedFaceRotationMatrices(6);
            expectedFaceRotationMatrices[0] = (Eigen::Matrix3d() << 9.9999999999999978e-01,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               9.9999999999999978e-01,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               1.0000000000000000e+00)
                                                  .finished();
            expectedFaceRotationMatrices[1] = (Eigen::Matrix3d() << 9.9999999999999978e-01,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               9.9999999999999978e-01,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               1.0000000000000000e+00)
                                                  .finished();
            expectedFaceRotationMatrices[2] = (Eigen::Matrix3d() << 0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               1.0000000000000000e+00,
                                               9.9999999999999989e-01,
                                               -1.1102230246251565e-16,
                                               0.0000000000000000e+00,
                                               1.1102230246251565e-16,
                                               9.9999999999999989e-01,
                                               0.0000000000000000e+00)
                                                  .finished();
            expectedFaceRotationMatrices[3] = (Eigen::Matrix3d() << 0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               1.0000000000000000e+00,
                                               9.9999999999999989e-01,
                                               -1.1102230246251565e-16,
                                               0.0000000000000000e+00,
                                               1.1102230246251565e-16,
                                               9.9999999999999989e-01,
                                               0.0000000000000000e+00)
                                                  .finished();
            expectedFaceRotationMatrices[4] = (Eigen::Matrix3d() << 9.9999999999999978e-01,
                                               0.0000000000000000e+00,
                                               4.9650683064945600e-17,
                                               4.9650683064945576e-17,
                                               2.4825341532472850e-17,
                                               -1.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               9.9999999999999978e-01,
                                               2.4825341532472825e-17)
                                                  .finished();
            expectedFaceRotationMatrices[5] = (Eigen::Matrix3d() << 9.9999999999999978e-01,
                                               0.0000000000000000e+00,
                                               4.9650683064945600e-17,
                                               4.9650683064945576e-17,
                                               2.4825341532472850e-17,
                                               -1.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               9.9999999999999978e-01,
                                               2.4825341532472825e-17)
                                                  .finished();

            ASSERT_EQ(expectedFaceRotationMatrices,
                      geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices, faceNormals, faceTranslations));
        }

        // check tetrahedron face rotation matrices
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));
            const vector<Eigen::MatrixXd> faceVertices =
                geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices, tetrahedron.Faces);

            const Eigen::Vector3d barycenter = geometryUtilities.PolyhedronBarycenter(tetrahedron.Vertices);
            const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
            const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);

            vector<Eigen::Matrix3d> expectedFaceRotationMatrices(4);
            expectedFaceRotationMatrices[0] = (Eigen::Matrix3d() << 1.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               1.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               1.0000000000000000e+00)
                                                  .finished();
            expectedFaceRotationMatrices[1] = (Eigen::Matrix3d() << 1.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               -1.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               1.0000000000000000e+00,
                                               0.0000000000000000e+00)
                                                  .finished();
            expectedFaceRotationMatrices[2] = (Eigen::Matrix3d() << 0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               1.0000000000000000e+00,
                                               1.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               0.0000000000000000e+00,
                                               1.0000000000000000e+00,
                                               0.0000000000000000e+00)
                                                  .finished();
            expectedFaceRotationMatrices[3] = (Eigen::Matrix3d() << -7.0710678118654724e-01,
                                               -4.0824829046386296e-01,
                                               5.7735026918962595e-01,
                                               7.0710678118654724e-01,
                                               -4.0824829046386313e-01,
                                               5.7735026918962573e-01,
                                               0.0000000000000000e+00,
                                               8.1649658092772581e-01,
                                               5.7735026918962651e-01)
                                                  .finished();

            ASSERT_EQ(expectedFaceRotationMatrices,
                      geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices, faceNormals, faceTranslations));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronCoordinateSystem)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check cube coordinate system
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));

            ASSERT_EQ(geometryUtilities.PolyhedronCoordinateSystem(cube.Vertices, cube.Edges),
                      vector<unsigned int>({0, 1, 3, 4}));
        }

        // check tetrahedron face normals
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));
            ASSERT_EQ(geometryUtilities.PolyhedronCoordinateSystem(tetrahedron.Vertices, tetrahedron.Edges),
                      vector<unsigned int>({0, 1, 2, 3}));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceBarycenters)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check cube face baycenters
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));
            const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices, cube.Faces);
            ASSERT_EQ(geometryUtilities.PolyhedronFaceBarycenter(faceVertices),
                      vector<Eigen::Vector3d>({Eigen::Vector3d(0.5, 0.5, 0.0),
                                               Eigen::Vector3d(0.5, 0.5, 1.0),
                                               Eigen::Vector3d(0.0, 0.5, 0.5),
                                               Eigen::Vector3d(1.0, 0.5, 0.5),
                                               Eigen::Vector3d(0.5, 0.0, 0.5),
                                               Eigen::Vector3d(0.5, 1.0, 0.5)}));
        }

        // check tetrahedron face baycenters
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));

            const vector<Eigen::MatrixXd> faceVertices =
                geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices, tetrahedron.Faces);

            ASSERT_EQ(geometryUtilities.PolyhedronFaceBarycenter(faceVertices),
                      vector<Eigen::Vector3d>({Eigen::Vector3d(1.0 / 3.0, 1.0 / 3.0, 0.0),
                                               Eigen::Vector3d(1.0 / 3.0, 0.0, 1.0 / 3.0),
                                               Eigen::Vector3d(0.0, 1.0 / 3.0, 1.0 / 3.0),
                                               Eigen::Vector3d(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)}));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceTriangulations)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check cube face triangulations
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));
            const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices, cube.Faces);
            const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);

            const vector<vector<unsigned int>> faceTriangulationsByFirstVertex =
                geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(cube.Faces, faceVertices);
            const vector<vector<unsigned int>> faceTriangulationsByInternalPoint =
                geometryUtilities.PolyhedronFaceTriangulationsByInternalPoint(cube.Vertices, cube.Faces, faceVertices, faceBarycenters);

            ASSERT_EQ(faceTriangulationsByFirstVertex,
                      vector<vector<unsigned int>>({vector<unsigned int>({0, 1, 2, 0, 2, 3}),
                                                    vector<unsigned int>({0, 1, 2, 0, 2, 3}),
                                                    vector<unsigned int>({0, 1, 2, 0, 2, 3}),
                                                    vector<unsigned int>({0, 1, 2, 0, 2, 3}),
                                                    vector<unsigned int>({0, 1, 2, 0, 2, 3}),
                                                    vector<unsigned int>({0, 1, 2, 0, 2, 3})}));
            ASSERT_EQ(
                geometryUtilities.PolyhedronFaceExtractTriangulationPoints(faceVertices, faceTriangulationsByFirstVertex),
                vector<vector<Eigen::Matrix3d>>(
                    {vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << cube.Vertices.col(0), cube.Vertices.col(1), cube.Vertices.col(2)).finished(),
                         (Eigen::Matrix3d() << cube.Vertices.col(0), cube.Vertices.col(2), cube.Vertices.col(3)).finished(),
                     }),
                     vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << cube.Vertices.col(4), cube.Vertices.col(5), cube.Vertices.col(6)).finished(),
                         (Eigen::Matrix3d() << cube.Vertices.col(4), cube.Vertices.col(6), cube.Vertices.col(7)).finished(),
                     }),
                     vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << cube.Vertices.col(0), cube.Vertices.col(3), cube.Vertices.col(7)).finished(),
                         (Eigen::Matrix3d() << cube.Vertices.col(0), cube.Vertices.col(7), cube.Vertices.col(4)).finished(),
                     }),
                     vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << cube.Vertices.col(1), cube.Vertices.col(2), cube.Vertices.col(6)).finished(),
                         (Eigen::Matrix3d() << cube.Vertices.col(1), cube.Vertices.col(6), cube.Vertices.col(5)).finished(),
                     }),
                     vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << cube.Vertices.col(0), cube.Vertices.col(1), cube.Vertices.col(5)).finished(),
                         (Eigen::Matrix3d() << cube.Vertices.col(0), cube.Vertices.col(5), cube.Vertices.col(4)).finished(),
                     }),
                     vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << cube.Vertices.col(3), cube.Vertices.col(2), cube.Vertices.col(6)).finished(),
                         (Eigen::Matrix3d() << cube.Vertices.col(3), cube.Vertices.col(6), cube.Vertices.col(7)).finished(),
                     })}));

            ASSERT_EQ(faceTriangulationsByInternalPoint,
                      vector<vector<unsigned int>>({vector<unsigned int>({4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0}),
                                                    vector<unsigned int>({4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0}),
                                                    vector<unsigned int>({4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0}),
                                                    vector<unsigned int>({4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0}),
                                                    vector<unsigned int>({4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0}),
                                                    vector<unsigned int>({4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0})}));

            ASSERT_EQ(
                geometryUtilities.PolyhedronFaceTriangulationPointsByInternalPoint(faceVertices, faceBarycenters, faceTriangulationsByInternalPoint),
                vector<vector<Eigen::Matrix3d>>(
                    {vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << faceBarycenters[0], cube.Vertices.col(0), cube.Vertices.col(1)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[0], cube.Vertices.col(1), cube.Vertices.col(2)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[0], cube.Vertices.col(2), cube.Vertices.col(3)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[0], cube.Vertices.col(3), cube.Vertices.col(0)).finished(),
                     }),
                     vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << faceBarycenters[1], cube.Vertices.col(4), cube.Vertices.col(5)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[1], cube.Vertices.col(5), cube.Vertices.col(6)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[1], cube.Vertices.col(6), cube.Vertices.col(7)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[1], cube.Vertices.col(7), cube.Vertices.col(4)).finished(),
                     }),
                     vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << faceBarycenters[2], cube.Vertices.col(0), cube.Vertices.col(3)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[2], cube.Vertices.col(3), cube.Vertices.col(7)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[2], cube.Vertices.col(7), cube.Vertices.col(4)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[2], cube.Vertices.col(4), cube.Vertices.col(0)).finished(),
                     }),
                     vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << faceBarycenters[3], cube.Vertices.col(1), cube.Vertices.col(2)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[3], cube.Vertices.col(2), cube.Vertices.col(6)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[3], cube.Vertices.col(6), cube.Vertices.col(5)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[3], cube.Vertices.col(5), cube.Vertices.col(1)).finished(),
                     }),
                     vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << faceBarycenters[4], cube.Vertices.col(0), cube.Vertices.col(1)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[4], cube.Vertices.col(1), cube.Vertices.col(5)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[4], cube.Vertices.col(5), cube.Vertices.col(4)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[4], cube.Vertices.col(4), cube.Vertices.col(0)).finished(),
                     }),
                     vector<Eigen::Matrix3d>({
                         (Eigen::Matrix3d() << faceBarycenters[5], cube.Vertices.col(3), cube.Vertices.col(2)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[5], cube.Vertices.col(2), cube.Vertices.col(6)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[5], cube.Vertices.col(6), cube.Vertices.col(7)).finished(),
                         (Eigen::Matrix3d() << faceBarycenters[5], cube.Vertices.col(7), cube.Vertices.col(3)).finished(),
                     })}));
        }

        // check tetrahedron face triangulations
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));

            const vector<Eigen::MatrixXd> faceVertices =
                geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices, tetrahedron.Faces);
            const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);

            ASSERT_EQ(geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(tetrahedron.Faces, faceVertices),
                      vector<vector<unsigned int>>({vector<unsigned int>({0, 1, 2}),
                                                    vector<unsigned int>({0, 1, 2}),
                                                    vector<unsigned int>({0, 1, 2}),
                                                    vector<unsigned int>({0, 1, 2})}));
            ASSERT_EQ(geometryUtilities.PolyhedronFaceTriangulationsByInternalPoint(tetrahedron.Vertices, tetrahedron.Faces, faceVertices, faceBarycenters),
                      vector<vector<unsigned int>>({vector<unsigned int>({3, 0, 1, 3, 1, 2, 3, 2, 0}),
                                                    vector<unsigned int>({3, 0, 1, 3, 1, 2, 3, 2, 0}),
                                                    vector<unsigned int>({3, 0, 1, 3, 1, 2, 3, 2, 0}),
                                                    vector<unsigned int>({3, 0, 1, 3, 1, 2, 3, 2, 0})}));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronTetrahedrons)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        std::string exportFolder = "./Export/TestPolyhedronTetrahedrons";
        Gedim::Output::CreateFolder(exportFolder);

        // check cube face triangulations
        {
            const Gedim::GeometryUtilities::Polyhedron cube =
                geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                 Eigen::Vector3d(0.0, 1.0, 0.0));
            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(cube.Vertices);
            const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices, cube.Faces);
            const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);
            const vector<vector<unsigned int>> faceTriangulations =
                geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(cube.Faces, faceVertices);
            const vector<vector<unsigned int>> faceTriangulationsByInternalPoint =
                geometryUtilities.PolyhedronFaceTriangulationsByInternalPoint(cube.Vertices, cube.Faces, faceVertices, faceBarycenters);
            const vector<unsigned int> tetrahedronList =
                geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(cube.Vertices, cube.Faces, faceTriangulations, polyhedronBarycenter);
            const vector<unsigned int> tetrahedronByInternalPointsList =
                geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(cube.Vertices,
                                                                             cube.Faces,
                                                                             faceTriangulationsByInternalPoint,
                                                                             faceBarycenters,
                                                                             polyhedronBarycenter);
            ASSERT_EQ(tetrahedronList,
                      vector<unsigned int>({0, 1, 2, 8, 0, 2, 3, 8, 4, 5, 6, 8, 4, 6, 7, 8, 0, 3, 7, 8, 0, 7, 4, 8,
                                            1, 2, 6, 8, 1, 6, 5, 8, 0, 1, 5, 8, 0, 5, 4, 8, 3, 2, 6, 8, 3, 6, 7, 8}));
            ASSERT_EQ(tetrahedronByInternalPointsList,
                      vector<unsigned int>({8,  0, 1, 14, 8,  1, 2, 14, 8,  2, 3, 14, 8,  3, 0, 14, 9,  4, 5, 14,
                                            9,  5, 6, 14, 9,  6, 7, 14, 9,  7, 4, 14, 10, 0, 3, 14, 10, 3, 7, 14,
                                            10, 7, 4, 14, 10, 4, 0, 14, 11, 1, 2, 14, 11, 2, 6, 14, 11, 6, 5, 14,
                                            11, 5, 1, 14, 12, 0, 1, 14, 12, 1, 5, 14, 12, 5, 4, 14, 12, 4, 0, 14,
                                            13, 3, 2, 14, 13, 2, 6, 14, 13, 6, 7, 14, 13, 7, 3, 14}));
            // Export tetrahedrons
            {
                vector<Eigen::MatrixXd> tetrahedrons =
                    geometryUtilities.ExtractTetrahedronPoints(cube.Vertices, polyhedronBarycenter, tetrahedronList);

                Gedim::VTKUtilities vtkExperter;
                for (unsigned int t = 0; t < tetrahedrons.size(); t++)
                {
                    Gedim::GeometryUtilities::Polyhedron subTetra =
                        geometryUtilities.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
                                                                        tetrahedrons[t].col(1),
                                                                        tetrahedrons[t].col(2),
                                                                        tetrahedrons[t].col(3));
                    vector<double> id(1, t);

                    vtkExperter.AddPolyhedron(
                        subTetra.Vertices,
                        subTetra.Edges,
                        subTetra.Faces,
                        {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()}});
                }

                vtkExperter.Export(exportFolder + "/Cube_Tetra_1.vtu", Gedim::VTKUtilities::Ascii);
            }

            {
                vector<Eigen::MatrixXd> tetrahedrons =
                    geometryUtilities.ExtractTetrahedronPoints(cube.Vertices, polyhedronBarycenter, faceBarycenters, tetrahedronByInternalPointsList);

                Gedim::VTKUtilities vtkExperter;
                for (unsigned int t = 0; t < tetrahedrons.size(); t++)
                {
                    Gedim::GeometryUtilities::Polyhedron subTetra =
                        geometryUtilities.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
                                                                        tetrahedrons[t].col(1),
                                                                        tetrahedrons[t].col(2),
                                                                        tetrahedrons[t].col(3));
                    vector<double> id(1, t);

                    vtkExperter.AddPolyhedron(
                        subTetra.Vertices,
                        subTetra.Edges,
                        subTetra.Faces,
                        {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()}});
                }

                vtkExperter.Export(exportFolder + "/Cube_Tetra_2.vtu", Gedim::VTKUtilities::Ascii);
            }
        }

        // check tetrahedron face triangulations
        {
            const Gedim::GeometryUtilities::Polyhedron tetrahedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));

            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(tetrahedron.Vertices);
            const vector<Eigen::MatrixXd> faceVertices =
                geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices, tetrahedron.Faces);
            const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);
            const vector<vector<unsigned int>> faceTriangulations =
                geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(tetrahedron.Faces, faceVertices);
            const vector<vector<unsigned int>> faceTriangulationsByInternalPoint =
                geometryUtilities.PolyhedronFaceTriangulationsByInternalPoint(tetrahedron.Vertices, tetrahedron.Faces, faceVertices, faceBarycenters);

            const vector<unsigned int> tetrahedronList =
                geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(tetrahedron.Vertices,
                                                                             tetrahedron.Faces,
                                                                             faceTriangulations,
                                                                             polyhedronBarycenter);
            const vector<unsigned int> tetrahedronByInternalPointsList =
                geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(tetrahedron.Vertices,
                                                                             tetrahedron.Faces,
                                                                             faceTriangulationsByInternalPoint,
                                                                             faceBarycenters,
                                                                             polyhedronBarycenter);
            ASSERT_EQ(tetrahedronList, vector<unsigned int>({0, 1, 2, 4, 0, 1, 3, 4, 0, 2, 3, 4, 1, 2, 3, 4}));
            ASSERT_EQ(tetrahedronByInternalPointsList,
                      vector<unsigned int>({4, 0, 1, 8, 4, 1, 2, 8, 4, 2, 0, 8, 5, 0, 1, 8, 5, 1, 3, 8, 5, 3, 0, 8,
                                            6, 0, 2, 8, 6, 2, 3, 8, 6, 3, 0, 8, 7, 1, 2, 8, 7, 2, 3, 8, 7, 3, 1, 8}));

            // Export tetrahedrons
            {
                vector<Eigen::MatrixXd> tetrahedrons =
                    geometryUtilities.ExtractTetrahedronPoints(tetrahedron.Vertices, polyhedronBarycenter, tetrahedronList);

                Gedim::VTKUtilities vtkExperter;
                for (unsigned int t = 0; t < tetrahedrons.size(); t++)
                {
                    Gedim::GeometryUtilities::Polyhedron subTetra =
                        geometryUtilities.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
                                                                        tetrahedrons[t].col(1),
                                                                        tetrahedrons[t].col(2),
                                                                        tetrahedrons[t].col(3));
                    vector<double> id(1, t);

                    vtkExperter.AddPolyhedron(
                        subTetra.Vertices,
                        subTetra.Edges,
                        subTetra.Faces,
                        {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()}});
                }

                vtkExperter.Export(exportFolder + "/Tetrahedrons_Tetra_1.vtu", Gedim::VTKUtilities::Ascii);
            }

            {
                vector<Eigen::MatrixXd> tetrahedrons =
                    geometryUtilities.ExtractTetrahedronPoints(tetrahedron.Vertices, polyhedronBarycenter, faceBarycenters, tetrahedronByInternalPointsList);

                Gedim::VTKUtilities vtkExperter;
                for (unsigned int t = 0; t < tetrahedrons.size(); t++)
                {
                    Gedim::GeometryUtilities::Polyhedron subTetra =
                        geometryUtilities.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
                                                                        tetrahedrons[t].col(1),
                                                                        tetrahedrons[t].col(2),
                                                                        tetrahedrons[t].col(3));
                    vector<double> id(1, t);

                    vtkExperter.AddPolyhedron(
                        subTetra.Vertices,
                        subTetra.Edges,
                        subTetra.Faces,
                        {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()}});
                }

                vtkExperter.Export(exportFolder + "/Tetrahedrons_Tetra_2.vtu", Gedim::VTKUtilities::Ascii);
            }
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronVolume)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-15;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    // check cube volume
    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron =
            geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::MatrixXd> polyhedronFace3DVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const vector<vector<unsigned int>> polyhedronFaceTriangulations =
            geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces, polyhedronFace3DVertices);

        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
        const vector<bool> polyhedronFaceNormalDirections =
            geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices, polyhedronBarycenter, polyhedronFaceNormals);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations =
            geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
            geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

        const vector<Eigen::MatrixXd> polyhedronFace2DVertices =
            geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices, polyhedronFaceTranslations, polyhedronFaceRotationMatrices);

        const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints =
            geometryUtilities.PolyhedronFaceExtractTriangulationPoints(polyhedronFace2DVertices, polyhedronFaceTriangulations);

        const double polyhedronVolume = geometryUtilities.PolyhedronVolumeByBoundaryIntegral(polyhedronFace2DTriangulationPoints,
                                                                                             polyhedronFaceNormals,
                                                                                             polyhedronFaceNormalDirections,
                                                                                             polyhedronFaceTranslations,
                                                                                             polyhedronFaceRotationMatrices);
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(1.0, polyhedronVolume, geometryUtilities.Tolerance3D()));

        const vector<unsigned int> polyhedronTetrahedrons =
            geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(polyhedron.Vertices,
                                                                         polyhedron.Faces,
                                                                         polyhedronFaceTriangulations,
                                                                         polyhedronBarycenter);
        const vector<Eigen::MatrixXd> polyhedronTetrahedronPoints =
            geometryUtilities.ExtractTetrahedronPoints(polyhedron.Vertices, polyhedronBarycenter, polyhedronTetrahedrons);
        const double polyhedronVolumeByInternal = geometryUtilities.PolyhedronVolumeByInternalIntegral(polyhedronTetrahedronPoints);

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(1.0, polyhedronVolumeByInternal, geometryUtilities.Tolerance3D()));
    }

    // check tetrahedron volume
    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron =
            geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                          Eigen::Vector3d(1.0, 0.0, 0.0),
                                                          Eigen::Vector3d(0.0, 0.0, 1.0),
                                                          Eigen::Vector3d(0.0, 1.0, 0.0));

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(1.0 / 6.0,
                                                     geometryUtilities.SimplexMeasure(polyhedron.Vertices),
                                                     geometryUtilities.Tolerance3D()));

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::MatrixXd> polyhedronFace3DVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const vector<vector<unsigned int>> polyhedronFaceTriangulations =
            geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces, polyhedronFace3DVertices);

        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
        const vector<bool> polyhedronFaceNormalDirections =
            geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices, polyhedronBarycenter, polyhedronFaceNormals);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations =
            geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
            geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

        const vector<Eigen::MatrixXd> polyhedronFace2DVertices =
            geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices, polyhedronFaceTranslations, polyhedronFaceRotationMatrices);

        const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints =
            geometryUtilities.PolyhedronFaceExtractTriangulationPoints(polyhedronFace2DVertices, polyhedronFaceTriangulations);

        const double polyhedronVolume = geometryUtilities.PolyhedronVolumeByBoundaryIntegral(polyhedronFace2DTriangulationPoints,
                                                                                             polyhedronFaceNormals,
                                                                                             polyhedronFaceNormalDirections,
                                                                                             polyhedronFaceTranslations,
                                                                                             polyhedronFaceRotationMatrices);

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(1.0 / 6.0, polyhedronVolume, geometryUtilities.Tolerance3D()));

        const vector<unsigned int> polyhedronTetrahedrons =
            geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(polyhedron.Vertices,
                                                                         polyhedron.Faces,
                                                                         polyhedronFaceTriangulations,
                                                                         polyhedronBarycenter);
        const vector<Eigen::MatrixXd> polyhedronTetrahedronPoints =
            geometryUtilities.ExtractTetrahedronPoints(polyhedron.Vertices, polyhedronBarycenter, polyhedronTetrahedrons);
        const double polyhedronVolumeByInternal = geometryUtilities.PolyhedronVolumeByInternalIntegral(polyhedronTetrahedronPoints);

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(1.0 / 6.0, polyhedronVolumeByInternal, geometryUtilities.Tolerance3D()));
    }
}

TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronCentroid)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-14;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    // check cube centroid
    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron =
            geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::MatrixXd> polyhedronFace3DVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const vector<vector<unsigned int>> polyhedronFaceTriangulations =
            geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces, polyhedronFace3DVertices);

        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
        const vector<bool> polyhedronFaceNormalDirections =
            geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices, polyhedronBarycenter, polyhedronFaceNormals);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations =
            geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
            geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

        const vector<Eigen::MatrixXd> polyhedronFace2DVertices =
            geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices, polyhedronFaceTranslations, polyhedronFaceRotationMatrices);

        const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints =
            geometryUtilities.PolyhedronFaceExtractTriangulationPoints(polyhedronFace2DVertices, polyhedronFaceTriangulations);

        const double polyhedronVolume = geometryUtilities.PolyhedronVolumeByBoundaryIntegral(polyhedronFace2DTriangulationPoints,
                                                                                             polyhedronFaceNormals,
                                                                                             polyhedronFaceNormalDirections,
                                                                                             polyhedronFaceTranslations,
                                                                                             polyhedronFaceRotationMatrices);

        const Eigen::Vector3d polyhedronCentroid = geometryUtilities.PolyhedronCentroid(polyhedronFace2DTriangulationPoints,
                                                                                        polyhedronFaceNormals,
                                                                                        polyhedronFaceNormalDirections,
                                                                                        polyhedronFaceTranslations,
                                                                                        polyhedronFaceRotationMatrices,
                                                                                        polyhedronVolume);

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronCentroid.x(), 0.5, geometryUtilities.Tolerance1D()));
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronCentroid.y(), 0.5, geometryUtilities.Tolerance1D()));
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronCentroid.z(), 0.5, geometryUtilities.Tolerance1D()));
    }

    // check tetrahedron volume
    {
        const Gedim::GeometryUtilities::Polyhedron polyhedron =
            geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                          Eigen::Vector3d(1.0, 0.0, 0.0),
                                                          Eigen::Vector3d(0.0, 0.0, 1.0),
                                                          Eigen::Vector3d(0.0, 1.0, 0.0));

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::MatrixXd> polyhedronFace3DVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const vector<vector<unsigned int>> polyhedronFaceTriangulations =
            geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces, polyhedronFace3DVertices);

        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
        const vector<bool> polyhedronFaceNormalDirections =
            geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices, polyhedronBarycenter, polyhedronFaceNormals);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations =
            geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
            geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

        const vector<Eigen::MatrixXd> polyhedronFace2DVertices =
            geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices, polyhedronFaceTranslations, polyhedronFaceRotationMatrices);

        const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints =
            geometryUtilities.PolyhedronFaceExtractTriangulationPoints(polyhedronFace2DVertices, polyhedronFaceTriangulations);

        const double polyhedronVolume = geometryUtilities.PolyhedronVolumeByBoundaryIntegral(polyhedronFace2DTriangulationPoints,
                                                                                             polyhedronFaceNormals,
                                                                                             polyhedronFaceNormalDirections,
                                                                                             polyhedronFaceTranslations,
                                                                                             polyhedronFaceRotationMatrices);

        const Eigen::Vector3d polyhedronCentroid = geometryUtilities.PolyhedronCentroid(polyhedronFace2DTriangulationPoints,
                                                                                        polyhedronFaceNormals,
                                                                                        polyhedronFaceNormalDirections,
                                                                                        polyhedronFaceTranslations,
                                                                                        polyhedronFaceRotationMatrices,
                                                                                        polyhedronVolume);

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronCentroid.x(), 1.0 / 4.0, geometryUtilities.Tolerance1D()));
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronCentroid.y(), 1.0 / 4.0, geometryUtilities.Tolerance1D()));
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronCentroid.z(), 1.0 / 4.0, geometryUtilities.Tolerance1D()));
    }
}

TEST(TestGeometryUtilities, TestPolyhedronIsConvex_Tetrahedron)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        std::string exportFolder = "./Export/TestPolyhedronIsConvex_Tetrahedron";
        Gedim::Output::CreateFolder(exportFolder);

        // check convex polyhedron
        {
            const Gedim::GeometryUtilities::Polyhedron polyhedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));

            geometryUtilities.ExportPolyhedronToVTU(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, exportFolder);

            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
            const vector<Eigen::MatrixXd> polyhedronFace3DVertices =
                geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
            const vector<Eigen::Vector3d> polyhedronFaceBarycenters =
                geometryUtilities.PolyhedronFaceBarycenter(polyhedronFace3DVertices);

            const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
            const vector<bool> polyhedronFaceNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices, polyhedronBarycenter, polyhedronFaceNormals);
            const vector<Eigen::Vector3d> polyhedronFaceTranslations =
                geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
            const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

            const vector<Eigen::MatrixXd> polyhedronFace2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices, polyhedronFaceTranslations, polyhedronFaceRotationMatrices);

            ASSERT_TRUE(geometryUtilities.PolyhedronIsConvex(polyhedronFace3DVertices,
                                                             polyhedronFace2DVertices,
                                                             polyhedronFaceBarycenters,
                                                             polyhedronFaceNormals,
                                                             polyhedronFaceNormalDirections,
                                                             polyhedronBarycenter));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedronIsConvex_Cube)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        std::string exportFolder = "./Export/TestPolyhedronIsConvex_Cube";
        Gedim::Output::CreateFolder(exportFolder);

        // check convex polyhedron
        {
            const Gedim::GeometryUtilities::Polyhedron polyhedron =
                geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);

            geometryUtilities.ExportPolyhedronToVTU(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, exportFolder);

            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
            const vector<Eigen::MatrixXd> polyhedronFace3DVertices =
                geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
            const vector<Eigen::Vector3d> polyhedronFaceBarycenters =
                geometryUtilities.PolyhedronFaceBarycenter(polyhedronFace3DVertices);

            const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
            const vector<bool> polyhedronFaceNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices, polyhedronBarycenter, polyhedronFaceNormals);
            const vector<Eigen::Vector3d> polyhedronFaceTranslations =
                geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
            const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

            const vector<Eigen::MatrixXd> polyhedronFace2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices, polyhedronFaceTranslations, polyhedronFaceRotationMatrices);

            ASSERT_TRUE(geometryUtilities.PolyhedronIsConvex(polyhedronFace3DVertices,
                                                             polyhedronFace2DVertices,
                                                             polyhedronFaceBarycenters,
                                                             polyhedronFaceNormals,
                                                             polyhedronFaceNormalDirections,
                                                             polyhedronBarycenter));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedronIsConvex_Parallellepiped)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        std::string exportFolder = "./Export/TestPolyhedronIsConvex_Parallellepiped";
        Gedim::Output::CreateFolder(exportFolder);

        // check convex polyhedron
        {
            Eigen::MatrixXd polygon(3, 4);
            polygon.col(0) << 0.0, 0.0, 0.0;
            polygon.col(1) << 1.0, 0.25, 0.0;
            polygon.col(2) << 0.75, 0.75, 0.0;
            polygon.col(3) << 0.25, 1.0, 0.0;

            const Gedim::GeometryUtilities::Polyhedron polyhedron =
                geometryUtilities.CreatePolyhedronWithExtrusion(polygon, Eigen::Vector3d(0.5, 0.25, 1.0));

            geometryUtilities.ExportPolyhedronToVTU(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, exportFolder);

            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
            const vector<Eigen::MatrixXd> polyhedronFace3DVertices =
                geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
            const vector<Eigen::Vector3d> polyhedronFaceBarycenters =
                geometryUtilities.PolyhedronFaceBarycenter(polyhedronFace3DVertices);

            const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
            const vector<bool> polyhedronFaceNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices, polyhedronBarycenter, polyhedronFaceNormals);
            const vector<Eigen::Vector3d> polyhedronFaceTranslations =
                geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
            const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

            const vector<Eigen::MatrixXd> polyhedronFace2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices, polyhedronFaceTranslations, polyhedronFaceRotationMatrices);

            ASSERT_TRUE(geometryUtilities.PolyhedronIsConvex(polyhedronFace3DVertices,
                                                             polyhedronFace2DVertices,
                                                             polyhedronFaceBarycenters,
                                                             polyhedronFaceNormals,
                                                             polyhedronFaceNormalDirections,
                                                             polyhedronBarycenter));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedronIsConvex_Concave)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        std::string exportFolder = "./Export/TestPolyhedronIsConvex_Concave";
        Gedim::Output::CreateFolder(exportFolder);

        // check concave polyhedron
        {
            Gedim::GeometryUtilities::Polyhedron polyhedron;

            // create vertices
            polyhedron.Vertices.resize(3, 5);
            polyhedron.Vertices.col(0) << 0.0, 0.0, 0.0;
            polyhedron.Vertices.col(1) << 1.0, 0.0, 0.0;
            polyhedron.Vertices.col(2) << 0.0, 1.0, 0.0;
            polyhedron.Vertices.col(3) << 0.0, 0.0, 1.0;
            polyhedron.Vertices.col(4) << 0.25, 0.25, 0.25;

            // create edges
            polyhedron.Edges.resize(2, 9);
            polyhedron.Edges.col(0) << 0, 1;
            polyhedron.Edges.col(1) << 0, 2;
            polyhedron.Edges.col(2) << 1, 2;
            polyhedron.Edges.col(3) << 0, 3;
            polyhedron.Edges.col(4) << 1, 3;
            polyhedron.Edges.col(5) << 2, 3;
            polyhedron.Edges.col(6) << 1, 4;
            polyhedron.Edges.col(7) << 2, 4;
            polyhedron.Edges.col(8) << 3, 4;

            // create faces
            polyhedron.Faces.resize(6, Eigen::MatrixXi::Zero(2, 3));
            polyhedron.Faces[0].row(0) << 0, 1, 2;
            polyhedron.Faces[1].row(0) << 0, 1, 3;
            polyhedron.Faces[2].row(0) << 0, 2, 3;
            polyhedron.Faces[3].row(0) << 1, 2, 4;
            polyhedron.Faces[4].row(0) << 2, 3, 4;
            polyhedron.Faces[5].row(0) << 3, 1, 4;

            polyhedron.Faces[0].row(1) << 0, 2, 1;
            polyhedron.Faces[1].row(1) << 0, 4, 3;
            polyhedron.Faces[2].row(1) << 1, 5, 3;
            polyhedron.Faces[3].row(1) << 2, 7, 6;
            polyhedron.Faces[4].row(1) << 5, 8, 7;
            polyhedron.Faces[5].row(1) << 4, 6, 8;

            geometryUtilities.ExportPolyhedronToVTU(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, exportFolder);

            const Eigen::Vector3d polyhedronPointInside(0.125, 0.125, 0.125);
            const vector<Eigen::MatrixXd> polyhedronFace3DVertices =
                geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
            const vector<Eigen::Vector3d> polyhedronFaceBarycenters =
                geometryUtilities.PolyhedronFaceBarycenter(polyhedronFace3DVertices);

            const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
            const vector<bool> polyhedronFaceNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices, polyhedronPointInside, polyhedronFaceNormals);
            const vector<Eigen::Vector3d> polyhedronFaceTranslations =
                geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
            const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

            const vector<Eigen::MatrixXd> polyhedronFace2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices, polyhedronFaceTranslations, polyhedronFaceRotationMatrices);

            ASSERT_FALSE(geometryUtilities.PolyhedronIsConvex(polyhedronFace3DVertices,
                                                              polyhedronFace2DVertices,
                                                              polyhedronFaceBarycenters,
                                                              polyhedronFaceNormals,
                                                              polyhedronFaceNormalDirections,
                                                              polyhedronPointInside));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestUnalignedPolyhedronPoints)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    std::string exportFolder = "./Export/TestUnalignedPolyhedronPoints";
    Gedim::Output::CreateFolder(exportFolder);

    {
        // Aligned thetrahedron
        Gedim::GeometryUtilities::Polyhedron polyhedron;

        // create vertices
        polyhedron.Vertices.resize(3, 11);
        polyhedron.Vertices.col(0) << 0.0, 0.0, 0.0;
        polyhedron.Vertices.col(1) << 0.1, 0.0, 0.0;
        polyhedron.Vertices.col(2) << 1.0, 0.0, 0.0;
        polyhedron.Vertices.col(3) << 0.75, 0.0, 0.25;
        polyhedron.Vertices.col(4) << 0.25, 0.0, 0.75;
        polyhedron.Vertices.col(5) << 0.0, 0.0, 1.0;
        polyhedron.Vertices.col(6) << 0.0, 0.0, 0.5;
        polyhedron.Vertices.col(7) << 0.0, 0.5, 0.0;
        polyhedron.Vertices.col(8) << 0.5, 0.5, 0.0;
        polyhedron.Vertices.col(9) << 0.0, 1.0, 0.0;
        polyhedron.Vertices.col(10) << 0.0, 0.5, 0.5;

        // create edges
        polyhedron.Edges.resize(2, 14);
        polyhedron.Edges.col(0) << 0, 1;
        polyhedron.Edges.col(1) << 1, 2;
        polyhedron.Edges.col(2) << 2, 3;
        polyhedron.Edges.col(3) << 3, 4;
        polyhedron.Edges.col(4) << 4, 5;
        polyhedron.Edges.col(5) << 5, 6;
        polyhedron.Edges.col(6) << 6, 0;
        polyhedron.Edges.col(7) << 0, 7;
        polyhedron.Edges.col(8) << 7, 9;
        polyhedron.Edges.col(9) << 9, 8;
        polyhedron.Edges.col(10) << 2, 8;
        polyhedron.Edges.col(11) << 5, 10;
        polyhedron.Edges.col(12) << 9, 10;
        polyhedron.Edges.col(13) << 7, 10;

        // create faces
        polyhedron.Faces.resize(5);

        polyhedron.Faces[0].resize(2, 7);
        polyhedron.Faces[0].row(0) << 0, 1, 2, 3, 4, 5, 6;
        polyhedron.Faces[0].row(1) << 0, 1, 2, 3, 4, 5, 6;

        polyhedron.Faces[1].resize(2, 7);
        polyhedron.Faces[1].row(0) << 2, 8, 9, 10, 5, 4, 3;
        polyhedron.Faces[1].row(1) << 10, 9, 12, 11, 4, 3, 2;

        polyhedron.Faces[2].resize(2, 6);
        polyhedron.Faces[2].row(0) << 0, 1, 2, 8, 9, 7;
        polyhedron.Faces[2].row(1) << 0, 1, 10, 9, 8, 7;

        polyhedron.Faces[3].resize(2, 5);
        polyhedron.Faces[3].row(0) << 0, 7, 10, 5, 6;
        polyhedron.Faces[3].row(1) << 7, 13, 11, 5, 6;

        polyhedron.Faces[4].resize(2, 3);
        polyhedron.Faces[4].row(0) << 7, 9, 10;
        polyhedron.Faces[4].row(1) << 8, 12, 13;

        geometryUtilities.ExportPolyhedronToVTU(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, exportFolder);

        const vector<Eigen::MatrixXd> faceVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
        const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);
        const vector<Eigen::Matrix3d> faceRotationMatrices =
            geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices, faceNormals, faceTranslations);

        const vector<Eigen::MatrixXd> face2DVertices =
            geometryUtilities.PolyhedronFaceRotatedVertices(faceVertices, faceTranslations, faceRotationMatrices);
        const std::vector<std::vector<unsigned int>> facesUnalignedPoints =
            geometryUtilities.PolyhedronFacesUnalignedVertices(face2DVertices);

        const std::vector<std::vector<unsigned int>> polyhedronUnaligedFaces = {{0}, {1}, {2}, {3, 4}};

        const std::vector<unsigned int> unalignedVertices =
            geometryUtilities.UnalignedPolyhedronPoints(polyhedron.Vertices,
                                                        polyhedron.Faces,
                                                        faceTranslations,
                                                        faceRotationMatrices,
                                                        polyhedronUnaligedFaces,
                                                        facesUnalignedPoints);

        ASSERT_EQ(std::vector<unsigned int>({0, 2, 5, 9}), unalignedVertices);
    }
}

TEST(TestGeometryUtilities, TestAlignedPolyhedronEdges)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::GraphUtilities graphUtilities;

    std::string exportFolder = "./Export/TestAlignedPolyhedronEdges";
    Gedim::Output::CreateFolder(exportFolder);

    {
        // Aligned thetrahedron
        Gedim::GeometryUtilities::Polyhedron polyhedron;

        // create vertices
        polyhedron.Vertices.resize(3, 11);
        polyhedron.Vertices.col(0) << 0.0, 0.0, 0.0;
        polyhedron.Vertices.col(1) << 0.1, 0.0, 0.0;
        polyhedron.Vertices.col(2) << 1.0, 0.0, 0.0;
        polyhedron.Vertices.col(3) << 0.75, 0.0, 0.25;
        polyhedron.Vertices.col(4) << 0.25, 0.0, 0.75;
        polyhedron.Vertices.col(5) << 0.0, 0.0, 1.0;
        polyhedron.Vertices.col(6) << 0.0, 0.0, 0.5;
        polyhedron.Vertices.col(7) << 0.0, 0.5, 0.0;
        polyhedron.Vertices.col(8) << 0.5, 0.5, 0.0;
        polyhedron.Vertices.col(9) << 0.0, 1.0, 0.0;
        polyhedron.Vertices.col(10) << 0.0, 0.5, 0.5;

        // create edges
        polyhedron.Edges.resize(2, 14);
        polyhedron.Edges.col(0) << 0, 1;
        polyhedron.Edges.col(1) << 1, 2;
        polyhedron.Edges.col(2) << 2, 3;
        polyhedron.Edges.col(3) << 3, 4;
        polyhedron.Edges.col(4) << 4, 5;
        polyhedron.Edges.col(5) << 5, 6;
        polyhedron.Edges.col(6) << 6, 0;
        polyhedron.Edges.col(7) << 0, 7;
        polyhedron.Edges.col(8) << 7, 9;
        polyhedron.Edges.col(9) << 9, 8;
        polyhedron.Edges.col(10) << 2, 8;
        polyhedron.Edges.col(11) << 5, 10;
        polyhedron.Edges.col(12) << 9, 10;
        polyhedron.Edges.col(13) << 7, 10;

        // create faces
        polyhedron.Faces.resize(5);

        polyhedron.Faces[0].resize(2, 7);
        polyhedron.Faces[0].row(0) << 0, 1, 2, 3, 4, 5, 6;
        polyhedron.Faces[0].row(1) << 0, 1, 2, 3, 4, 5, 6;

        polyhedron.Faces[1].resize(2, 7);
        polyhedron.Faces[1].row(0) << 2, 8, 9, 10, 5, 4, 3;
        polyhedron.Faces[1].row(1) << 10, 9, 12, 11, 4, 3, 2;

        polyhedron.Faces[2].resize(2, 6);
        polyhedron.Faces[2].row(0) << 0, 1, 2, 8, 9, 7;
        polyhedron.Faces[2].row(1) << 0, 1, 10, 9, 8, 7;

        polyhedron.Faces[3].resize(2, 5);
        polyhedron.Faces[3].row(0) << 0, 7, 10, 5, 6;
        polyhedron.Faces[3].row(1) << 7, 13, 11, 5, 6;

        polyhedron.Faces[4].resize(2, 3);
        polyhedron.Faces[4].row(0) << 7, 9, 10;
        polyhedron.Faces[4].row(1) << 8, 12, 13;

        geometryUtilities.ExportPolyhedronToVTU(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, exportFolder);

        const Eigen::MatrixXd edgesTangent = geometryUtilities.PolyhedronEdgeTangents(polyhedron.Vertices, polyhedron.Edges);

        const Eigen::VectorXd edgesLength = geometryUtilities.PolyhedronEdgesLength(polyhedron.Vertices, polyhedron.Edges);

        const Gedim::GraphUtilities::GraphAdjacencyData edgesAdjacency =
            graphUtilities.GraphConnectivityToGraphAdjacency(polyhedron.Vertices.cols(), polyhedron.Edges, false);

        const Gedim::GeometryUtilities::AlignedPolyhedronEdgesResult alignedEdges =
            geometryUtilities.AlignedPolyhedronEdges(polyhedron.Vertices,
                                                     edgesAdjacency.GraphAdjacencyVertices,
                                                     edgesAdjacency.GraphAdjacencyEdges,
                                                     edgesAdjacency.GraphAdjacencyVerticesMap,
                                                     edgesTangent,
                                                     edgesLength.array().square());

        ASSERT_EQ(7, alignedEdges.AlignedEdgesVertices.size());
        ASSERT_EQ(7, alignedEdges.AlignedEdgesEdges.size());
        ASSERT_EQ(std::vector<unsigned int>({0, 1, 2}), alignedEdges.AlignedEdgesVertices[0]);
        ASSERT_EQ(std::vector<unsigned int>({0, 1}), alignedEdges.AlignedEdgesEdges[0]);
        ASSERT_EQ(std::vector<unsigned int>({5, 6, 0}), alignedEdges.AlignedEdgesVertices[1]);
        ASSERT_EQ(std::vector<unsigned int>({5, 6}), alignedEdges.AlignedEdgesEdges[1]);
        ASSERT_EQ(std::vector<unsigned int>({0, 7, 9}), alignedEdges.AlignedEdgesVertices[2]);
        ASSERT_EQ(std::vector<unsigned int>({7, 8}), alignedEdges.AlignedEdgesEdges[2]);
        ASSERT_EQ(std::vector<unsigned int>({7, 10}), alignedEdges.AlignedEdgesVertices[3]);
        ASSERT_EQ(std::vector<unsigned int>({13}), alignedEdges.AlignedEdgesEdges[3]);
        ASSERT_EQ(std::vector<unsigned int>({2, 3, 4, 5}), alignedEdges.AlignedEdgesVertices[4]);
        ASSERT_EQ(std::vector<unsigned int>({2, 3, 4}), alignedEdges.AlignedEdgesEdges[4]);
        ASSERT_EQ(std::vector<unsigned int>({2, 8, 9}), alignedEdges.AlignedEdgesVertices[5]);
        ASSERT_EQ(std::vector<unsigned int>({10, 9}), alignedEdges.AlignedEdgesEdges[5]);
        ASSERT_EQ(std::vector<unsigned int>({5, 10, 9}), alignedEdges.AlignedEdgesVertices[6]);
        ASSERT_EQ(std::vector<unsigned int>({11, 12}), alignedEdges.AlignedEdgesEdges[6]);
    }
}

TEST(TestGeometryUtilities, TestUnalignedPolyhedronPoints_Two)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    std::string exportFolder = "./Export/TestUnalignedPolyhedronPoints_Two";
    Gedim::Output::CreateFolder(exportFolder);

    {
        // Aligned thetrahedron
        Gedim::GeometryUtilities::Polyhedron polyhedron;

        // create vertices
        polyhedron.Vertices.resize(3, 6);
        polyhedron.Vertices.row(0) << 0.0000000000000000e+00, 0.0000000000000000e+00, 5.0000000000000000e-01,
            0.0000000000000000e+00, 2.5000000000000000e-01, 1.6666666666666666e-01;
        polyhedron.Vertices.row(1) << 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, -5.0000000000000000e-01, -3.3333333333333331e-01;
        polyhedron.Vertices.row(2) << 1.0000000000000000e+00, 0.0000000000000000e+00, 5.0000000000000000e-01,
            5.0000000000000000e-01, 7.5000000000000000e-01, 5.0000000000000000e-01;

        // create edges
        polyhedron.Edges.resize(2, 9);
        polyhedron.Edges.row(0) << 2, 1, 4, 2, 3, 1, 0, 5, 5;
        polyhedron.Edges.row(1) << 4, 2, 5, 0, 0, 3, 4, 1, 0;

        // create faces
        polyhedron.Faces.resize(5);

        polyhedron.Faces[0].resize(2, 4);
        polyhedron.Faces[0].row(0) << 1, 2, 0, 3;
        polyhedron.Faces[0].row(1) << 1, 3, 4, 5;

        polyhedron.Faces[1].resize(2, 3);
        polyhedron.Faces[1].row(0) << 0, 4, 2;
        polyhedron.Faces[1].row(1) << 6, 0, 3;

        polyhedron.Faces[2].resize(2, 4);
        polyhedron.Faces[2].row(0) << 3, 1, 5, 0;
        polyhedron.Faces[2].row(1) << 5, 7, 8, 4;

        polyhedron.Faces[3].resize(2, 4);
        polyhedron.Faces[3].row(0) << 4, 2, 1, 5;
        polyhedron.Faces[3].row(1) << 0, 1, 7, 2;

        polyhedron.Faces[4].resize(2, 3);
        polyhedron.Faces[4].row(0) << 5, 4, 0;
        polyhedron.Faces[4].row(1) << 2, 6, 8;

        geometryUtilities.ExportPolyhedronToVTU(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, exportFolder);

        const vector<Eigen::MatrixXd> faceVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
        const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);
        const vector<Eigen::Matrix3d> faceRotationMatrices =
            geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices, faceNormals, faceTranslations);

        const vector<Eigen::MatrixXd> face2DVertices =
            geometryUtilities.PolyhedronFaceRotatedVertices(faceVertices, faceTranslations, faceRotationMatrices);
        const std::vector<std::vector<unsigned int>> facesUnalignedPoints =
            geometryUtilities.PolyhedronFacesUnalignedVertices(face2DVertices);

        const std::vector<std::vector<unsigned int>> polyhedronUnaligedFaces = {{1}, {2, 4}, {0}, {3}};

        const std::vector<unsigned int> unalignedVertices =
            geometryUtilities.UnalignedPolyhedronPoints(polyhedron.Vertices,
                                                        polyhedron.Faces,
                                                        faceTranslations,
                                                        faceRotationMatrices,
                                                        polyhedronUnaligedFaces,
                                                        facesUnalignedPoints);

        ASSERT_EQ(std::vector<unsigned int>({0, 1, 2, 4}), unalignedVertices);
    }
}

TEST(TestGeometryUtilities, TestPolyhedronInertia_ReferenceTetra)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-14;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check inertia of reference tetra 3D
        {
            const auto polyhedron = geometryUtilities.CreateTetrahedronWithVertices(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                    Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                    Eigen::Vector3d(0.0, 1.0, 0.0),
                                                                                    Eigen::Vector3d(0.0, 0.0, 1.0));

            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
            const vector<Eigen::MatrixXd> polyhedronFace3DVertices =
                geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);

            const vector<vector<unsigned int>> polyhedronFaceTriangulations =
                geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces, polyhedronFace3DVertices);

            const vector<unsigned int> polyhedronTetrahedronList =
                geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(polyhedron.Vertices,
                                                                             polyhedron.Faces,
                                                                             polyhedronFaceTriangulations,
                                                                             polyhedronBarycenter);
            const vector<Eigen::MatrixXd> polyhedronTetrahedronPoints =
                geometryUtilities.ExtractTetrahedronPoints(polyhedron.Vertices, polyhedronBarycenter, polyhedronTetrahedronList);

            const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
            const vector<bool> polyhedronFaceNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices, polyhedronBarycenter, polyhedronFaceNormals);
            const vector<Eigen::Vector3d> polyhedronFaceTranslations =
                geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
            const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

            const vector<Eigen::MatrixXd> polyhedronFace2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices, polyhedronFaceTranslations, polyhedronFaceRotationMatrices);

            const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints =
                geometryUtilities.PolyhedronFaceExtractTriangulationPoints(polyhedronFace2DVertices, polyhedronFaceTriangulations);

            const double polyhedronVolume =
                geometryUtilities.PolyhedronVolumeByBoundaryIntegral(polyhedronFace2DTriangulationPoints,
                                                                     polyhedronFaceNormals,
                                                                     polyhedronFaceNormalDirections,
                                                                     polyhedronFaceTranslations,
                                                                     polyhedronFaceRotationMatrices);

            const Eigen::Vector3d polyhedronCentroid = geometryUtilities.PolyhedronCentroid(polyhedronFace2DTriangulationPoints,
                                                                                            polyhedronFaceNormals,
                                                                                            polyhedronFaceNormalDirections,
                                                                                            polyhedronFaceTranslations,
                                                                                            polyhedronFaceRotationMatrices,
                                                                                            polyhedronVolume);

            const Eigen::Matrix3d polyhedronInertia =
                geometryUtilities.PolyhedronInertia(polyhedronCentroid, polyhedronTetrahedronPoints);

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronInertia(0, 0), +1.0 / 80.0, geometryUtilities.Tolerance1D()));
            ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronInertia(1, 1), +1.0 / 80.0, geometryUtilities.Tolerance1D()));
            ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronInertia(2, 2), +1.0 / 80.0, geometryUtilities.Tolerance1D()));
            ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronInertia(0, 1), +1.0 / 480.0, geometryUtilities.Tolerance1D()));
            ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronInertia(1, 0), +1.0 / 480.0, geometryUtilities.Tolerance1D()));
            ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronInertia(0, 2), +1.0 / 480.0, geometryUtilities.Tolerance1D()));
            ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronInertia(2, 0), +1.0 / 480.0, geometryUtilities.Tolerance1D()));
            ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronInertia(1, 2), +1.0 / 480.0, geometryUtilities.Tolerance1D()));
            ASSERT_TRUE(geometryUtilities.AreValuesEqual(polyhedronInertia(2, 1), +1.0 / 480.0, geometryUtilities.Tolerance1D()));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestPolyhedronInRadius)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    // check in radius of reference tetrahedron
    {
        const auto polyhedron = geometryUtilities.CreateTetrahedronWithVertices(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                Eigen::Vector3d(0.0, 1.0, 0.0),
                                                                                Eigen::Vector3d(0.0, 0.0, 1.0));
        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::MatrixXd> polyhedronFace3DVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);

        const auto polyhedronCentroidFacesDistance =
            geometryUtilities.PolyhedronCentroidFacesDistance(polyhedronBarycenter, polyhedronFaceNormals, polyhedronFace3DVertices);

        ASSERT_DOUBLE_EQ(sqrt(3.0) / 12.0, geometryUtilities.PolyhedronInRadius(polyhedronCentroidFacesDistance));
    }
}

TEST(TestGeometryUtilities, Test_Export_Polyhedron)
{
    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;

    const Gedim::PlatonicSolid platonicSolid = Gedim::PlatonicSolid(geometryUtilities, meshUtilities);

    const Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.dodecahedron();

    const Eigen::Vector3d polyhedron_centroid = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
    const vector<Eigen::MatrixXd> polyhedron_faces_vertices =
        geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
    const vector<Eigen::Vector3d> polyhedron_faces_centroid = geometryUtilities.PolyhedronFaceBarycenter(polyhedron_faces_vertices);
    const vector<vector<unsigned int>> polyhedron_faces_triangulation =
        geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces, polyhedron_faces_vertices);
    const vector<unsigned int> polyhedron_tetra =
        geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(polyhedron.Vertices,
                                                                     polyhedron.Faces,
                                                                     polyhedron_faces_triangulation,
                                                                     polyhedron_centroid);

    const vector<Eigen::Vector3d> polyhedron_faces_normal = geometryUtilities.PolyhedronFaceNormals(polyhedron_faces_vertices);
    const vector<Eigen::Vector3d> polyhedron_faces_translation =
        geometryUtilities.PolyhedronFaceTranslations(polyhedron_faces_vertices);
    const vector<Eigen::Matrix3d> polyhedron_faces_rotation =
        geometryUtilities.PolyhedronFaceRotationMatrices(polyhedron_faces_vertices, polyhedron_faces_normal, polyhedron_faces_translation);
    const vector<bool> polyhedron_faces_normal_direction =
        geometryUtilities.PolyhedronFaceNormalDirections(polyhedron_faces_vertices, polyhedron_centroid, polyhedron_faces_normal);

    const vector<Eigen::MatrixXd> polyhedron_faces_rotated_vertices =
        geometryUtilities.PolyhedronFaceRotatedVertices(polyhedron_faces_vertices, polyhedron_faces_translation, polyhedron_faces_rotation);
    const vector<Eigen::Vector3d> polyhedron_faces_rotated_centroid =
        geometryUtilities.PolyhedronFaceBarycenter(polyhedron_faces_rotated_vertices);
    const std::vector<std::vector<Eigen::Matrix3d>> polyhedron_faces_rotated_triangulation =
        geometryUtilities.PolyhedronFaceExtractTriangulationPoints(polyhedron_faces_rotated_vertices, polyhedron_faces_triangulation);

    const vector<Eigen::MatrixXd> tetrahedrons =
        geometryUtilities.ExtractTetrahedronPoints(polyhedron.Vertices, polyhedron_centroid, polyhedron_tetra);
    const double polyhedron_volume = geometryUtilities.PolyhedronVolumeByBoundaryIntegral(polyhedron_faces_rotated_triangulation,
                                                                                          polyhedron_faces_normal,
                                                                                          polyhedron_faces_normal_direction,
                                                                                          polyhedron_faces_translation,
                                                                                          polyhedron_faces_rotation);

    const std::vector<std::vector<Eigen::Matrix3d>> polyhedron_faces_triangulation_vertices =
        geometryUtilities.PolyhedronFaceExtractTriangulationPoints(polyhedron_faces_vertices, polyhedron_faces_triangulation);

    const auto tetra_reference_points =
        Gedim::Quadrature::Quadrature_Gauss3D_Tetrahedron_PositiveWeights::FillPointsAndWeights(2);
    Eigen::MatrixXd polyhedron_quadrature_points;

    const unsigned int polyhedron_num_tetra = tetrahedrons.size();

    const unsigned int tetra_num_reference_points = tetra_reference_points.Points.cols();
    const unsigned int polyhedron_num_quadrature_points = polyhedron_num_tetra * tetra_num_reference_points;

    Eigen::VectorXd quadrature_points_tetra_id(polyhedron_num_quadrature_points);
    polyhedron_quadrature_points.setZero(3, polyhedron_num_quadrature_points);

    Gedim::MapTetrahedron mapTetrahedron(geometryUtilities);

    for (unsigned int t = 0; t < polyhedron_num_tetra; t++)
    {
        const Eigen::MatrixXd &tetrahedronVertices = tetrahedrons[t];

        Gedim::MapTetrahedron::MapTetrahedronData mapTetrahedronData = mapTetrahedron.Compute(tetrahedronVertices);
        polyhedron_quadrature_points.block(0, tetra_num_reference_points * t, 3, tetra_num_reference_points) =
            mapTetrahedron.F(mapTetrahedronData, tetra_reference_points.Points);

        quadrature_points_tetra_id.segment(tetra_num_reference_points * t, tetra_num_reference_points).setConstant(t);
    }

    // Export to VTK
    std::string exportFolder = "./Export/TestGeometryUtilities/Test_Export_Polyhedron";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities vtpUtilities;

        //  original polyhedron
        vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

        vtpUtilities.Export(exportFolder + "/polyhedron.vtu", Gedim::VTKUtilities::Ascii);
    }

    {
        const unsigned int polyhedron_face_index = 0;

        const std::string export_polygon_folder = exportFolder + "/Polyhedron_Face_" + std::to_string(polyhedron_face_index);
        Gedim::Output::CreateFolder(export_polygon_folder);

        const auto polygon_edges_centroid =
            geometryUtilities.PolygonEdgesCentroid(polyhedron_faces_rotated_vertices.at(polyhedron_face_index));
        const auto polygon_edges_normal =
            geometryUtilities.PolygonEdgeNormals(polyhedron_faces_rotated_vertices.at(polyhedron_face_index));
        const auto polygon_edges_length =
            geometryUtilities.PolygonEdgeLengths(polyhedron_faces_rotated_vertices.at(polyhedron_face_index));

        geometryUtilities.ExportPolygonToVTU(polyhedron_face_index,
                                             polyhedron_faces_rotated_vertices.at(polyhedron_face_index),
                                             polyhedron_faces_rotated_triangulation.at(polyhedron_face_index),
                                             0.0,
                                             polyhedron_faces_rotated_centroid.at(polyhedron_face_index),
                                             polygon_edges_centroid,
                                             polygon_edges_normal,
                                             std::vector<bool>(polygon_edges_normal.cols(), true),
                                             export_polygon_folder);
    }

    {
        Gedim::VTKUtilities exporter;
        exporter.AddPoint(polyhedron_centroid);
        exporter.Export(exportFolder + "/polyhedron_centroid.vtu");
    }

    {
        Gedim::VTKUtilities vtkExperter;
        for (unsigned int t = 0; t < tetrahedrons.size(); t++)
        {
            Gedim::GeometryUtilities::Polyhedron subTetra =
                geometryUtilities.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
                                                                tetrahedrons[t].col(1),
                                                                tetrahedrons[t].col(2),
                                                                tetrahedrons[t].col(3));
            vector<double> id(1, t);

            vtkExperter.AddPolyhedron(subTetra.Vertices,
                                      subTetra.Edges,
                                      subTetra.Faces,
                                      {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()}});
        }

        vtkExperter.Export(exportFolder + "/polyhedron_tetra.vtu", Gedim::VTKUtilities::Ascii);
    }

    {
        Gedim::VTKUtilities exporter;
        exporter.AddPoints(polyhedron_quadrature_points,
                           {{"Id",
                             Gedim::VTPProperty::Formats::Cells,
                             static_cast<unsigned int>(quadrature_points_tetra_id.size()),
                             quadrature_points_tetra_id.data()}});
        exporter.Export(exportFolder + "/Polyhedron_quadrature.vtu");
    }

    const std::string export_polyhedron_folder = exportFolder + "/Polyhedron";
    Gedim::Output::CreateFolder(export_polyhedron_folder);
    geometryUtilities.ExportPolyhedronToVTU(0,
                                            polyhedron.Vertices,
                                            polyhedron.Edges,
                                            polyhedron.Faces,
                                            tetrahedrons,
                                            polyhedron_volume,
                                            polyhedron_centroid,
                                            polyhedron_faces_vertices,
                                            std::vector<double>(polyhedron.Faces.size(), 0.0),
                                            polyhedron_faces_rotated_centroid,
                                            polyhedron_faces_translation,
                                            polyhedron_faces_rotation,
                                            polyhedron_faces_triangulation_vertices,
                                            polyhedron_faces_centroid,
                                            polyhedron_faces_normal,
                                            polyhedron_faces_normal_direction,
                                            export_polyhedron_folder);
}

} // namespace GedimUnitTesting

#endif // __TEST_GEOMETRY_POLYHEDRON_H
