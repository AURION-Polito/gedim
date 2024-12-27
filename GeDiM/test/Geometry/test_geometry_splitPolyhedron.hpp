#ifndef __TEST_GEOMETRY_SPLITPOLYHEDRON_H
#define __TEST_GEOMETRY_SPLITPOLYHEDRON_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"
#include "RefinementUtilities.hpp"
#include "VTKUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

TEST(TestGeometryUtilities, TestSplitPolyhedronWithPlane_TetraOne)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // split tetrahedron with plane
        {
            const Gedim::GeometryUtilities::Polyhedron polyhedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));

            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
            const Eigen::MatrixXd polyhedronEdgeTangents =
                geometryUtilities.PolyhedronEdgeTangents(polyhedron.Vertices, polyhedron.Edges);
            const vector<Eigen::MatrixXd> polyhedronFaceVertices =
                geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
            const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFaceVertices);
            const vector<bool> polyhedronFaceNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFaceVertices, polyhedronBarycenter, polyhedronFaceNormals);
            const vector<vector<bool>> polyhedronFaceEdgeDirections =
                geometryUtilities.PolyhedronFaceEdgeDirections(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
            const vector<Eigen::MatrixXd> polyhedronFaceEdgesTangents =
                geometryUtilities.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                             polyhedron.Edges,
                                                             polyhedron.Faces,
                                                             polyhedronFaceEdgeDirections,
                                                             polyhedronEdgeTangents);
            const vector<Eigen::Vector3d> polyhedronFaceTranslations =
                geometryUtilities.PolyhedronFaceTranslations(polyhedronFaceVertices);
            const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFaceVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

            const Eigen::Vector3d planeOrigin(0.0, 0.0, 0.5);
            const Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
            const Eigen::MatrixXd planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal);
            const Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

            // Export to VTK
            std::string exportFolder = "./Export/TestSplitPolyhedronWithPlane/TetraOne";
            Gedim::Output::CreateFolder(exportFolder);

            {
                Gedim::VTKUtilities vtpUtilities;

                //  original polyhedron
                vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

                vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
            }

            const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result =
                geometryUtilities.SplitPolyhedronWithPlane(polyhedron.Vertices,
                                                           polyhedron.Edges,
                                                           polyhedron.Faces,
                                                           polyhedronFaceVertices,
                                                           polyhedronFaceEdgesTangents,
                                                           polyhedronFaceTranslations,
                                                           polyhedronFaceRotationMatrices,
                                                           planeNormal,
                                                           planeOrigin,
                                                           planeRotationMatrix,
                                                           planeTranslation);

            const vector<Gedim::GeometryUtilities::Polyhedron> splitPolyhedra =
                geometryUtilities.SplitPolyhedronWithPlaneResultToPolyhedra(result);

            ASSERT_EQ(splitPolyhedra.size(), 2);

            {
                Gedim::VTKUtilities vtpUtilities;

                //  positive polyhedron
                vtpUtilities.AddPolyhedron(splitPolyhedra[0].Vertices, splitPolyhedra[0].Edges, splitPolyhedra[0].Faces);

                vtpUtilities.Export(exportFolder + "/Positive.vtu", Gedim::VTKUtilities::Ascii);
            }

            {
                Gedim::VTKUtilities vtpUtilities;

                //  positive polyhedron
                vtpUtilities.AddPolyhedron(splitPolyhedra[1].Vertices, splitPolyhedra[1].Edges, splitPolyhedra[1].Faces);

                vtpUtilities.Export(exportFolder + "/Negative.vtu", Gedim::VTKUtilities::Ascii);
            }

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split);
            ASSERT_EQ(result.Vertices.Vertices,
                      (Eigen::MatrixXd(3, 7) << 0, 1, 0, 0, 0.5, 0, 0, 0, 0, 1, 0, 0, 0, 0.5, 0, 0, 0, 1, 0.5, 0.5, 0.5).finished());
            ASSERT_EQ(result.Vertices.NewVerticesOriginalEdge, std::vector<unsigned int>({4, 3, 5}));
            ASSERT_EQ(result.Edges.Edges,
                      (Eigen::MatrixXi(2, 12) << 4, 6, 5, 4, 3, 6, 1, 0, 0, 1, 5, 2, 6, 5, 4, 3, 5, 3, 2, 2, 1, 4, 0, 6).finished());
            ASSERT_EQ(result.Edges.NewEdgesOriginalEdges, std::vector<int>({-1, -1, -1, 4, 3, 5, 2, 1, 0, 4, 3, 5}));
            ASSERT_EQ(result.Faces.Faces.size(), 8);

            ASSERT_EQ(result.Faces.Faces[0], (Eigen::MatrixXi(2, 3) << 4, 3, 5, 3, 4, 2).finished());
            ASSERT_EQ(result.Faces.Faces[1], (Eigen::MatrixXi(2, 3) << 6, 3, 5, 5, 4, 1).finished());
            ASSERT_EQ(result.Faces.Faces[2], (Eigen::MatrixXi(2, 3) << 6, 3, 4, 5, 3, 0).finished());
            ASSERT_EQ(result.Faces.Faces[3], (Eigen::MatrixXi(2, 3) << 1, 2, 0, 6, 7, 8).finished());
            ASSERT_EQ(result.Faces.Faces[4], (Eigen::MatrixXi(2, 4) << 1, 4, 5, 0, 9, 2, 10, 8).finished());
            ASSERT_EQ(result.Faces.Faces[5], (Eigen::MatrixXi(2, 4) << 2, 6, 5, 0, 11, 1, 10, 7).finished());
            ASSERT_EQ(result.Faces.Faces[6], (Eigen::MatrixXi(2, 4) << 2, 6, 4, 1, 11, 0, 9, 6).finished());
            ASSERT_EQ(result.Faces.Faces[7], (Eigen::MatrixXi(2, 3) << 4, 6, 5, 0, 1, 2).finished());
            ASSERT_EQ(result.Faces.NewFacesOriginalFaces, std::vector<int>({1, 2, 3, 0, 1, 2, 3, -1}));

            ASSERT_EQ(result.PositivePolyhedron.Vertices, std::vector<unsigned int>({3, 4, 5, 6}));
            ASSERT_EQ(result.PositivePolyhedron.Edges, std::vector<unsigned int>({0, 1, 2, 3, 4, 5}));
            ASSERT_EQ(result.PositivePolyhedron.Faces, std::vector<unsigned int>({0, 1, 2, 7}));

            ASSERT_EQ(result.NegativePolyhedron.Vertices, std::vector<unsigned int>({0, 1, 2, 4, 5, 6}));
            ASSERT_EQ(result.NegativePolyhedron.Edges, std::vector<unsigned int>({0, 1, 2, 6, 7, 8, 9, 10, 11}));
            ASSERT_EQ(result.NegativePolyhedron.Faces, std::vector<unsigned int>({3, 4, 5, 6, 7}));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestSplitPolyhedronWithPlane_TetraTwo)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // split tetrahedron with plane
        {
            const Gedim::GeometryUtilities::Polyhedron polyhedron =
                geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                              Eigen::Vector3d(1.0, 0.0, 0.0),
                                                              Eigen::Vector3d(0.0, 0.0, 1.0),
                                                              Eigen::Vector3d(0.0, 1.0, 0.0));

            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
            const Eigen::MatrixXd polyhedronEdgeTangents =
                geometryUtilities.PolyhedronEdgeTangents(polyhedron.Vertices, polyhedron.Edges);
            const vector<Eigen::MatrixXd> polyhedronFaceVertices =
                geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
            const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFaceVertices);
            const vector<bool> polyhedronFaceNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFaceVertices, polyhedronBarycenter, polyhedronFaceNormals);
            const vector<vector<bool>> polyhedronFaceEdgeDirections =
                geometryUtilities.PolyhedronFaceEdgeDirections(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
            const vector<Eigen::MatrixXd> polyhedronFaceEdgesTangents =
                geometryUtilities.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                             polyhedron.Edges,
                                                             polyhedron.Faces,
                                                             polyhedronFaceEdgeDirections,
                                                             polyhedronEdgeTangents);
            const vector<Eigen::Vector3d> polyhedronFaceTranslations =
                geometryUtilities.PolyhedronFaceTranslations(polyhedronFaceVertices);
            const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFaceVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

            const Eigen::Vector3d planeOrigin(0.0, 0.0, 0.0);
            const Eigen::Vector3d planeNormal(-1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0));
            const Eigen::MatrixXd planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal);
            const Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

            const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result =
                geometryUtilities.SplitPolyhedronWithPlane(polyhedron.Vertices,
                                                           polyhedron.Edges,
                                                           polyhedron.Faces,
                                                           polyhedronFaceVertices,
                                                           polyhedronFaceEdgesTangents,
                                                           polyhedronFaceTranslations,
                                                           polyhedronFaceRotationMatrices,
                                                           planeNormal,
                                                           planeOrigin,
                                                           planeRotationMatrix,
                                                           planeTranslation);

            const vector<Gedim::GeometryUtilities::Polyhedron> splitPolyhedra =
                geometryUtilities.SplitPolyhedronWithPlaneResultToPolyhedra(result);

            // Export to VTK
            std::string exportFolder = "./Export/TestSplitPolyhedronWithPlane/TetraTwo";
            Gedim::Output::CreateFolder(exportFolder);

            {
                Gedim::VTKUtilities vtpUtilities;

                //  original polyhedron
                vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

                vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
            }

            {
                Gedim::VTKUtilities vtpUtilities;

                //  positive polyhedron
                vtpUtilities.AddPolyhedron(splitPolyhedra[0].Vertices, splitPolyhedra[0].Edges, splitPolyhedra[0].Faces);

                vtpUtilities.Export(exportFolder + "/Positive.vtu", Gedim::VTKUtilities::Ascii);
            }

            {
                Gedim::VTKUtilities vtpUtilities;

                //  positive polyhedron
                vtpUtilities.AddPolyhedron(splitPolyhedra[1].Vertices, splitPolyhedra[1].Edges, splitPolyhedra[1].Faces);

                vtpUtilities.Export(exportFolder + "/Negative.vtu", Gedim::VTKUtilities::Ascii);
            }

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split);
            ASSERT_EQ(result.Vertices.Vertices,
                      (Eigen::MatrixXd(3, 6) << 0, 1, 0, 0, 0.5, 0, 0, 0, 1, 0, 0, 0.5, 0, 0, 0, 1, 0.5, 0.5).finished());
            ASSERT_EQ(result.Vertices.NewVerticesOriginalEdge, std::vector<unsigned int>({4, 5}));
            ASSERT_EQ(result.Edges.Edges,
                      (Eigen::MatrixXi(2, 11) << 4, 5, 0, 4, 0, 5, 1, 0, 0, 1, 2, 5, 0, 4, 3, 3, 3, 2, 2, 1, 4, 5).finished());

            ASSERT_EQ(result.Edges.NewEdgesOriginalEdges, std::vector<int>({-1, -1, -1, 4, 3, 5, 2, 1, 0, 4, 5}));
            ASSERT_EQ(result.Faces.Faces.size(), 8);
            ASSERT_EQ(result.Faces.Faces[0], (Eigen::MatrixXi(2, 3) << 4, 3, 0, 3, 4, 2).finished());
            ASSERT_EQ(result.Faces.Faces[1], (Eigen::MatrixXi(2, 3) << 5, 3, 0, 5, 4, 1).finished());
            ASSERT_EQ(result.Faces.Faces[2], (Eigen::MatrixXi(2, 3) << 5, 3, 4, 5, 3, 0).finished());
            ASSERT_EQ(result.Faces.Faces[3], (Eigen::MatrixXi(2, 3) << 1, 2, 0, 6, 7, 8).finished());
            ASSERT_EQ(result.Faces.Faces[4], (Eigen::MatrixXi(2, 3) << 1, 4, 0, 9, 2, 8).finished());
            ASSERT_EQ(result.Faces.Faces[5], (Eigen::MatrixXi(2, 3) << 2, 5, 0, 10, 1, 7).finished());
            ASSERT_EQ(result.Faces.Faces[6], (Eigen::MatrixXi(2, 4) << 2, 5, 4, 1, 10, 0, 9, 6).finished());
            ASSERT_EQ(result.Faces.Faces[7], (Eigen::MatrixXi(2, 3) << 4, 5, 0, 0, 1, 2).finished());
            ASSERT_EQ(result.Faces.NewFacesOriginalFaces, std::vector<int>({1, 2, 3, 0, 1, 2, 3, -1}));

            ASSERT_EQ(result.PositivePolyhedron.Vertices, std::vector<unsigned int>({0, 3, 4, 5}));
            ASSERT_EQ(result.PositivePolyhedron.Edges, std::vector<unsigned int>({0, 1, 2, 3, 4, 5}));
            ASSERT_EQ(result.PositivePolyhedron.Faces, std::vector<unsigned int>({0, 1, 2, 7}));

            ASSERT_EQ(result.NegativePolyhedron.Vertices, std::vector<unsigned int>({0, 1, 2, 4, 5}));
            ASSERT_EQ(result.NegativePolyhedron.Edges, std::vector<unsigned int>({0, 1, 2, 6, 7, 8, 9, 10}));
            ASSERT_EQ(result.NegativePolyhedron.Faces, std::vector<unsigned int>({3, 4, 5, 6, 7}));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestSplitPolyhedronWithPlane_CubeOne)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // split cube with plane
        {
            const Gedim::GeometryUtilities::Polyhedron polyhedron =
                geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);

            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
            const Eigen::MatrixXd polyhedronEdgeTangents =
                geometryUtilities.PolyhedronEdgeTangents(polyhedron.Vertices, polyhedron.Edges);
            const vector<Eigen::MatrixXd> polyhedronFaceVertices =
                geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
            const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFaceVertices);
            const vector<bool> polyhedronFaceNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFaceVertices, polyhedronBarycenter, polyhedronFaceNormals);
            const vector<vector<bool>> polyhedronFaceEdgeDirections =
                geometryUtilities.PolyhedronFaceEdgeDirections(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
            const vector<Eigen::MatrixXd> polyhedronFaceEdgesTangents =
                geometryUtilities.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                             polyhedron.Edges,
                                                             polyhedron.Faces,
                                                             polyhedronFaceEdgeDirections,
                                                             polyhedronEdgeTangents);
            const vector<Eigen::Vector3d> polyhedronFaceTranslations =
                geometryUtilities.PolyhedronFaceTranslations(polyhedronFaceVertices);
            const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFaceVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

            const Eigen::Vector3d planeOrigin(0.0, 0.0, 0.0);
            const Eigen::Vector3d planeNormal(-1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0));
            const Eigen::MatrixXd planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal);
            const Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

            const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result =
                geometryUtilities.SplitPolyhedronWithPlane(polyhedron.Vertices,
                                                           polyhedron.Edges,
                                                           polyhedron.Faces,
                                                           polyhedronFaceVertices,
                                                           polyhedronFaceEdgesTangents,
                                                           polyhedronFaceTranslations,
                                                           polyhedronFaceRotationMatrices,
                                                           planeNormal,
                                                           planeOrigin,
                                                           planeRotationMatrix,
                                                           planeTranslation);

            const vector<Gedim::GeometryUtilities::Polyhedron> splitPolyhedra =
                geometryUtilities.SplitPolyhedronWithPlaneResultToPolyhedra(result);

            // Export to VTK
            std::string exportFolder = "./Export/TestSplitPolyhedronWithPlane/CubeOne";
            Gedim::Output::CreateFolder(exportFolder);

            {
                Gedim::VTKUtilities vtpUtilities;

                //  original polyhedron
                vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

                vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
            }

            {
                Gedim::VTKUtilities vtpUtilities;

                //  positive polyhedron
                vtpUtilities.AddPolyhedron(splitPolyhedra[0].Vertices, splitPolyhedra[0].Edges, splitPolyhedra[0].Faces);

                vtpUtilities.Export(exportFolder + "/Positive.vtu", Gedim::VTKUtilities::Ascii);
            }

            {
                Gedim::VTKUtilities vtpUtilities;

                //  positive polyhedron
                vtpUtilities.AddPolyhedron(splitPolyhedra[1].Vertices, splitPolyhedra[1].Edges, splitPolyhedra[1].Faces);

                vtpUtilities.Export(exportFolder + "/Negative.vtu", Gedim::VTKUtilities::Ascii);
            }

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split);
            ASSERT_EQ(result.Vertices.Vertices,
                      (Eigen::MatrixXd(3, 8) << 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1).finished());
            ASSERT_EQ(result.Vertices.NewVerticesOriginalEdge, std::vector<unsigned int>({}));
            ASSERT_EQ(result.Edges.Edges,
                      (Eigen::MatrixXi(2, 15) << 5, 7, 0, 7, 4, 0, 1, 2, 3, 0, 6, 5, 3, 2, 1, 7, 0, 5, 4, 5, 4, 2, 3, 0, 1, 7, 6, 7, 6, 5)
                          .finished());
            ASSERT_EQ(result.Edges.NewEdgesOriginalEdges, std::vector<int>({-1, -1, -1, 7, 4, 8, 1, 2, 3, 0, 6, 5, 11, 10, 9}));
            ASSERT_EQ(result.Faces.Faces.size(), 10);
            ASSERT_EQ(result.Faces.Faces[0], (Eigen::MatrixXi(2, 3) << 5, 7, 4, 0, 3, 4).finished());
            ASSERT_EQ(result.Faces.Faces[1], (Eigen::MatrixXi(2, 3) << 7, 4, 0, 3, 5, 1).finished());
            ASSERT_EQ(result.Faces.Faces[2], (Eigen::MatrixXi(2, 3) << 5, 4, 0, 4, 5, 2).finished());
            ASSERT_EQ(result.Faces.Faces[3], (Eigen::MatrixXi(2, 4) << 1, 2, 3, 0, 6, 7, 8, 9).finished());
            ASSERT_EQ(result.Faces.Faces[4], (Eigen::MatrixXi(2, 3) << 6, 7, 5, 10, 0, 11).finished());
            ASSERT_EQ(result.Faces.Faces[5], (Eigen::MatrixXi(2, 3) << 3, 7, 0, 12, 1, 8).finished());
            ASSERT_EQ(result.Faces.Faces[6], (Eigen::MatrixXi(2, 4) << 2, 6, 5, 1, 13, 11, 14, 6).finished());
            ASSERT_EQ(result.Faces.Faces[7], (Eigen::MatrixXi(2, 3) << 1, 5, 0, 14, 2, 9).finished());
            ASSERT_EQ(result.Faces.Faces[8], (Eigen::MatrixXi(2, 4) << 2, 6, 7, 3, 13, 10, 12, 7).finished());
            ASSERT_EQ(result.Faces.Faces[9], (Eigen::MatrixXi(2, 3) << 5, 7, 0, 0, 1, 2).finished());
            ASSERT_EQ(result.Faces.NewFacesOriginalFaces, std::vector<int>({1, 2, 4, 0, 1, 2, 3, 4, 5, -1}));

            ASSERT_EQ(result.PositivePolyhedron.Vertices, std::vector<unsigned int>({0, 4, 5, 7}));
            ASSERT_EQ(result.PositivePolyhedron.Edges, std::vector<unsigned int>({0, 1, 2, 3, 4, 5}));
            ASSERT_EQ(result.PositivePolyhedron.Faces, std::vector<unsigned int>({0, 1, 2, 9}));

            ASSERT_EQ(result.NegativePolyhedron.Vertices, std::vector<unsigned int>({0, 1, 2, 3, 5, 6, 7}));
            ASSERT_EQ(result.NegativePolyhedron.Edges, std::vector<unsigned int>({0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14}));
            ASSERT_EQ(result.NegativePolyhedron.Faces, std::vector<unsigned int>({3, 4, 5, 6, 7, 8, 9}));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestSplitPolyhedronWithPlane_CubeTwo)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // split cube with plane
        {
            const Gedim::GeometryUtilities::Polyhedron polyhedron =
                geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);

            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
            const Eigen::MatrixXd polyhedronEdgeTangents =
                geometryUtilities.PolyhedronEdgeTangents(polyhedron.Vertices, polyhedron.Edges);
            const vector<Eigen::MatrixXd> polyhedronFaceVertices =
                geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
            const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFaceVertices);
            const vector<bool> polyhedronFaceNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFaceVertices, polyhedronBarycenter, polyhedronFaceNormals);
            const vector<vector<bool>> polyhedronFaceEdgeDirections =
                geometryUtilities.PolyhedronFaceEdgeDirections(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
            const vector<Eigen::MatrixXd> polyhedronFaceEdgesTangents =
                geometryUtilities.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                             polyhedron.Edges,
                                                             polyhedron.Faces,
                                                             polyhedronFaceEdgeDirections,
                                                             polyhedronEdgeTangents);
            const vector<Eigen::Vector3d> polyhedronFaceTranslations =
                geometryUtilities.PolyhedronFaceTranslations(polyhedronFaceVertices);
            const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFaceVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

            const Eigen::Vector3d planeOrigin(1.0, 0.0, 0.0);
            const Eigen::Vector3d planeNormal(1.0 / sqrt(2.0), 0.0, 1.0 / sqrt(2.0));
            const Eigen::MatrixXd planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal);
            const Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

            const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result =
                geometryUtilities.SplitPolyhedronWithPlane(polyhedron.Vertices,
                                                           polyhedron.Edges,
                                                           polyhedron.Faces,
                                                           polyhedronFaceVertices,
                                                           polyhedronFaceEdgesTangents,
                                                           polyhedronFaceTranslations,
                                                           polyhedronFaceRotationMatrices,
                                                           planeNormal,
                                                           planeOrigin,
                                                           planeRotationMatrix,
                                                           planeTranslation);

            const vector<Gedim::GeometryUtilities::Polyhedron> splitPolyhedra =
                geometryUtilities.SplitPolyhedronWithPlaneResultToPolyhedra(result);

            // Export to VTK
            std::string exportFolder = "./Export/TestSplitPolyhedronWithPlane/CubeTwo";
            Gedim::Output::CreateFolder(exportFolder);

            {
                Gedim::VTKUtilities vtpUtilities;

                //  original polyhedron
                vtpUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);

                vtpUtilities.Export(exportFolder + "/Original.vtu", Gedim::VTKUtilities::Ascii);
            }

            {
                Gedim::VTKUtilities vtpUtilities;

                //  positive polyhedron
                vtpUtilities.AddPolyhedron(splitPolyhedra[0].Vertices, splitPolyhedra[0].Edges, splitPolyhedra[0].Faces);

                vtpUtilities.Export(exportFolder + "/Positive.vtu", Gedim::VTKUtilities::Ascii);
            }

            {
                Gedim::VTKUtilities vtpUtilities;

                //  positive polyhedron
                vtpUtilities.AddPolyhedron(splitPolyhedra[1].Vertices, splitPolyhedra[1].Edges, splitPolyhedra[1].Faces);

                vtpUtilities.Export(exportFolder + "/Negative.vtu", Gedim::VTKUtilities::Ascii);
            }

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split);
            ASSERT_EQ(result.Vertices.Vertices,
                      (Eigen::MatrixXd(3, 8) << 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1).finished());
            ASSERT_EQ(result.Vertices.NewVerticesOriginalEdge, std::vector<unsigned int>({}));
            ASSERT_EQ(result.Edges.Edges,
                      (Eigen::MatrixXi(2, 14) << 2, 7, 4, 1, 5, 6, 4, 2, 1, 2, 3, 0, 3, 0, 7, 4, 1, 2, 6, 7, 5, 6, 5, 3, 0, 1, 7, 4)
                          .finished());
            ASSERT_EQ(result.Edges.NewEdgesOriginalEdges, std::vector<int>({-1, 7, -1, 1, 5, 6, 4, 10, 9, 2, 3, 0, 11, 8}));

            ASSERT_EQ(result.Faces.Faces.size(), 9);
            ASSERT_EQ(result.Faces.Faces[0], (Eigen::MatrixXi(2, 4) << 5, 6, 7, 4, 4, 5, 1, 6).finished());
            ASSERT_EQ(result.Faces.Faces[1], (Eigen::MatrixXi(2, 4) << 2, 6, 5, 1, 7, 4, 8, 3).finished());
            ASSERT_EQ(result.Faces.Faces[2], (Eigen::MatrixXi(2, 3) << 5, 4, 1, 6, 2, 8).finished());
            ASSERT_EQ(result.Faces.Faces[3], (Eigen::MatrixXi(2, 3) << 6, 7, 2, 5, 0, 7).finished());
            ASSERT_EQ(result.Faces.Faces[4], (Eigen::MatrixXi(2, 4) << 1, 2, 3, 0, 3, 9, 10, 11).finished());
            ASSERT_EQ(result.Faces.Faces[5], (Eigen::MatrixXi(2, 4) << 3, 7, 4, 0, 12, 1, 13, 10).finished());
            ASSERT_EQ(result.Faces.Faces[6], (Eigen::MatrixXi(2, 3) << 1, 4, 0, 2, 13, 11).finished());
            ASSERT_EQ(result.Faces.Faces[7], (Eigen::MatrixXi(2, 3) << 2, 7, 3, 0, 12, 9).finished());
            ASSERT_EQ(result.Faces.Faces[8], (Eigen::MatrixXi(2, 4) << 2, 7, 4, 1, 0, 1, 2, 3).finished());
            ASSERT_EQ(result.Faces.NewFacesOriginalFaces, std::vector<int>({1, 3, 4, 5, 0, 2, 4, 5, -1}));

            ASSERT_EQ(result.PositivePolyhedron.Vertices, std::vector<unsigned int>({1, 2, 4, 5, 6, 7}));
            ASSERT_EQ(result.PositivePolyhedron.Edges, std::vector<unsigned int>({0, 1, 2, 3, 4, 5, 6, 7, 8}));
            ASSERT_EQ(result.PositivePolyhedron.Faces, std::vector<unsigned int>({0, 1, 2, 3, 8}));

            ASSERT_EQ(result.NegativePolyhedron.Vertices, std::vector<unsigned int>({0, 1, 2, 3, 4, 7}));
            ASSERT_EQ(result.NegativePolyhedron.Edges, std::vector<unsigned int>({0, 1, 2, 3, 9, 10, 11, 12, 13}));
            ASSERT_EQ(result.NegativePolyhedron.Faces, std::vector<unsigned int>({4, 5, 6, 7, 8}));
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestSplitPolyhedronWithPlane_AlignedTetrahedron)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities, meshUtilities);

    std::string exportFolder = "./Export/TestSplitPolyhedronWithPlane_AlignedTetrahedron";
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

        Gedim::Output::CreateFolder(exportFolder + "/Original_Polyhedron");
        geometryUtilities.ExportPolyhedronToVTU(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, exportFolder + "/Original_Polyhedron");

        const vector<Eigen::MatrixXd> polyhedronFacesVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFacesNormal = geometryUtilities.PolyhedronFaceNormals(polyhedronFacesVertices);
        const vector<Eigen::Vector3d> polyhedronFacesTranslation =
            geometryUtilities.PolyhedronFaceTranslations(polyhedronFacesVertices);
        const vector<Eigen::Matrix3d> polyhedronFacesRotationMatrix =
            geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFacesVertices, polyhedronFacesNormal, polyhedronFacesTranslation);

        const vector<Eigen::MatrixXd> polyhedronFaces2DVertices =
            geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFacesVertices, polyhedronFacesTranslation, polyhedronFacesRotationMatrix);

        const Eigen::MatrixXd polyhedronEdgesTangent =
            geometryUtilities.PolyhedronEdgeTangents(polyhedron.Vertices, polyhedron.Edges);
        const vector<vector<bool>> polyhedronFacesEdgesDirection =
            geometryUtilities.PolyhedronFaceEdgeDirections(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
        const vector<Eigen::MatrixXd> polyhedronFacesEdgesTangent =
            geometryUtilities.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                         polyhedron.Edges,
                                                         polyhedron.Faces,
                                                         polyhedronFacesEdgesDirection,
                                                         polyhedronEdgesTangent);

        const std::vector<std::vector<unsigned int>> facesUnalignedPoints =
            geometryUtilities.PolyhedronFacesUnalignedVertices(polyhedronFaces2DVertices);

        const std::vector<std::vector<unsigned int>> polyhedronUnaligedFaces = {{0}, {1}, {2}, {3, 4}};

        const std::vector<unsigned int> unalignedVertices =
            geometryUtilities.UnalignedPolyhedronPoints(polyhedron.Vertices,
                                                        polyhedron.Faces,
                                                        polyhedronFacesTranslation,
                                                        polyhedronFacesRotationMatrix,
                                                        polyhedronUnaligedFaces,
                                                        facesUnalignedPoints);

        // create unaligned tetrahedron
        const Gedim::GeometryUtilities::Polyhedron unalignedTetrahedron =
            geometryUtilities.CreateTetrahedronWithVertices(polyhedron.Vertices.col(unalignedVertices.at(0)),
                                                            polyhedron.Vertices.col(unalignedVertices.at(1)),
                                                            polyhedron.Vertices.col(unalignedVertices.at(2)),
                                                            polyhedron.Vertices.col(unalignedVertices.at(3)));

        Gedim::Output::CreateFolder(exportFolder + "/Unaligned_Polyhedron");
        geometryUtilities.ExportPolyhedronToVTU(unalignedTetrahedron.Vertices,
                                                unalignedTetrahedron.Edges,
                                                unalignedTetrahedron.Faces,
                                                exportFolder + "/Unaligned_Polyhedron");

        const Eigen::VectorXd unalignedEdgesLength =
            geometryUtilities.PolyhedronEdgesLength(unalignedTetrahedron.Vertices, unalignedTetrahedron.Edges);

        const Gedim::RefinementUtilities::TetrahedronMaxEdgeDirection direction =
            refinementUtilities.ComputeTetrahedronMaxEdgeDirection(unalignedTetrahedron.Edges, unalignedEdgesLength);
        EXPECT_EQ(2, direction.MaxEdgeIndex);
        EXPECT_EQ(0, direction.OppositeVerticesIndex[0]);
        EXPECT_EQ(3, direction.OppositeVerticesIndex[1]);

        const Eigen::Vector3d planeOrigin =
            0.5 * (unalignedTetrahedron.Vertices.col(unalignedTetrahedron.Edges(0, direction.MaxEdgeIndex)) +
                   unalignedTetrahedron.Vertices.col(unalignedTetrahedron.Edges(1, direction.MaxEdgeIndex)));

        Eigen::Matrix3d planeTriangle;
        planeTriangle.col(0) << planeOrigin;
        planeTriangle.col(1) << unalignedTetrahedron.Vertices.col(direction.OppositeVerticesIndex[1]);
        planeTriangle.col(2) << unalignedTetrahedron.Vertices.col(direction.OppositeVerticesIndex[0]);

        const Eigen::Vector3d planeNormal = geometryUtilities.PolygonNormal(planeTriangle);
        const Eigen::Matrix3d planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal);
        const Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

        const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result =
            geometryUtilities.SplitPolyhedronWithPlane(polyhedron.Vertices,
                                                       polyhedron.Edges,
                                                       polyhedron.Faces,
                                                       polyhedronFacesVertices,
                                                       polyhedronFacesEdgesTangent,
                                                       polyhedronFacesTranslation,
                                                       polyhedronFacesRotationMatrix,
                                                       planeNormal,
                                                       planeOrigin,
                                                       planeRotationMatrix,
                                                       planeTranslation);
        {
            Gedim::VTKUtilities exporter;
            exporter.AddPoints(result.Vertices.Vertices);
            exporter.Export(exportFolder + "/Split_Vertices.vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddSegments(result.Vertices.Vertices, result.Edges.Edges);
            exporter.Export(exportFolder + "/Split_Edges.vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            for (unsigned int f = 0; f < result.Faces.Faces.size(); f++)
            {
                std::vector<unsigned int> faceVertices(result.Faces.Faces[f].cols());
                for (unsigned int v = 0; v < faceVertices.size(); v++)
                    faceVertices[v] = result.Faces.Faces[f](0, v);

                exporter.AddPolygon(result.Vertices.Vertices, faceVertices);
            }
            exporter.Export(exportFolder + "/Split_Faces.vtu");
        }

        const vector<Gedim::GeometryUtilities::Polyhedron> splitPolyhedra =
            geometryUtilities.SplitPolyhedronWithPlaneResultToPolyhedra(result);

        Gedim::Output::CreateFolder(exportFolder + "/Positive_Polyhedron");
        geometryUtilities.ExportPolyhedronToVTU(splitPolyhedra[0].Vertices,
                                                splitPolyhedra[0].Edges,
                                                splitPolyhedra[0].Faces,
                                                exportFolder + "/Positive_Polyhedron");
        Gedim::Output::CreateFolder(exportFolder + "/Negative_Polyhedron");
        geometryUtilities.ExportPolyhedronToVTU(splitPolyhedra[1].Vertices,
                                                splitPolyhedra[1].Edges,
                                                splitPolyhedra[1].Faces,
                                                exportFolder + "/Negative_Polyhedron");

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split);
        ASSERT_EQ(
            result.Vertices.Vertices,
            (Eigen::MatrixXd(3, 12) << 0, 0.1, 1, 0.75, 0.25, 0, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 1, 0.5, 0, 0, 0, 0, 0.25, 0.75, 1, 0.5, 0, 0, 0, 0.5, 0.5)
                .finished());
        ASSERT_EQ(result.Vertices.NewVerticesOriginalEdge, std::vector<unsigned int>({3}));
        ASSERT_EQ(
            result.Edges.Edges,
            (Eigen::MatrixXi(2, 17) << 0, 0, 11, 7, 11, 4, 5, 6, 5, 9, 7, 1, 2, 3, 0, 9, 2, 7, 11, 9, 9, 4, 5, 6, 0, 10, 10, 10, 2, 3, 11, 1, 8, 8)
                .finished());
        ASSERT_EQ(result.Edges.NewEdgesOriginalEdges, std::vector<int>({7, -1, -1, 8, 3, 4, 5, 6, 11, 12, 13, 1, 2, 3, 0, 9, 10}));

        ASSERT_EQ(result.Faces.Faces.size(), 8);
        ASSERT_EQ(result.Faces.Faces[0], (Eigen::MatrixXi(2, 5) << 11, 4, 5, 6, 0, 4, 5, 6, 7, 1).finished());
        ASSERT_EQ(result.Faces.Faces[1], (Eigen::MatrixXi(2, 5) << 10, 5, 4, 11, 9, 8, 5, 4, 2, 9).finished());
        ASSERT_EQ(result.Faces.Faces[2], (Eigen::MatrixXi(2, 5) << 7, 10, 5, 6, 0, 10, 8, 6, 7, 0).finished());
        ASSERT_EQ(result.Faces.Faces[3], (Eigen::MatrixXi(2, 3) << 9, 10, 7, 9, 10, 3).finished());
        ASSERT_EQ(result.Faces.Faces[4], (Eigen::MatrixXi(2, 5) << 1, 2, 3, 11, 0, 11, 12, 13, 1, 14).finished());
        ASSERT_EQ(result.Faces.Faces[5], (Eigen::MatrixXi(2, 5) << 8, 9, 11, 3, 2, 15, 2, 13, 12, 16).finished());
        ASSERT_EQ(result.Faces.Faces[6], (Eigen::MatrixXi(2, 6) << 1, 2, 8, 9, 7, 0, 11, 16, 15, 3, 0, 14).finished());
        ASSERT_EQ(result.Faces.Faces[7], (Eigen::MatrixXi(2, 4) << 7, 0, 11, 9, 0, 1, 2, 3).finished());
        ASSERT_EQ(result.Faces.NewFacesOriginalFaces, std::vector<int>({0, 1, 3, 4, 0, 1, 2, -1}));

        ASSERT_EQ(result.PositivePolyhedron.Vertices, std::vector<unsigned int>({0, 4, 5, 6, 7, 9, 10, 11}));
        ASSERT_EQ(result.PositivePolyhedron.Edges, std::vector<unsigned int>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
        ASSERT_EQ(result.PositivePolyhedron.Faces, std::vector<unsigned int>({0, 1, 2, 3, 7}));

        ASSERT_EQ(result.NegativePolyhedron.Vertices, std::vector<unsigned int>({0, 1, 2, 3, 7, 8, 9, 11}));
        ASSERT_EQ(result.NegativePolyhedron.Edges, std::vector<unsigned int>({0, 1, 2, 3, 11, 12, 13, 14, 15, 16}));
        ASSERT_EQ(result.NegativePolyhedron.Faces, std::vector<unsigned int>({4, 5, 6, 7}));
    }
}
} // namespace GedimUnitTesting

#endif // __TEST_GEOMETRY_SPLIT_H
