#ifndef __TEST_GEOMETRY_SPLITPOLYHEDRON_H
#define __TEST_GEOMETRY_SPLITPOLYHEDRON_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "VTKUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

  TEST(TestGeometryUtilities, TestSplitPolyhedronWithPlane_TetraOne)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // split tetrahedron with plane
      {
        const Gedim::GeometryUtilities::Polyhedron polyhedron =  geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                               Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                                               Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                                                               Eigen::Vector3d(0.0, 1.0, 0.0));

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const Eigen::MatrixXd polyhedronEdgeTangents = geometryUtilities.PolyhedronEdgeTangents(polyhedron.Vertices,
                                                                                                polyhedron.Edges);
        const vector<Eigen::MatrixXd> polyhedronFaceVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                        polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFaceVertices);
        const vector<bool> polyhedronFaceNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFaceVertices,
                                                                                                             polyhedronBarycenter,
                                                                                                             polyhedronFaceNormals);
        const vector<vector<bool>> polyhedronFaceEdgeDirections = geometryUtilities.PolyhedronFaceEdgeDirections(polyhedron.Vertices,
                                                                                                                 polyhedron.Edges,
                                                                                                                 polyhedron.Faces);
        const vector<Eigen::MatrixXd> polyhedronFaceEdgesTangents = geometryUtilities.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                                                                                 polyhedron.Edges,
                                                                                                                 polyhedron.Faces,
                                                                                                                 polyhedronFaceEdgeDirections,
                                                                                                                 polyhedronEdgeTangents);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities.PolyhedronFaceTranslations(polyhedronFaceVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFaceVertices,
                                                                                                                        polyhedronFaceNormals,
                                                                                                                        polyhedronFaceTranslations);

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
          vtpUtilities.AddPolyhedron(polyhedron.Vertices,
                                     polyhedron.Edges,
                                     polyhedron.Faces);


          vtpUtilities.Export(exportFolder + "/Original.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result = geometryUtilities.SplitPolyhedronWithPlane(polyhedron.Vertices,
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

        const vector<Gedim::GeometryUtilities::Polyhedron> splitPolyhedra = geometryUtilities.SplitPolyhedronWithPlaneResultToPolyhedra(result);

        ASSERT_EQ(splitPolyhedra.size(), 2);

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[0].Vertices,
              splitPolyhedra[0].Edges,
              splitPolyhedra[0].Faces);


          vtpUtilities.Export(exportFolder + "/Positive.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[1].Vertices,
              splitPolyhedra[1].Edges,
              splitPolyhedra[1].Faces);


          vtpUtilities.Export(exportFolder + "/Negative.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split);
        ASSERT_EQ(result.Vertices.Vertices, (Eigen::MatrixXd(3, 7)<<
                                             0,   1,   0,   0, 0.5,   0,   0,
                                             0,   0,   1,   0,   0,   0, 0.5,
                                             0,   0,   0,   1, 0.5, 0.5, 0.5).finished());
        ASSERT_EQ(result.Vertices.NewVerticesOriginalEdge, std::vector<unsigned int>({ 4, 3, 5 }));
        ASSERT_EQ(result.Edges.Edges, (Eigen::MatrixXi(2, 12)<<
                                       4, 6, 5, 4, 3, 6, 1, 0, 0, 1, 5, 2,
                                       6, 5, 4, 3, 5, 3, 2, 2, 1, 4, 0, 6).finished());
        ASSERT_EQ(result.Edges.NewEdgesOriginalEdges, std::vector<int>({ -1, -1, -1, 4, 3, 5, 2, 1, 0, 4, 3, 5 }));
        ASSERT_EQ(result.Faces.Faces.size(), 8);

        ASSERT_EQ(result.Faces.Faces[0], (Eigen::MatrixXi(2, 3)<<
                                          4, 3, 5,
                                          3, 4, 2).finished());
        ASSERT_EQ(result.Faces.Faces[1], (Eigen::MatrixXi(2, 3)<<
                                          6, 3, 5,
                                          5, 4, 1).finished());
        ASSERT_EQ(result.Faces.Faces[2], (Eigen::MatrixXi(2, 3)<<
                                          6, 3, 4,
                                          5, 3, 0).finished());
        ASSERT_EQ(result.Faces.Faces[3], (Eigen::MatrixXi(2, 3)<<
                                          1, 2, 0,
                                          6, 7, 8).finished());
        ASSERT_EQ(result.Faces.Faces[4], (Eigen::MatrixXi(2, 4)<<
                                          1,  4,  5,  0,
                                          9,  2, 10,  8).finished());
        ASSERT_EQ(result.Faces.Faces[5], (Eigen::MatrixXi(2, 4)<<
                                          2,  6,  5,  0,
                                          11,  1, 10,  7).finished());
        ASSERT_EQ(result.Faces.Faces[6], (Eigen::MatrixXi(2, 4)<<
                                          2,  6,  4,  1,
                                          11,  0,  9,  6).finished());
        ASSERT_EQ(result.Faces.Faces[7], (Eigen::MatrixXi(2, 3)<<
                                          4, 6, 5,
                                          0, 1, 2).finished());
        ASSERT_EQ(result.Faces.NewFacesOriginalFaces, std::vector<int>({ 1,2,3,0,1,2,3,-1 }));

        ASSERT_EQ(result.PositivePolyhedron.Vertices, std::vector<unsigned int>({ 3,4,5,6 }));
        ASSERT_EQ(result.PositivePolyhedron.Edges, std::vector<unsigned int>({ 0,1,2,3,4,5 }));
        ASSERT_EQ(result.PositivePolyhedron.Faces, std::vector<unsigned int>({ 0,1,2,7 }));

        ASSERT_EQ(result.NegativePolyhedron.Vertices, std::vector<unsigned int>({ 0,1,2,4,5,6 }));
        ASSERT_EQ(result.NegativePolyhedron.Edges, std::vector<unsigned int>({ 0,1,2,6,7,8,9,10,11 }));
        ASSERT_EQ(result.NegativePolyhedron.Faces, std::vector<unsigned int>({ 3,4,5,6,7 }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
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
        const Gedim::GeometryUtilities::Polyhedron polyhedron =  geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                               Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                                               Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                                                               Eigen::Vector3d(0.0, 1.0, 0.0));

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const Eigen::MatrixXd polyhedronEdgeTangents = geometryUtilities.PolyhedronEdgeTangents(polyhedron.Vertices,
                                                                                                polyhedron.Edges);
        const vector<Eigen::MatrixXd> polyhedronFaceVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                        polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFaceVertices);
        const vector<bool> polyhedronFaceNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFaceVertices,
                                                                                                             polyhedronBarycenter,
                                                                                                             polyhedronFaceNormals);
        const vector<vector<bool>> polyhedronFaceEdgeDirections = geometryUtilities.PolyhedronFaceEdgeDirections(polyhedron.Vertices,
                                                                                                                 polyhedron.Edges,
                                                                                                                 polyhedron.Faces);
        const vector<Eigen::MatrixXd> polyhedronFaceEdgesTangents = geometryUtilities.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                                                                                 polyhedron.Edges,
                                                                                                                 polyhedron.Faces,
                                                                                                                 polyhedronFaceEdgeDirections,
                                                                                                                 polyhedronEdgeTangents);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities.PolyhedronFaceTranslations(polyhedronFaceVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFaceVertices,
                                                                                                                        polyhedronFaceNormals,
                                                                                                                        polyhedronFaceTranslations);

        const Eigen::Vector3d planeOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d planeNormal(-1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0));
        const Eigen::MatrixXd planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal);
        const Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

        const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result = geometryUtilities.SplitPolyhedronWithPlane(polyhedron.Vertices,
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

        const vector<Gedim::GeometryUtilities::Polyhedron> splitPolyhedra = geometryUtilities.SplitPolyhedronWithPlaneResultToPolyhedra(result);

        // Export to VTK
        std::string exportFolder = "./Export/TestSplitPolyhedronWithPlane/TetraTwo";
        Gedim::Output::CreateFolder(exportFolder);

        {
          Gedim::VTKUtilities vtpUtilities;

          //  original polyhedron
          vtpUtilities.AddPolyhedron(polyhedron.Vertices,
                                     polyhedron.Edges,
                                     polyhedron.Faces);


          vtpUtilities.Export(exportFolder + "/Original.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[0].Vertices,
              splitPolyhedra[0].Edges,
              splitPolyhedra[0].Faces);


          vtpUtilities.Export(exportFolder + "/Positive.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[1].Vertices,
              splitPolyhedra[1].Edges,
              splitPolyhedra[1].Faces);


          vtpUtilities.Export(exportFolder + "/Negative.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split);
        ASSERT_EQ(result.Vertices.Vertices, (Eigen::MatrixXd(3, 6)<<
                                             0,   1,   0,   0, 0.5,   0,
                                             0,   0,   1,   0,   0, 0.5,
                                             0,   0,   0,   1, 0.5, 0.5).finished());
        ASSERT_EQ(result.Vertices.NewVerticesOriginalEdge, std::vector<unsigned int>({ 4, 5 }));
        ASSERT_EQ(result.Edges.Edges, (Eigen::MatrixXi(2, 11)<<
                                       4, 5, 0, 4, 0, 5, 1, 0, 0, 1, 2,
                                       5, 0, 4, 3, 3, 3, 2, 2, 1, 4, 5).finished());

        ASSERT_EQ(result.Edges.NewEdgesOriginalEdges, std::vector<int>({ -1,-1,-1,4,3,5,2,1,0,4,5 }));
        ASSERT_EQ(result.Faces.Faces.size(), 8);
        ASSERT_EQ(result.Faces.Faces[0], (Eigen::MatrixXi(2, 3)<<
                                          4, 3, 0,
                                          3, 4, 2).finished());
        ASSERT_EQ(result.Faces.Faces[1], (Eigen::MatrixXi(2, 3)<<
                                          5, 3, 0,
                                          5, 4, 1).finished());
        ASSERT_EQ(result.Faces.Faces[2], (Eigen::MatrixXi(2, 3)<<
                                          5, 3, 4,
                                          5, 3, 0).finished());
        ASSERT_EQ(result.Faces.Faces[3], (Eigen::MatrixXi(2, 3)<<
                                          1, 2, 0,
                                          6, 7, 8).finished());
        ASSERT_EQ(result.Faces.Faces[4], (Eigen::MatrixXi(2, 3)<<
                                          1, 4, 0,
                                          9, 2, 8).finished());
        ASSERT_EQ(result.Faces.Faces[5], (Eigen::MatrixXi(2, 3)<<
                                          2, 5, 0,
                                          10, 1, 7).finished());
        ASSERT_EQ(result.Faces.Faces[6], (Eigen::MatrixXi(2, 4)<<
                                          2,  5,  4,  1,
                                          10,  0,  9,  6).finished());
        ASSERT_EQ(result.Faces.Faces[7], (Eigen::MatrixXi(2, 3)<<
                                          4, 5, 0,
                                          0, 1, 2).finished());
        ASSERT_EQ(result.Faces.NewFacesOriginalFaces, std::vector<int>({ 1,2,3,0,1,2,3,-1 }));

        ASSERT_EQ(result.PositivePolyhedron.Vertices, std::vector<unsigned int>({ 0,3,4,5 }));
        ASSERT_EQ(result.PositivePolyhedron.Edges, std::vector<unsigned int>({ 0,1,2,3,4,5 }));
        ASSERT_EQ(result.PositivePolyhedron.Faces, std::vector<unsigned int>({ 0,1,2,7 }));

        ASSERT_EQ(result.NegativePolyhedron.Vertices, std::vector<unsigned int>({ 0,1,2,4,5 }));
        ASSERT_EQ(result.NegativePolyhedron.Edges, std::vector<unsigned int>({ 0,1,2,6,7,8,9,10 }));
        ASSERT_EQ(result.NegativePolyhedron.Faces, std::vector<unsigned int>({ 3,4,5,6,7 }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
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
        const Gedim::GeometryUtilities::Polyhedron polyhedron =  geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                        1.0);

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const Eigen::MatrixXd polyhedronEdgeTangents = geometryUtilities.PolyhedronEdgeTangents(polyhedron.Vertices,
                                                                                                polyhedron.Edges);
        const vector<Eigen::MatrixXd> polyhedronFaceVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                        polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFaceVertices);
        const vector<bool> polyhedronFaceNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFaceVertices,
                                                                                                             polyhedronBarycenter,
                                                                                                             polyhedronFaceNormals);
        const vector<vector<bool>> polyhedronFaceEdgeDirections = geometryUtilities.PolyhedronFaceEdgeDirections(polyhedron.Vertices,
                                                                                                                 polyhedron.Edges,
                                                                                                                 polyhedron.Faces);
        const vector<Eigen::MatrixXd> polyhedronFaceEdgesTangents = geometryUtilities.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                                                                                 polyhedron.Edges,
                                                                                                                 polyhedron.Faces,
                                                                                                                 polyhedronFaceEdgeDirections,
                                                                                                                 polyhedronEdgeTangents);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities.PolyhedronFaceTranslations(polyhedronFaceVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFaceVertices,
                                                                                                                        polyhedronFaceNormals,
                                                                                                                        polyhedronFaceTranslations);

        const Eigen::Vector3d planeOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d planeNormal(-1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0));
        const Eigen::MatrixXd planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal);
        const Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

        const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result = geometryUtilities.SplitPolyhedronWithPlane(polyhedron.Vertices,
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

        const vector<Gedim::GeometryUtilities::Polyhedron> splitPolyhedra = geometryUtilities.SplitPolyhedronWithPlaneResultToPolyhedra(result);

        // Export to VTK
        std::string exportFolder = "./Export/TestSplitPolyhedronWithPlane/CubeOne";
        Gedim::Output::CreateFolder(exportFolder);

        {
          Gedim::VTKUtilities vtpUtilities;

          //  original polyhedron
          vtpUtilities.AddPolyhedron(polyhedron.Vertices,
                                     polyhedron.Edges,
                                     polyhedron.Faces);


          vtpUtilities.Export(exportFolder + "/Original.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[0].Vertices,
              splitPolyhedra[0].Edges,
              splitPolyhedra[0].Faces);


          vtpUtilities.Export(exportFolder + "/Positive.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[1].Vertices,
              splitPolyhedra[1].Edges,
              splitPolyhedra[1].Faces);


          vtpUtilities.Export(exportFolder + "/Negative.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split);
        ASSERT_EQ(result.Vertices.Vertices, (Eigen::MatrixXd(3, 8)<<
                                             0, 1, 1, 0, 0, 1, 1, 0,
                                             0, 0, 1, 1, 0, 0, 1, 1,
                                             0, 0, 0, 0, 1, 1, 1, 1).finished());
        ASSERT_EQ(result.Vertices.NewVerticesOriginalEdge, std::vector<unsigned int>({ }));
        ASSERT_EQ(result.Edges.Edges, (Eigen::MatrixXi(2, 15)<<
                                       5, 7, 0, 7, 4, 0, 1, 2, 3, 0, 6, 5, 3, 2, 1,
                                       7, 0, 5, 4, 5, 4, 2, 3, 0, 1, 7, 6, 7, 6, 5).finished());
        ASSERT_EQ(result.Edges.NewEdgesOriginalEdges, std::vector<int>({ -1, -1, -1, 7, 4, 8, 1, 2, 3, 0, 6, 5, 11, 10, 9 }));
        ASSERT_EQ(result.Faces.Faces.size(), 10);
        ASSERT_EQ(result.Faces.Faces[0], (Eigen::MatrixXi(2, 3)<<
                                          5, 7, 4,
                                          0, 3, 4).finished());
        ASSERT_EQ(result.Faces.Faces[1], (Eigen::MatrixXi(2, 3)<<
                                          7, 4, 0,
                                          3, 5, 1).finished());
        ASSERT_EQ(result.Faces.Faces[2], (Eigen::MatrixXi(2, 3)<<
                                          5, 4, 0,
                                          4, 5, 2).finished());
        ASSERT_EQ(result.Faces.Faces[3], (Eigen::MatrixXi(2, 4)<<
                                          1, 2, 3, 0,
                                          6, 7, 8, 9).finished());
        ASSERT_EQ(result.Faces.Faces[4], (Eigen::MatrixXi(2, 3)<<
                                          6 , 7,  5,
                                          10, 0, 11).finished());
        ASSERT_EQ(result.Faces.Faces[5], (Eigen::MatrixXi(2, 3)<<
                                          3,  7,  0,
                                          12,  1, 8).finished());
        ASSERT_EQ(result.Faces.Faces[6], (Eigen::MatrixXi(2, 4)<<
                                          2 , 6 , 5 , 1,
                                          13, 11,  14, 6).finished());
        ASSERT_EQ(result.Faces.Faces[7], (Eigen::MatrixXi(2, 3)<<
                                          1,  5,  0,
                                          14,  2,  9).finished());
        ASSERT_EQ(result.Faces.Faces[8], (Eigen::MatrixXi(2, 4)<<
                                          2,  6,  7,  3,
                                          13, 10, 12,  7).finished());
        ASSERT_EQ(result.Faces.Faces[9], (Eigen::MatrixXi(2, 3)<<
                                          5, 7, 0,
                                          0, 1, 2).finished());
        ASSERT_EQ(result.Faces.NewFacesOriginalFaces, std::vector<int>({ 1,2,4,0,1,2,3,4,5,-1 }));

        ASSERT_EQ(result.PositivePolyhedron.Vertices, std::vector<unsigned int>({ 0,4,5,7 }));
        ASSERT_EQ(result.PositivePolyhedron.Edges, std::vector<unsigned int>({ 0,1,2,3,4,5 }));
        ASSERT_EQ(result.PositivePolyhedron.Faces, std::vector<unsigned int>({ 0,1,2,9 }));

        ASSERT_EQ(result.NegativePolyhedron.Vertices, std::vector<unsigned int>({ 0,1,2,3,5,6,7 }));
        ASSERT_EQ(result.NegativePolyhedron.Edges, std::vector<unsigned int>({ 0,1,2,6,7,8,9,10,11,12,13,14 }));
        ASSERT_EQ(result.NegativePolyhedron.Faces, std::vector<unsigned int>({ 3,4,5,6,7,8,9 }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
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
        const Gedim::GeometryUtilities::Polyhedron polyhedron =  geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                        1.0);

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const Eigen::MatrixXd polyhedronEdgeTangents = geometryUtilities.PolyhedronEdgeTangents(polyhedron.Vertices,
                                                                                                polyhedron.Edges);
        const vector<Eigen::MatrixXd> polyhedronFaceVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                        polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFaceVertices);
        const vector<bool> polyhedronFaceNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFaceVertices,
                                                                                                             polyhedronBarycenter,
                                                                                                             polyhedronFaceNormals);
        const vector<vector<bool>> polyhedronFaceEdgeDirections = geometryUtilities.PolyhedronFaceEdgeDirections(polyhedron.Vertices,
                                                                                                                 polyhedron.Edges,
                                                                                                                 polyhedron.Faces);
        const vector<Eigen::MatrixXd> polyhedronFaceEdgesTangents = geometryUtilities.PolyhedronFaceEdgeTangents(polyhedron.Vertices,
                                                                                                                 polyhedron.Edges,
                                                                                                                 polyhedron.Faces,
                                                                                                                 polyhedronFaceEdgeDirections,
                                                                                                                 polyhedronEdgeTangents);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities.PolyhedronFaceTranslations(polyhedronFaceVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFaceVertices,
                                                                                                                        polyhedronFaceNormals,
                                                                                                                        polyhedronFaceTranslations);

        const Eigen::Vector3d planeOrigin(1.0, 0.0, 0.0);
        const Eigen::Vector3d planeNormal(1.0 / sqrt(2.0), 0.0, 1.0 / sqrt(2.0));
        const Eigen::MatrixXd planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal);
        const Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

        const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result = geometryUtilities.SplitPolyhedronWithPlane(polyhedron.Vertices,
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

        const vector<Gedim::GeometryUtilities::Polyhedron> splitPolyhedra = geometryUtilities.SplitPolyhedronWithPlaneResultToPolyhedra(result);

        // Export to VTK
        std::string exportFolder = "./Export/TestSplitPolyhedronWithPlane/CubeTwo";
        Gedim::Output::CreateFolder(exportFolder);

        {
          Gedim::VTKUtilities vtpUtilities;

          //  original polyhedron
          vtpUtilities.AddPolyhedron(polyhedron.Vertices,
                                     polyhedron.Edges,
                                     polyhedron.Faces);


          vtpUtilities.Export(exportFolder + "/Original.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[0].Vertices,
              splitPolyhedra[0].Edges,
              splitPolyhedra[0].Faces);


          vtpUtilities.Export(exportFolder + "/Positive.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        {
          Gedim::VTKUtilities vtpUtilities;

          //  positive polyhedron
          vtpUtilities.AddPolyhedron(splitPolyhedra[1].Vertices,
              splitPolyhedra[1].Edges,
              splitPolyhedra[1].Faces);


          vtpUtilities.Export(exportFolder + "/Negative.vtu",
                              Gedim::VTKUtilities::Ascii);
        }

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split);
        ASSERT_EQ(result.Vertices.Vertices, (Eigen::MatrixXd(3, 8)<<
                                             0, 1, 1, 0, 0, 1, 1, 0,
                                             0, 0, 1, 1, 0, 0, 1, 1,
                                             0, 0, 0, 0, 1, 1, 1, 1).finished());
        ASSERT_EQ(result.Vertices.NewVerticesOriginalEdge, std::vector<unsigned int>({ }));
        ASSERT_EQ(result.Edges.Edges, (Eigen::MatrixXi(2, 14)<<
                                       2, 7, 4, 1, 5, 6, 4, 2, 1, 2, 3, 0, 3, 0,
                                       7, 4, 1, 2, 6, 7, 5, 6, 5, 3, 0, 1, 7, 4).finished());
        ASSERT_EQ(result.Edges.NewEdgesOriginalEdges, std::vector<int>({ -1, 7, -1, 1, 5, 6, 4, 10, 9, 2, 3, 0, 11, 8 }));

        ASSERT_EQ(result.Faces.Faces.size(), 9);
        ASSERT_EQ(result.Faces.Faces[0], (Eigen::MatrixXi(2, 4)<<
                                          5, 6, 7, 4,
                                          4, 5, 1, 6).finished());
        ASSERT_EQ(result.Faces.Faces[1], (Eigen::MatrixXi(2, 4)<<
                                          2, 6, 5, 1,
                                          7, 4, 8, 3).finished());
        ASSERT_EQ(result.Faces.Faces[2], (Eigen::MatrixXi(2, 3)<<
                                          5, 4, 1,
                                          6, 2, 8).finished());
        ASSERT_EQ(result.Faces.Faces[3], (Eigen::MatrixXi(2, 3)<<
                                          6, 7, 2,
                                          5, 0, 7).finished());
        ASSERT_EQ(result.Faces.Faces[4], (Eigen::MatrixXi(2, 4)<<
                                          1,  2,  3,  0,
                                          3,  9, 10, 11).finished());
        ASSERT_EQ(result.Faces.Faces[5], (Eigen::MatrixXi(2, 4)<<
                                          3,  7,  4,  0,
                                          12,  1, 13, 10).finished());
        ASSERT_EQ(result.Faces.Faces[6], (Eigen::MatrixXi(2, 3)<<
                                          1,  4,  0,
                                          2, 13, 11).finished());
        ASSERT_EQ(result.Faces.Faces[7], (Eigen::MatrixXi(2, 3)<<
                                          2,  7,  3,
                                          0, 12,  9).finished());
        ASSERT_EQ(result.Faces.Faces[8], (Eigen::MatrixXi(2, 4)<<
                                          2, 7, 4, 1,
                                          0, 1, 2, 3).finished());
        ASSERT_EQ(result.Faces.NewFacesOriginalFaces, std::vector<int>({ 1,3,4,5,0,2,4,5,-1 }));

        ASSERT_EQ(result.PositivePolyhedron.Vertices, std::vector<unsigned int>({ 1,2,4,5,6,7 }));
        ASSERT_EQ(result.PositivePolyhedron.Edges, std::vector<unsigned int>({ 0,1,2,3,4,5,6,7,8 }));
        ASSERT_EQ(result.PositivePolyhedron.Faces, std::vector<unsigned int>({ 0,1,2,3,8 }));

        ASSERT_EQ(result.NegativePolyhedron.Vertices, std::vector<unsigned int>({ 0,1,2,3,4,7 }));
        ASSERT_EQ(result.NegativePolyhedron.Edges, std::vector<unsigned int>({ 0,1,2,3,9,10,11,12,13 }));
        ASSERT_EQ(result.NegativePolyhedron.Faces, std::vector<unsigned int>({ 4,5,6,7,8 }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_SPLIT_H
