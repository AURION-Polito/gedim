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

#ifndef __TEST_MEDIT_Utilities_H
#define __TEST_MEDIT_Utilities_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "MeshMatrices_3D_329Cells_Mock.hpp"

#include "GeometryUtilities.hpp"
#include "IOUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "MEDIT_Utilities.hpp"

namespace GedimUnitTesting
{
// ***************************************************************************
TEST(Test_MEDIT_Utilities, MEDIT_Utilities_Test0Ds)
{
    std::string exportFolder = "./Export/Test_MEDIT_Utilities";
    Gedim::Output::CreateFolder(exportFolder);

    const unsigned int numGeometries = 4;

    Gedim::MEDIT_Utilities exporter;

    // Export to UCD
    Eigen::MatrixXd points(3, numGeometries);
    std::vector<unsigned int> references_id(numGeometries);

    for (unsigned int g = 0; g < numGeometries; g++)
    {
        points.col(g) << 1.0 + g, 0.0 + g, 0.0 + g;

        references_id[g] = g + 1;
    }

    exporter.ExportPoints(exportFolder + "/Geometry0Ds.mesh",
                          points,
                          references_id);
}
// ***************************************************************************
TEST(Test_MEDIT_Utilities, MEDIT_Utilities_Test1Ds)
{
    std::string exportFolder = "./Export/Test_MEDIT_Utilities";
    Gedim::Output::CreateFolder(exportFolder);

    const unsigned int numPoints = 4;
    const unsigned int numGeometries = 5;

    Gedim::MEDIT_Utilities exporter;

    // Export to UCD
    const Eigen::MatrixXd points =
        (Eigen::MatrixXd(3, numPoints) << 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0).finished();
    const Eigen::MatrixXi edges = (Eigen::MatrixXi(2, 5) << 0, 1, 2, 3, 0, 1, 2, 3, 0, 2).finished();

    std::vector<unsigned int> points_references_id(numPoints);
    std::vector<unsigned int> cells_references_id(numGeometries);

    for (unsigned int p = 0; p < numPoints; p++)
        points_references_id[p] = p + 1;

    for (unsigned int g = 0; g < numGeometries; g++)
        cells_references_id[g] = g + 1;

    exporter.ExportSegments(exportFolder + "/Geometry1Ds.inp",
                            points,
                            edges,
                            points_references_id,
                            cells_references_id);
}
// ***************************************************************************
TEST(Test_MEDIT_Utilities, MEDIT_Utilities_Test2Ds)
{
    std::string exportFolder = "./Export/Test_MEDIT_Utilities";
    Gedim::Output::CreateFolder(exportFolder);

    const unsigned int numPoints = 8;
    const unsigned int numGeometries = 3;

    Gedim::MEDIT_Utilities exporter;

    // Export to UCD
    const Eigen::MatrixXd points =
        (Eigen::MatrixXd(3, numPoints) << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 4.0, 4.0, 4.0, 4.0)
            .finished();
    const vector<vector<unsigned int>> polygons = {{0, 1, 2}, {0, 2, 3}, {4, 5, 6, 7}};

    std::vector<unsigned int> points_references_id(numPoints);
    std::vector<unsigned int> cells_references_id(numGeometries);

    for (unsigned int p = 0; p < numPoints; p++)
        points_references_id[p] = p + 1;

    for (unsigned int g = 0; g < numGeometries; g++)
        cells_references_id[g] = g + 1;

    exporter.ExportPolygons(exportFolder + "/Geometry2Ds.inp",
                            points,
                            polygons,
                            points_references_id,
                            cells_references_id);
}
// ***************************************************************************
TEST(Test_MEDIT_Utilities, MEDIT_Utilities_Test3D)
{
    std::string exportFolder = "./Export/Test_MEDIT_Utilities";
    Gedim::Output::CreateFolder(exportFolder);

    const unsigned int numPoints = 9;
    const unsigned int numFaces = 8;
    const unsigned int numGeometries = 2;

    Gedim::MEDIT_Utilities exporter;

    // Export to UCD
    const Eigen::MatrixXd points =
        (Eigen::MatrixXd(3, numPoints) << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0)
            .finished();
    const vector<vector<unsigned int>> faces = {{0, 2, 5}, {2, 1, 5}, { 1, 5, 0 }, { 0, 2, 1 }, { 5, 4, 8 }, { 4, 7, 8 }, { 7, 8, 5 }, { 5, 4, 7 }};
    const vector<vector<unsigned int>> polyhedra = {{0, 2, 1, 5}, {5, 4, 7, 8}};

    std::vector<unsigned int> points_references_id(numPoints);
    std::vector<unsigned int> faces_references_id(numGeometries);
    std::vector<unsigned int> cells_references_id(numGeometries);

    for (unsigned int p = 0; p < numPoints; p++)
        points_references_id[p] = p + 1;

    for (unsigned int f = 0; f < numFaces; f++)
        faces_references_id[f] = f + 1;

    for (unsigned int g = 0; g < numGeometries; g++)
        cells_references_id[g] = g + 1;

    exporter.ExportPolyhedra(exportFolder + "/Geometry3Ds.inp",
                             points,
                             faces,
                             polyhedra,
                             points_references_id,
                             faces_references_id,
                             cells_references_id);
}
// ***************************************************************************
TEST(Test_MEDIT_Utilities, MEDIT_Utilities_TestMesh3D)
{
    std::string exportFolder = "./Export/Test_MEDIT_Utilities/TestMesh3D";
    Gedim::Output::CreateFolder(exportFolder);

    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh(mockMesh.Mesh);

//    Gedim::MeshUtilities meshUtilities;
//    meshUtilities.ExportMeshToUCD(mesh, exportFolder, "Mesh3D", false);
}
// ***************************************************************************

} // namespace GedimUnitTesting

#endif // __TEST_MEDIT_Utilities_H
