#ifndef __TEST_MESH_H
#define __TEST_MESH_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MeshMatrices.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "MeshMatrices_2D_4Cells_Mock.hpp"
#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "MeshMatricesDAO.hpp"

#include "MeshDAOExporterToCsv.hpp"
#include "MeshFromCsvUtilities.hpp"
#include "MeshDAOImporterFromCsv.hpp"
#include "FileTextReader.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

  TEST(TestMesh, TestMeshMatricesDAO)
  {
    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    EXPECT_NO_THROW(meshDao.InitializeDimension(3));
    EXPECT_EQ(meshDao.Dimension(), 3);

    EXPECT_NO_THROW(meshDao.Cell0DsInitialize(2));
    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 2);
    EXPECT_NO_THROW(meshDao.Cell0DSetState(1, true));
    EXPECT_EQ(meshDao.Cell0DIsActive(0), false);
    EXPECT_EQ(meshDao.Cell0DIsActive(1), true);
    EXPECT_NO_THROW(meshDao.Cell0DSetMarker(1, 1));
    EXPECT_EQ(meshDao.Cell0DMarker(0), 0);
    EXPECT_EQ(meshDao.Cell0DMarker(1), 1);
    EXPECT_NO_THROW(meshDao.Cell0DInitializeDoubleProperties(1));
    EXPECT_EQ(meshDao.Cell0DAddDoubleProperty("Test"), 0);
    EXPECT_NO_THROW(meshDao.Cell0DInitializeDoublePropertyValues(0, 0, 3));
    EXPECT_NO_THROW(meshDao.Cell0DInsertDoublePropertyValue(0, 0, 0, 5.2));
    EXPECT_NO_THROW(meshDao.Cell0DInsertDoublePropertyValue(0, 0, 1, 7.4));
    EXPECT_NO_THROW(meshDao.Cell0DInsertDoublePropertyValue(0, 0, 2, 3));
    EXPECT_NO_THROW(meshDao.Cell0DInitializeDoublePropertyValues(1, 0, 2));
    EXPECT_NO_THROW(meshDao.Cell0DInsertDoublePropertyValue(1, 0, 0, 15.2));
    EXPECT_NO_THROW(meshDao.Cell0DInsertDoublePropertyValue(1, 0, 1, 17.4));
    EXPECT_EQ(meshDao.Cell0DNumberDoubleProperties(), 1);
    EXPECT_TRUE(meshDao.Cell0DDoublePropertyExists("Test"));
    EXPECT_EQ(meshDao.Cell0DDoublePropertyIndex("Test"), 0);
    EXPECT_FALSE(meshDao.Cell0DDoublePropertyExists("NoProperty"));
    EXPECT_EQ(meshDao.Cell0DDoublePropertyId(0), "Test");
    EXPECT_EQ(meshDao.Cell0DDoublePropertySize(0, 0), 3);
    EXPECT_EQ(meshDao.Cell0DDoublePropertyValue(0, 0, 0), 5.2);
    EXPECT_EQ(meshDao.Cell0DDoublePropertyValue(0, 0, 1), 7.4);
    EXPECT_EQ(meshDao.Cell0DDoublePropertyValue(0, 0, 2), 3);
    EXPECT_EQ(meshDao.Cell0DDoublePropertySize(1, 0), 2);
    EXPECT_EQ(meshDao.Cell0DDoublePropertyValue(1, 0, 0), 15.2);
    EXPECT_EQ(meshDao.Cell0DDoublePropertyValue(1, 0, 1), 17.4);
    EXPECT_NO_THROW(meshDao.Cell0DAppend(2));
    EXPECT_NO_THROW(meshDao.Cell0DInitializeDoublePropertyValues(2, 0, 1));
    EXPECT_NO_THROW(meshDao.Cell0DInsertDoublePropertyValue(2, 0, 0, 75.2));
    EXPECT_NO_THROW(meshDao.Cell0DRemove(1));

    EXPECT_NO_THROW(meshDao.Cell1DsInitialize(2));
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 2);
    EXPECT_NO_THROW(meshDao.Cell1DInsertExtremes(1, 1, 0));
    EXPECT_EQ(meshDao.Cell1DOrigin(0), 0);
    EXPECT_EQ(meshDao.Cell1DEnd(0), 0);
    EXPECT_EQ(meshDao.Cell1DOrigin(1), 1);
    EXPECT_EQ(meshDao.Cell1DEnd(1), 0);
    EXPECT_EQ(meshDao.Cell1DFindExtreme(1, 0), 1);
    EXPECT_EQ(meshDao.Cell1DFindExtreme(1, 1), 0);
    EXPECT_THROW(meshDao.Cell1DFindExtreme(1, 2),
                 std::runtime_error);
    EXPECT_NO_THROW(meshDao.Cell1DSetState(1, true));
    EXPECT_EQ(meshDao.Cell1DIsActive(0), false);
    EXPECT_EQ(meshDao.Cell1DIsActive(1), true);
    EXPECT_NO_THROW(meshDao.Cell1DSetMarker(1, 1));
    EXPECT_EQ(meshDao.Cell1DMarker(0), 0);
    EXPECT_EQ(meshDao.Cell1DMarker(1), 1);
    EXPECT_NO_THROW(meshDao.Cell1DInitializeDoubleProperties(1));
    EXPECT_EQ(meshDao.Cell1DAddDoubleProperty("Test"), 0);
    EXPECT_NO_THROW(meshDao.Cell1DInitializeDoublePropertyValues(0, 0, 3));
    EXPECT_NO_THROW(meshDao.Cell1DInsertDoublePropertyValue(0, 0, 0, 5.2));
    EXPECT_NO_THROW(meshDao.Cell1DInsertDoublePropertyValue(0, 0, 1, 7.4));
    EXPECT_NO_THROW(meshDao.Cell1DInsertDoublePropertyValue(0, 0, 2, 3));
    EXPECT_NO_THROW(meshDao.Cell1DInitializeDoublePropertyValues(1, 0, 2));
    EXPECT_NO_THROW(meshDao.Cell1DInsertDoublePropertyValue(1, 0, 0, 15.2));
    EXPECT_NO_THROW(meshDao.Cell1DInsertDoublePropertyValue(1, 0, 1, 17.4));
    EXPECT_EQ(meshDao.Cell1DNumberDoubleProperties(), 1);
    EXPECT_TRUE(meshDao.Cell1DDoublePropertyExists("Test"));
    EXPECT_EQ(meshDao.Cell1DDoublePropertyIndex("Test"), 0);
    EXPECT_FALSE(meshDao.Cell1DDoublePropertyExists("NoProperty"));
    EXPECT_EQ(meshDao.Cell1DDoublePropertyId(0), "Test");
    EXPECT_EQ(meshDao.Cell1DDoublePropertySize(0, 0), 3);
    EXPECT_EQ(meshDao.Cell1DDoublePropertyValue(0, 0, 0), 5.2);
    EXPECT_EQ(meshDao.Cell1DDoublePropertyValue(0, 0, 1), 7.4);
    EXPECT_EQ(meshDao.Cell1DDoublePropertyValue(0, 0, 2), 3);
    EXPECT_EQ(meshDao.Cell1DDoublePropertySize(1, 0), 2);
    EXPECT_EQ(meshDao.Cell1DDoublePropertyValue(1, 0, 0), 15.2);
    EXPECT_EQ(meshDao.Cell1DDoublePropertyValue(1, 0, 1), 17.4);
    EXPECT_NO_THROW(meshDao.Cell1DAppend(2));
    EXPECT_NO_THROW(meshDao.Cell1DInitializeDoublePropertyValues(2, 0, 1));
    EXPECT_NO_THROW(meshDao.Cell1DInsertDoublePropertyValue(2, 0, 0, 75.2));
    EXPECT_NO_THROW(meshDao.Cell1DRemove(1));

    EXPECT_NO_THROW(meshDao.Cell2DsInitialize(2));
    EXPECT_EQ(meshDao.Cell2DTotalNumber(), 2);
    EXPECT_NO_THROW(meshDao.Cell2DAddVertices(1, {1, 0}));
    EXPECT_EQ(meshDao.Cell2DNumberVertices(0), 0);
    EXPECT_EQ(meshDao.Cell2DNumberVertices(1), 2);
    EXPECT_EQ(meshDao.Cell2DFindVertex(1, 0), 1);
    EXPECT_EQ(meshDao.Cell2DFindVertex(1, 1), 0);
    EXPECT_THROW(meshDao.Cell2DFindVertex(1, 3),
                 std::runtime_error);
    EXPECT_EQ(meshDao.Cell2DVertex(1, 0), 1);
    EXPECT_EQ(meshDao.Cell2DVertex(1, 1), 0);
    EXPECT_NO_THROW(meshDao.Cell2DAddEdges(1, {1, 0}));
    EXPECT_EQ(meshDao.Cell2DFindEdge(1, 0), 1);
    EXPECT_EQ(meshDao.Cell2DFindEdge(1, 1), 0);
    EXPECT_THROW(meshDao.Cell2DFindEdge(1, 3),
                 std::runtime_error);
    EXPECT_EQ(meshDao.Cell2DNumberEdges(0), 0);
    EXPECT_EQ(meshDao.Cell2DNumberEdges(1), 2);
    EXPECT_EQ(meshDao.Cell2DEdge(1, 0), 1);
    EXPECT_EQ(meshDao.Cell2DEdge(1, 1), 0);
    EXPECT_NO_THROW(meshDao.Cell2DSetState(1, true));
    EXPECT_EQ(meshDao.Cell2DIsActive(0), false);
    EXPECT_EQ(meshDao.Cell2DIsActive(1), true);
    EXPECT_NO_THROW(meshDao.Cell2DSetMarker(1, 1));
    EXPECT_EQ(meshDao.Cell2DMarker(0), 0);
    EXPECT_EQ(meshDao.Cell2DMarker(1), 1);
    EXPECT_NO_THROW(meshDao.Cell2DInitializeDoubleProperties(1));
    EXPECT_EQ(meshDao.Cell2DAddDoubleProperty("Test"), 0);
    EXPECT_NO_THROW(meshDao.Cell2DInitializeDoublePropertyValues(0, 0, 3));
    EXPECT_NO_THROW(meshDao.Cell2DInsertDoublePropertyValue(0, 0, 0, 5.2));
    EXPECT_NO_THROW(meshDao.Cell2DInsertDoublePropertyValue(0, 0, 1, 7.4));
    EXPECT_NO_THROW(meshDao.Cell2DInsertDoublePropertyValue(0, 0, 2, 3));
    EXPECT_NO_THROW(meshDao.Cell2DInitializeDoublePropertyValues(1, 0, 2));
    EXPECT_NO_THROW(meshDao.Cell2DInsertDoublePropertyValue(1, 0, 0, 15.2));
    EXPECT_NO_THROW(meshDao.Cell2DInsertDoublePropertyValue(1, 0, 1, 17.4));
    EXPECT_EQ(meshDao.Cell2DNumberDoubleProperties(), 1);
    EXPECT_TRUE(meshDao.Cell2DDoublePropertyExists("Test"));
    EXPECT_EQ(meshDao.Cell2DDoublePropertyIndex("Test"), 0);
    EXPECT_FALSE(meshDao.Cell2DDoublePropertyExists("NoProperty"));
    EXPECT_EQ(meshDao.Cell2DDoublePropertyId(0), "Test");
    EXPECT_EQ(meshDao.Cell2DDoublePropertySize(0, 0), 3);
    EXPECT_EQ(meshDao.Cell2DDoublePropertyValue(0, 0, 0), 5.2);
    EXPECT_EQ(meshDao.Cell2DDoublePropertyValue(0, 0, 1), 7.4);
    EXPECT_EQ(meshDao.Cell2DDoublePropertyValue(0, 0, 2), 3);
    EXPECT_EQ(meshDao.Cell2DDoublePropertySize(1, 0), 2);
    EXPECT_EQ(meshDao.Cell2DDoublePropertyValue(1, 0, 0), 15.2);
    EXPECT_EQ(meshDao.Cell2DDoublePropertyValue(1, 0, 1), 17.4);
    EXPECT_FALSE(meshDao.Cell2DHasOriginalCell2D(1));
    EXPECT_NO_THROW(meshDao.Cell2DAppend(2));
    EXPECT_EQ(meshDao.Cell2DTotalNumber(), 4);
    EXPECT_FALSE(meshDao.Cell2DHasOriginalCell2D(1));
    EXPECT_NO_THROW(meshDao.Cell2DInsertUpdatedCell2D(1, 2));
    EXPECT_NO_THROW(meshDao.Cell2DInsertUpdatedCell2D(1, 3));
    EXPECT_TRUE(meshDao.Cell2DHasOriginalCell2D(3));
    EXPECT_EQ(meshDao.Cell2DOriginalCell2D(3), 1);
    EXPECT_NO_THROW(meshDao.Cell2DInitializeDoublePropertyValues(2, 0, 1));
    EXPECT_NO_THROW(meshDao.Cell2DInsertDoublePropertyValue(2, 0, 0, 75.2));
    EXPECT_NO_THROW(meshDao.Cell2DRemove(1));
    EXPECT_FALSE(meshDao.Cell2DHasOriginalCell2D(2));

    EXPECT_NO_THROW(meshDao.Cell3DsInitialize(2));
    EXPECT_EQ(meshDao.Cell3DTotalNumber(), 2);
    EXPECT_NO_THROW(meshDao.Cell3DAddVertices(1, {1, 0}));
    EXPECT_EQ(meshDao.Cell3DNumberVertices(0), 0);
    EXPECT_EQ(meshDao.Cell3DNumberVertices(1), 2);
    EXPECT_EQ(meshDao.Cell3DVertex(1, 0), 1);
    EXPECT_EQ(meshDao.Cell3DVertex(1, 1), 0);
    EXPECT_NO_THROW(meshDao.Cell3DAddEdges(1, {1, 0}));
    EXPECT_EQ(meshDao.Cell3DNumberEdges(0), 0);
    EXPECT_EQ(meshDao.Cell3DNumberEdges(1), 2);
    EXPECT_EQ(meshDao.Cell3DEdge(1, 0), 1);
    EXPECT_EQ(meshDao.Cell3DEdge(1, 1), 0);
    EXPECT_NO_THROW(meshDao.Cell3DAddFaces(1, {1, 0}));
    EXPECT_EQ(meshDao.Cell3DNumberFaces(0), 0);
    EXPECT_EQ(meshDao.Cell3DNumberFaces(1), 2);
    EXPECT_EQ(meshDao.Cell3DFace(1, 0), 1);
    EXPECT_EQ(meshDao.Cell3DFace(1, 1), 0);
    EXPECT_NO_THROW(meshDao.Cell3DSetState(1, true));
    EXPECT_EQ(meshDao.Cell3DIsActive(0), false);
    EXPECT_EQ(meshDao.Cell3DIsActive(1), true);
    EXPECT_NO_THROW(meshDao.Cell3DSetMarker(1, 1));
    EXPECT_EQ(meshDao.Cell3DMarker(0), 0);
    EXPECT_EQ(meshDao.Cell3DMarker(1), 1);
    EXPECT_NO_THROW(meshDao.Cell3DInitializeDoubleProperties(1));
    EXPECT_EQ(meshDao.Cell3DAddDoubleProperty("Test"), 0);
    EXPECT_NO_THROW(meshDao.Cell3DInitializeDoublePropertyValues(0, 0, 3));
    EXPECT_NO_THROW(meshDao.Cell3DInsertDoublePropertyValue(0, 0, 0, 5.2));
    EXPECT_NO_THROW(meshDao.Cell3DInsertDoublePropertyValue(0, 0, 1, 7.4));
    EXPECT_NO_THROW(meshDao.Cell3DInsertDoublePropertyValue(0, 0, 2, 3));
    EXPECT_NO_THROW(meshDao.Cell3DInitializeDoublePropertyValues(1, 0, 2));
    EXPECT_NO_THROW(meshDao.Cell3DInsertDoublePropertyValue(1, 0, 0, 15.2));
    EXPECT_NO_THROW(meshDao.Cell3DInsertDoublePropertyValue(1, 0, 1, 17.4));
    EXPECT_EQ(meshDao.Cell3DNumberDoubleProperties(), 1);
    EXPECT_TRUE(meshDao.Cell3DDoublePropertyExists("Test"));
    EXPECT_EQ(meshDao.Cell3DDoublePropertyIndex("Test"), 0);
    EXPECT_FALSE(meshDao.Cell3DDoublePropertyExists("NoProperty"));
    EXPECT_EQ(meshDao.Cell3DDoublePropertyId(0), "Test");
    EXPECT_EQ(meshDao.Cell3DDoublePropertySize(0, 0), 3);
    EXPECT_EQ(meshDao.Cell3DDoublePropertyValue(0, 0, 0), 5.2);
    EXPECT_EQ(meshDao.Cell3DDoublePropertyValue(0, 0, 1), 7.4);
    EXPECT_EQ(meshDao.Cell3DDoublePropertyValue(0, 0, 2), 3);
    EXPECT_EQ(meshDao.Cell3DDoublePropertySize(1, 0), 2);
    EXPECT_EQ(meshDao.Cell3DDoublePropertyValue(1, 0, 0), 15.2);
    EXPECT_EQ(meshDao.Cell3DDoublePropertyValue(1, 0, 1), 17.4);
    EXPECT_NO_THROW(meshDao.Cell3DAppend(2));
    EXPECT_NO_THROW(meshDao.Cell3DInitializeDoublePropertyValues(2, 0, 1));
    EXPECT_NO_THROW(meshDao.Cell3DInsertDoublePropertyValue(2, 0, 0, 75.2));
    EXPECT_NO_THROW(meshDao.Cell3DRemove(1));

    EXPECT_NO_THROW(meshDao.Cell1DInitializeNeighbourCell2Ds(1, 2));
    EXPECT_NO_THROW(meshDao.Cell1DInsertNeighbourCell2D(1, 0, 1));
    EXPECT_NO_THROW(meshDao.Cell1DInsertNeighbourCell2D(1, 1, 0));
    EXPECT_EQ(meshDao.Cell1DNumberNeighbourCell2D(0), 0);
    EXPECT_EQ(meshDao.Cell1DNumberNeighbourCell2D(1), 2);
    EXPECT_EQ(meshDao.Cell1DNeighbourCell2D(1, 0), 1);
    EXPECT_EQ(meshDao.Cell1DNeighbourCell2D(1, 1), 0);

    EXPECT_NO_THROW(meshDao.Cell0DRemove(0));
    EXPECT_NO_THROW(meshDao.Cell0DRemove(0));
    EXPECT_NO_THROW(meshDao.Cell0DRemove(0));
    EXPECT_NO_THROW(meshDao.Cell1DRemove(0));
    EXPECT_NO_THROW(meshDao.Cell1DRemove(0));
    EXPECT_NO_THROW(meshDao.Cell1DRemove(0));
    EXPECT_NO_THROW(meshDao.Cell2DRemove(0));
    EXPECT_NO_THROW(meshDao.Cell2DRemove(0));
    EXPECT_NO_THROW(meshDao.Cell2DRemove(0));
    EXPECT_NO_THROW(meshDao.Cell3DRemove(0));
    EXPECT_NO_THROW(meshDao.Cell3DRemove(0));
    EXPECT_NO_THROW(meshDao.Cell3DRemove(0));
    EXPECT_NO_THROW(meshDao.Compress());

    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 0);
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 0);
    EXPECT_EQ(meshDao.Cell2DTotalNumber(), 0);
    EXPECT_EQ(meshDao.Cell3DTotalNumber(), 0);
  }

  TEST(TestMesh, TestImportExportMesh2D)
  {
    MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    string exportFolder = "./Export";
    Gedim::Output::CreateFolder(exportFolder);
    exportFolder = exportFolder + "/TestImportExportMesh2D";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::MeshFromCsvUtilities meshFromCsvUtilities;

    Gedim::MeshFromCsvUtilities::Configuration exportConfiguration;
    exportConfiguration.Folder = exportFolder;
    Gedim::MeshDAOExporterToCsv exporter(meshFromCsvUtilities);
    EXPECT_NO_THROW(exporter.Export(exportConfiguration,
                                    meshDao));

    Gedim::MeshMatrices importedMesh;
    Gedim::MeshMatricesDAO importedMeshDao(importedMesh);

    Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
    meshImporterConfiguration.Folder = exportFolder;
    Gedim::MeshDAOImporterFromCsv importer(meshFromCsvUtilities);

    EXPECT_NO_THROW(importer.Import(meshImporterConfiguration,
                                    importedMeshDao));
    ASSERT_EQ(mesh.Mesh.Dimension, importedMesh.Dimension);
    ASSERT_EQ(mesh.Mesh.NumberCell0D, importedMesh.NumberCell0D);
    ASSERT_EQ(mesh.Mesh.Cell0DCoordinates, importedMesh.Cell0DCoordinates);
    ASSERT_EQ(mesh.Mesh.Cell0DMarkers, importedMesh.Cell0DMarkers);
    ASSERT_EQ(mesh.Mesh.ActiveCell0D, importedMesh.ActiveCell0D);
    ASSERT_EQ(mesh.Mesh.Cell0DDoublePropertyIds, importedMesh.Cell0DDoublePropertyIds);
    ASSERT_EQ(mesh.Mesh.NumberCell1D, importedMesh.NumberCell1D);
    ASSERT_EQ(mesh.Mesh.Cell1DMarkers, importedMesh.Cell1DMarkers);
    ASSERT_EQ(mesh.Mesh.ActiveCell1D, importedMesh.ActiveCell1D);
    ASSERT_EQ(mesh.Mesh.Cell1DDoublePropertyIds, importedMesh.Cell1DDoublePropertyIds);
    ASSERT_EQ(mesh.Mesh.NumberCell2D, importedMesh.NumberCell2D);
    ASSERT_EQ(mesh.Mesh.NumberCell2DVertices, importedMesh.NumberCell2DVertices);
    ASSERT_EQ(mesh.Mesh.Cell2DVertices, importedMesh.Cell2DVertices);
    ASSERT_EQ(mesh.Mesh.NumberCell2DEdges, importedMesh.NumberCell2DEdges);
    ASSERT_EQ(mesh.Mesh.Cell2DEdges, importedMesh.Cell2DEdges);
    ASSERT_EQ(mesh.Mesh.NumberCell2DSubdivision, importedMesh.NumberCell2DSubdivision);
    ASSERT_EQ(mesh.Mesh.Cell2DSubdivision, importedMesh.Cell2DSubdivision);
    ASSERT_EQ(mesh.Mesh.Cell2DMarkers, importedMesh.Cell2DMarkers);
    ASSERT_EQ(mesh.Mesh.ActiveCell2D, importedMesh.ActiveCell2D);
    ASSERT_EQ(mesh.Mesh.Cell2DDoublePropertyIds, importedMesh.Cell2DDoublePropertyIds);
    ASSERT_EQ(mesh.Mesh.NumberCell3D, importedMesh.NumberCell3D);
    ASSERT_EQ(mesh.Mesh.Cell3DMarkers, importedMesh.Cell3DMarkers);
    ASSERT_EQ(mesh.Mesh.ActiveCell3D, importedMesh.ActiveCell3D);
    ASSERT_EQ(mesh.Mesh.Cell3DDoublePropertyIds, importedMesh.Cell3DDoublePropertyIds);
  }

  TEST(TestMesh, MeshDAOImporter2DFromCsv)
  {
    MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    string exportFolder = "./Export";
    Gedim::Output::CreateFolder(exportFolder);
    exportFolder = exportFolder + "/MeshDAOImporter2DFromCsv";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::MeshFromCsvUtilities meshFromCsvUtilities;

    Gedim::MeshFromCsvUtilities::Configuration exportConfiguration;
    exportConfiguration.Folder = exportFolder;
    Gedim::MeshDAOExporterToCsv exporter(meshFromCsvUtilities);
    EXPECT_NO_THROW(exporter.Export(exportConfiguration,
                                    meshDao));

    Gedim::MeshMatrices importedMesh;
    Gedim::MeshMatricesDAO importedMeshDao(importedMesh);

    Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
    meshImporterConfiguration.Folder = exportFolder;
    Gedim::MeshDAOImporterFromCsv importer(meshFromCsvUtilities);

    EXPECT_NO_THROW(importer.ImportMesh2D(meshImporterConfiguration,
                                          importedMeshDao));

    ASSERT_EQ(mesh.Mesh.Dimension, importedMesh.Dimension);
    ASSERT_EQ(mesh.Mesh.NumberCell0D, importedMesh.NumberCell0D);
    ASSERT_EQ(mesh.Mesh.Cell0DCoordinates, importedMesh.Cell0DCoordinates);
    ASSERT_EQ(mesh.Mesh.Cell0DMarkers, importedMesh.Cell0DMarkers);
    ASSERT_EQ(mesh.Mesh.ActiveCell0D, importedMesh.ActiveCell0D);
    ASSERT_EQ(mesh.Mesh.Cell0DDoublePropertyIds, importedMesh.Cell0DDoublePropertyIds);
    ASSERT_EQ(mesh.Mesh.NumberCell1D, importedMesh.NumberCell1D);
    ASSERT_EQ(mesh.Mesh.Cell1DMarkers, importedMesh.Cell1DMarkers);
    ASSERT_EQ(mesh.Mesh.ActiveCell1D, importedMesh.ActiveCell1D);
    ASSERT_EQ(mesh.Mesh.Cell1DDoublePropertyIds, importedMesh.Cell1DDoublePropertyIds);
    ASSERT_EQ(mesh.Mesh.NumberCell2D, importedMesh.NumberCell2D);
    ASSERT_EQ(mesh.Mesh.Cell2DMarkers, importedMesh.Cell2DMarkers);
    ASSERT_EQ(mesh.Mesh.ActiveCell2D, importedMesh.ActiveCell2D);
    ASSERT_EQ(mesh.Mesh.Cell2DDoublePropertyIds, importedMesh.Cell2DDoublePropertyIds);
    ASSERT_EQ(mesh.Mesh.NumberCell3D, importedMesh.NumberCell3D);
    ASSERT_EQ(mesh.Mesh.Cell3DMarkers, importedMesh.Cell3DMarkers);
    ASSERT_EQ(mesh.Mesh.ActiveCell3D, importedMesh.ActiveCell3D);
    ASSERT_EQ(mesh.Mesh.Cell3DDoublePropertyIds, importedMesh.Cell3DDoublePropertyIds);
  }
}

#endif // __TEST_MESH_H
