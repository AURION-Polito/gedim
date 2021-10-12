#ifndef __TEST_ExportConformedMesh_H
#define __TEST_ExportConformedMesh_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "IntersectorMesh2DSegment.hpp"
#include "UnionMeshSegment.hpp"
#include "ConformerMeshSegment.hpp"
#include "ConformerMeshPolygon.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "MeshMatrices_2D_4Cells_Mock.hpp"
#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "MeshDAOExporterToCsv.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
  TEST(TestExportConformedMesh, TestExportConformedMesh2D)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // conform of simple 2 points mesh with one segment
      {
        GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockOriginalMesh;
        Gedim::MeshMatricesDAO domainMesh(mockOriginalMesh.Mesh);

        Vector3d segmentOrigin(0.25, 0.25, 0.0);
        Vector3d segmentEnd(0.35, 0.35, 0.0);

        // Intersect mesh 2D with segment
        Gedim::IntersectorMesh2DSegment intersectorMeshSegment(domainMesh,
                                                               geometryUtilities);

        Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMesh;
        intersectorMeshSegment.CreateIntersectionMesh(segmentOrigin,
                                                      segmentEnd,
                                                      intersectionMesh);

        // Create mesh union
        vector<double> curvilinearCoordinatesMesh;
        Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMesh,
                                                                  curvilinearCoordinatesMesh);
        Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

        Gedim::UnionMeshSegment::UnionMesh unionMesh;
        unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMesh,
                                         curvilinearCoordinatesMesh,
                                         unionMesh);

        // Conform the mesh 1D
        Gedim::ConformerMeshSegment conformMeshSegment(geometryUtilities);

        Gedim::ConformerMeshSegment::ConformMesh conformMesh;
        ASSERT_NO_THROW(conformMeshSegment.CreateConformMesh(intersectionMesh,
                                                             unionMesh,
                                                             0,
                                                             conformMesh));

        // Conform the mesh 2D
        Gedim::ConformerMeshPolygon conformMeshPolygon(geometryUtilities);
        Gedim::ConformerMeshPolygon::ConformMesh domainConformedMeshData;

        ASSERT_NO_THROW(conformMeshPolygon.CreateConformMesh(segmentOrigin,
                                                             segmentEnd,
                                                             conformMesh,
                                                             domainMesh,
                                                             domainConformedMeshData));

        // Clean mesh 2D
        Gedim::MeshUtilities meshUtilities;
        Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
        ASSERT_NO_THROW(meshUtilities.ExtractActiveMesh(domainMesh,
                                                        extractionData));

        // Update mesh 1D with cleaned mesh 2D
        ASSERT_NO_THROW(conformMeshSegment.UpdateWithActiveMesh2D(extractionData,
                                                                  conformMesh));

        EXPECT_EQ(mockOriginalMesh.Mesh.NumberCell0D, 7);
        EXPECT_EQ(mockOriginalMesh.Mesh.NumberCell1D, 9);
        EXPECT_EQ(mockOriginalMesh.Mesh.NumberCell2D, 3);

        for (const auto& mesh1Dpoint : conformMesh.Points)
          EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
        for (const auto& mesh1Dsegment : conformMesh.Segments)
          EXPECT_EQ(mesh1Dsegment.Edge2DIds.size(), 1);

        // Add mesh 2D properties
        domainMesh.Cell0DInitializeDoubleProperties(1);
        domainMesh.Cell0DAddDoubleProperty("flag");
        for (unsigned int p = 0; p < domainMesh.Cell0DTotalNumber(); p++)
        {
          domainMesh.Cell0DInitializeDoublePropertyValues(p, 0, 1);
          domainMesh.Cell0DInsertDoublePropertyValue(p, 0, 0, 0.0);
        }
        for (const pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>& point : conformMesh.Points)
        {
          for (const unsigned int& v : point.second.Vertex2DIds)
            domainMesh.Cell0DInsertDoublePropertyValue(v, 0, 0, 1.0);
        }

        domainMesh.Cell1DInitializeDoubleProperties(1);
        domainMesh.Cell1DAddDoubleProperty("flag");
        for (unsigned int e = 0; e < domainMesh.Cell1DTotalNumber(); e++)
        {
          domainMesh.Cell1DInitializeDoublePropertyValues(e, 0, 1);
          domainMesh.Cell1DInsertDoublePropertyValue(e, 0, 0, 0.0);
        }
        for (const Gedim::ConformerMeshSegment::ConformMesh::ConformMeshSegment& segment : conformMesh.Segments)
        {
          for (const unsigned int& e : segment.Edge2DIds)
            domainMesh.Cell1DInsertDoublePropertyValue(e, 0, 0, 1.0);
        }

        // Export the resulting mesh
        string exportFolder = "./Export";
        Gedim::Output::CreateFolder(exportFolder);
        exportFolder = exportFolder + "/TestExportConformedMesh2D";
        Gedim::Output::CreateFolder(exportFolder);

        Gedim::MeshFromCsvUtilities meshFromCsvUtilities;
        Gedim::MeshFromCsvUtilities::Configuration exportConfiguration;
        exportConfiguration.Folder = exportFolder;
        Gedim::MeshDAOExporterToCsv exporter(meshFromCsvUtilities);
        EXPECT_NO_THROW(exporter.Export(exportConfiguration,
                                        domainMesh));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_ExportConformedMesh_H
