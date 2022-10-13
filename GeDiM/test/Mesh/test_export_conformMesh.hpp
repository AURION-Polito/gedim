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
#include "MashMatrices_2D_CleanTest_Mock.hpp"
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

        const Vector3d segmentOrigin(0.25, 0.25, 0.0);
        const Vector3d segmentEnd(0.35, 0.35, 0.0);
        const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin,
                                                                                segmentEnd);
        const Eigen::Vector3d segmentBarycenter = geometryUtilities.SegmentBarycenter(segmentOrigin,
                                                                                      segmentEnd);
        const double segmentLength = geometryUtilities.SegmentLength(segmentOrigin,
                                                                     segmentEnd);
        const double segmentSquaredLength = segmentLength * segmentLength;

        // Intersect mesh 2D with segment
        Gedim::IntersectorMesh2DSegment intersectorMeshSegment(domainMesh,
                                                               geometryUtilities);

        Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMesh;
        intersectorMeshSegment.CreateIntersectionMesh(segmentOrigin,
                                                      segmentEnd,
                                                      segmentTangent,
                                                      segmentBarycenter,
                                                      segmentLength,
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
                                                             segmentTangent,
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

        // Check missing mesh 2D cell0Ds
        ASSERT_NO_THROW(conformMeshSegment.AddMissingMesh2DCell0Ds(segmentOrigin,
                                                                   segmentTangent,
                                                                   segmentSquaredLength,
                                                                   domainMesh,
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

  TEST(TestExportConformedMesh, TestConformComplexMesh2D)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      {
        GedimUnitTesting::MashMatrices_2D_CleanTest_Mock mockOriginalMesh;
        Gedim::MeshMatricesDAO domainMesh(mockOriginalMesh.Mesh);

        const Eigen::Vector3d segmentOrigin(1.7745275237876366e+00, 2.5306770050657718e-01, 0.0000000000000000e+00);
        const Eigen::Vector3d segmentEnd(3.5355337963926114e+00, 1.7400074306985638e+00, 0.0000000000000000e+00);
        const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin,
                                                                                segmentEnd);

        Gedim::ConformerMeshSegment::ConformMesh segmentMesh;

        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(0.0000000000000000e+00, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(2.7221089123105585e-02, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(1.5081657587860950e-01, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(2.0636899211510995e-01, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(2.2032877240731363e-01, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(2.7084084009418125e-01, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(2.8618524662562461e-01, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(2.9736416870129134e-01, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(3.5907750767017105e-01, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(4.0661516876623027e-01, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(4.3798874569320156e-01, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(6.6913757867999368e-01, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        segmentMesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(1.0000000000000000e+00, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));

        segmentMesh.Points.at(0.0000000000000000e+00).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original;
        segmentMesh.Points.at(2.7221089123105585e-02).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        segmentMesh.Points.at(1.5081657587860950e-01).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        segmentMesh.Points.at(2.0636899211510995e-01).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        segmentMesh.Points.at(2.2032877240731363e-01).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        segmentMesh.Points.at(2.7084084009418125e-01).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        segmentMesh.Points.at(2.8618524662562461e-01).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        segmentMesh.Points.at(2.9736416870129134e-01).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        segmentMesh.Points.at(3.5907750767017105e-01).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original;
        segmentMesh.Points.at(4.0661516876623027e-01).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        segmentMesh.Points.at(4.3798874569320156e-01).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        segmentMesh.Points.at(6.6913757867999368e-01).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        segmentMesh.Points.at(1.0000000000000000e+00).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original;

        segmentMesh.Points.at(0.0000000000000000e+00).Cell2DIds.insert(3);
        segmentMesh.Points.at(2.7221089123105585e-02).Cell2DIds.insert(3);
        segmentMesh.Points.at(1.5081657587860950e-01).Cell2DIds.insert(3);
        segmentMesh.Points.at(2.0636899211510995e-01).Cell2DIds.insert(3);
        segmentMesh.Points.at(2.2032877240731363e-01).Cell2DIds.insert(3);
        segmentMesh.Points.at(2.7084084009418125e-01).Cell2DIds.insert(3);
        segmentMesh.Points.at(2.8618524662562461e-01).Cell2DIds.insert(3);
        segmentMesh.Points.at(2.9736416870129134e-01).Cell2DIds.insert(3);
        segmentMesh.Points.at(3.5907750767017105e-01).Vertex2DIds.insert(40);
        segmentMesh.Points.at(3.5907750767017105e-01).Edge2DIds.insert(43);
        segmentMesh.Points.at(3.5907750767017105e-01).Edge2DIds.insert(44);
        segmentMesh.Points.at(3.5907750767017105e-01).Cell2DIds.insert(3);
        segmentMesh.Points.at(3.5907750767017105e-01).Cell2DIds.insert(4);
        segmentMesh.Points.at(4.0661516876623027e-01).Cell2DIds.insert(4);
        segmentMesh.Points.at(4.3798874569320156e-01).Cell2DIds.insert(4);
        segmentMesh.Points.at(6.6913757867999368e-01).Cell2DIds.insert(4);
        segmentMesh.Points.at(1.0000000000000000e+00).Edge2DIds.insert(21);
        segmentMesh.Points.at(1.0000000000000000e+00).Cell2DIds.insert(4);

        segmentMesh.Segments.resize(12);
        segmentMesh.Segments[0].Points = vector<double> { 0.0000000000000000e+00, 2.7221089123105585e-02 };
        segmentMesh.Segments[1].Points = vector<double> { 2.7221089123105585e-02, 1.5081657587860950e-01 };
        segmentMesh.Segments[2].Points = vector<double> { 1.5081657587860950e-01, 2.0636899211510995e-01 };
        segmentMesh.Segments[3].Points = vector<double> { 2.0636899211510995e-01, 2.2032877240731363e-01 };
        segmentMesh.Segments[4].Points = vector<double> { 2.2032877240731363e-01, 2.7084084009418125e-01 };
        segmentMesh.Segments[5].Points = vector<double> { 2.7084084009418125e-01, 2.8618524662562461e-01 };
        segmentMesh.Segments[6].Points = vector<double> { 2.8618524662562461e-01, 2.9736416870129134e-01 };
        segmentMesh.Segments[7].Points = vector<double> { 2.9736416870129134e-01, 3.5907750767017105e-01 };
        segmentMesh.Segments[8].Points = vector<double> { 3.5907750767017105e-01, 4.0661516876623027e-01 };
        segmentMesh.Segments[9].Points = vector<double> { 4.0661516876623027e-01, 4.3798874569320156e-01 };
        segmentMesh.Segments[10].Points = vector<double> { 4.3798874569320156e-01, 6.6913757867999368e-01 };
        segmentMesh.Segments[11].Points = vector<double> { 6.6913757867999368e-01, 1.0000000000000000e+00 };

        segmentMesh.Segments[0].Cell2DIds.insert(3);
        segmentMesh.Segments[1].Cell2DIds.insert(3);
        segmentMesh.Segments[2].Cell2DIds.insert(3);
        segmentMesh.Segments[3].Cell2DIds.insert(3);
        segmentMesh.Segments[4].Cell2DIds.insert(3);
        segmentMesh.Segments[5].Cell2DIds.insert(3);
        segmentMesh.Segments[6].Cell2DIds.insert(3);
        segmentMesh.Segments[7].Cell2DIds.insert(3);
        segmentMesh.Segments[8].Cell2DIds.insert(4);
        segmentMesh.Segments[9].Cell2DIds.insert(4);
        segmentMesh.Segments[10].Cell2DIds.insert(4);
        segmentMesh.Segments[11].Cell2DIds.insert(4);

        Gedim::ConformerMeshPolygon::ConformMesh segmentConformMeshInfos;

        Gedim::ConformerMeshPolygon conformerMeshDomain(geometryUtilities);
        conformerMeshDomain.CreateConformMesh(segmentOrigin,
                                              segmentEnd,
                                              segmentTangent,
                                              segmentMesh,
                                              domainMesh,
                                              segmentConformMeshInfos);

        // Clean mesh 2D
        Gedim::MeshUtilities meshUtilities;
        Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
        meshUtilities.ExtractActiveMesh(domainMesh,
                                        extractionData);

        EXPECT_EQ(mockOriginalMesh.Mesh.NumberCell2D, 9);

        // Export the resulting mesh
        string exportFolder = "./Export";
        Gedim::Output::CreateFolder(exportFolder);
        exportFolder = exportFolder + "/TestCleanMesh2D";
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
