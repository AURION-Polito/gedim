#ifndef __TEST_METISUTILITIES_H
#define __TEST_METISUTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "MeshUtilities.hpp"
#include "Macro.hpp"
#include "MeshMatricesDAO.hpp"
#include "MetisUtilities.hpp"
#include "VTKUtilities.hpp"
#include "IOUtilities.hpp"

#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "MeshMatrices_2D_26Cells_Mock.hpp"

namespace UnitTesting
{
  TEST(TestMetisUtilities, TestNetworkPartition)
  {
    Gedim::MetisUtilities metisUtilities;

    const unsigned int numberVertices = 6;
    Eigen::MatrixXi edges(2, 6);
    edges.col(0)<< 0, 1;
    edges.col(1)<< 1, 2;
    edges.col(2)<< 1, 3;
    edges.col(3)<< 1, 5;
    edges.col(4)<< 2, 4;
    edges.col(5)<< 2, 5;

    const Gedim::MetisUtilities::Network network = metisUtilities.Mesh2DToGraph(numberVertices,
                                                                                edges,
                                                                                true);

    ASSERT_EQ(std::vector<unsigned int>({ 0,1,5,8,9,10,12 }), network.AdjacencyRows);
    ASSERT_EQ(std::vector<unsigned int>({ 1,0,2,3,5,1,4,5,1,2,1,2 }), network.AdjacencyCols);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partition = metisUtilities.NetworkPartition(partitionOptions,
                                                                                network);

    ASSERT_EQ(std::vector<unsigned int>({ 1, 1, 0, 1, 0, 0 }), partition);
  }

  TEST(TestMetisUtilities, TestNetworkPartition_Mesh2D_Graph)
  {
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    const Eigen::MatrixXi edges = meshDAO.Cell1DsExtremes();

    const Gedim::MetisUtilities::Network network = metisUtilities.Mesh2DToGraph(meshDAO.Cell0DTotalNumber(),
                                                                                edges,
                                                                                true);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partition = metisUtilities.NetworkPartition(partitionOptions,
                                                                                network);


    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_Graph";
    Gedim::Output::CreateFolder(exportFolder);

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(partition.size());
      property.assign(partition.begin(), partition.end());

      exporter.AddSegments(meshDAO.Cell0DsCoordinates(),
                           edges,
                           {
                             {
                               "Partition",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(property.size()),
                               property.data()
                             }
                           });

      exporter.Export(exportFolder + "/Cell1Ds.vtu");
    }
  }

  TEST(TestMetisUtilities, TestNetworkPartition_Mesh2D_DualGraph)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    const Gedim::MeshUtilities::MeshGeometricData2D geometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                          meshDAO);

    const Gedim::MetisUtilities::Network network = metisUtilities.Mesh2DToDualGraph(meshDAO);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partition = metisUtilities.NetworkPartition(partitionOptions,
                                                                                network);


    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_DualGraph";
    Gedim::Output::CreateFolder(exportFolder);

    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell2DsCentroids[c];

      const Eigen::MatrixXi graphEdges = metisUtilities.GraphToConnectivityMatrix(network);

      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(partition.size());
      property.assign(partition.begin(), partition.end());

      exporter.AddSegments(graphVertices,
                           graphEdges,
                           {
                             {
                               "Partition",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(property.size()),
                               property.data()
                             }
                           });

      exporter.Export(exportFolder + "/Graph.vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(partition.size());
      property.assign(partition.begin(), partition.end());

      exporter.AddPolygons(meshDAO.Cell0DsCoordinates(),
                           meshDAO.Cell2DsVertices(),
                           {
                             {
                               "Partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(property.size()),
                               property.data()
                             }
                           });

      exporter.Export(exportFolder + "/Cell2Ds.vtu");
    }
  }
}

#endif // __TEST_METISUTILITIES_H
