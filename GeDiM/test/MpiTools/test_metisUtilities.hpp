#ifndef __TEST_METISUTILITIES_H
#define __TEST_METISUTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "FileTextReader.hpp"
#include "GeometryUtilities.hpp"
#include "MeshDAOImporterFromCsv.hpp"
#include "MeshUtilities.hpp"
#include "Macro.hpp"
#include "MeshMatricesDAO.hpp"
#include "MetisUtilities.hpp"
#include "VTKUtilities.hpp"
#include "IOUtilities.hpp"

#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "MeshMatrices_2D_4Cells_Mock.hpp"

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

#if ENABLE_METIS == 1
    ASSERT_EQ(std::vector<unsigned int>({ 1, 1, 0, 1, 0, 0 }), partition);
#else
    ASSERT_EQ(std::vector<unsigned int>({ 0, 0, 0, 0, 0, 0 }), partition);
#endif
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
    partitionOptions.NumberOfParts = 3;

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

  TEST(TestMetisUtilities, TestNetworkPartition_Mesh2D_DualGraph_Constraints)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    const Gedim::MeshUtilities::MeshGeometricData2D geometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                          meshDAO);

    std::vector<bool> edgeConstrained = std::vector<bool>(meshDAO.Cell1DTotalNumber(), false);
    edgeConstrained[2] = true;
    edgeConstrained[5] = true;

    const Gedim::MetisUtilities::Network network = metisUtilities.Mesh2DToDualGraph(meshDAO,
                                                                                    edgeConstrained);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partition = metisUtilities.NetworkPartition(partitionOptions,
                                                                                network);


    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_DualGraph_Weights";
    Gedim::Output::CreateFolder(exportFolder);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh");

    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell2DsCentroids[c];

      const Eigen::MatrixXi graphEdges = metisUtilities.GraphToConnectivityMatrix(network);

      std::vector<double> weights;
      weights.reserve(network.EdgeWeights.size());
      weights.assign(network.EdgeWeights.begin(),
                     network.EdgeWeights.end());

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
                             },
                             {
                               "Weights",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(weights.size()),
                               weights.data()
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

      exporter.Export(exportFolder + "/Partition.vtu");
    }
  }

  TEST(TestMetisUtilities, TestNetworkPartition_ImportMesh2D_DualGraph)
  {
    GTEST_SKIP();

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_ImportMesh2D_DualGraph";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::MetisUtilities metisUtilities;

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh);

    Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
    meshImporterConfiguration.Folder = "/home/geoscore/Dropbox/Polito/Articles/GE_MESH/Poj_MeshQuality_ToGe/NumericalTests/METIS/SingleDomain/M1";
    Gedim::MeshFromCsvUtilities meshFromCsvUtilities;
    Gedim::MeshDAOImporterFromCsv importer(meshFromCsvUtilities);
    importer.Import(meshImporterConfiguration,
                    meshDAO);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "ImportedMesh");


    const Gedim::MeshUtilities::MeshGeometricData2D geometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                          meshDAO);

    Eigen::SparseMatrix<unsigned int> weights(meshDAO.Cell2DTotalNumber(),
                                              meshDAO.Cell2DTotalNumber());

    {
      list<Eigen::Triplet<unsigned int>> triplets;

      Gedim::FileReader fileReader("/home/geoscore/Dropbox/Polito/Articles/GE_MESH/Poj_MeshQuality_ToGe/NumericalTests/METIS/SingleDomain/M1/Weights.txt");
      fileReader.Open();

      std::vector<string> lines;
      fileReader.GetAllLines(lines);
      fileReader.Close();

      Gedim::Output::Assert(lines.size() == meshDAO.Cell2DTotalNumber());

      double weight = 0.0;
      for (unsigned int l = 0; l < lines.size(); l++)
      {
        istringstream converter(lines[l]);

        for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
        {
          converter >> weight;

          if (weight > 0.0 && weight < meshDAO.Cell2DTotalNumber())
          {
            triplets.push_back(Eigen::Triplet<unsigned int>(l,
                                                            c,
                                                            round(meshDAO.Cell2DTotalNumber() *
                                                                  meshDAO.Cell2DTotalNumber() /
                                                                  weight)));
          }
        }
      }

      weights.setFromTriplets(triplets.begin(), triplets.end());
      weights.makeCompressed();
    }

    std::vector<bool> edgeIsConstrained = std::vector<bool>(meshDAO.Cell1DTotalNumber(),
                                                            false);
    {
      const unsigned int constrainIndex = meshDAO.Cell1DDoublePropertyIndex("marked");
      for (unsigned int e = 0; e < meshDAO.Cell1DTotalNumber(); e++)
        edgeIsConstrained[e] = meshDAO.Cell1DDoublePropertyValue(e, constrainIndex, 0) == 1.0;
    }

    const Gedim::MetisUtilities::Network network = metisUtilities.Mesh2DToDualGraph(meshDAO,
                                                                                    edgeIsConstrained,
                                                                                    weights);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 50;

    const std::vector<unsigned int> partition = metisUtilities.NetworkPartition(partitionOptions,
                                                                                network);


    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell2DsCentroids[c];

      const Eigen::MatrixXi graphEdges = metisUtilities.GraphToConnectivityMatrix(network);

      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(partition.size());
      property.assign(partition.begin(), partition.end());

      std::vector<double> weights;
      weights.reserve(network.EdgeWeights.size());
      weights.assign(network.EdgeWeights.begin(),
                     network.EdgeWeights.end());

      exporter.AddSegments(graphVertices,
                           graphEdges,
                           {
                             {
                               "Partition",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(property.size()),
                               property.data()
                             },
                             {
                               "Weights",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(weights.size()),
                               weights.data()
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

      exporter.Export(exportFolder + "/Partitioned.vtu");
    }
  }
}

#endif // __TEST_METISUTILITIES_H
