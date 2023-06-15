#ifndef __TEST_METISUTILITIES_H
#define __TEST_METISUTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "FileTextReader.hpp"
#include "GeometryUtilities.hpp"
#include "GraphUtilities.hpp"
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

#include "MeshMatrices_3D_22Cells_Mock.hpp"
#include "MeshMatrices_3D_329Cells_Mock.hpp"

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

    const Gedim::MetisUtilities::MetisNetwork network = metisUtilities.Mesh2DToGraph(numberVertices,
                                                                                     edges,
                                                                                     true);

    ASSERT_EQ(std::vector<unsigned int>({ 0,1,5,8,9,10,12 }), network.Adjacency.Rows);
    ASSERT_EQ(std::vector<unsigned int>({ 1,0,2,3,5,1,4,5,1,2,1,2 }), network.Adjacency.Cols);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions,
                                                                                 network);

    const std::vector<unsigned int> fix_constraints_partitions = metisUtilities.PartitionCheckConstraints(network,
                                                                                                          partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions = metisUtilities.PartitionCheckConnectedComponents(network,
                                                                                                                          fix_constraints_partitions);

#if ENABLE_METIS == 1
    ASSERT_EQ(std::vector<unsigned int>({ 1, 1, 0, 1, 0, 0 }), partitions);
#else
    ASSERT_EQ(std::vector<unsigned int>({ 0, 0, 0, 0, 0, 0 }), partitions);
#endif
  }

  TEST(TestMetisUtilities, TestNetworkPartition_Mesh2D_Graph)
  {
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    const Eigen::MatrixXi edges = meshDAO.Cell1DsExtremes();

    const Gedim::MetisUtilities::MetisNetwork network = metisUtilities.Mesh2DToGraph(meshDAO.Cell0DTotalNumber(),
                                                                                     edges,
                                                                                     true);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions,
                                                                                 network);

    const std::vector<unsigned int> fix_constraints_partitions = metisUtilities.PartitionCheckConstraints(network,
                                                                                                          partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions = metisUtilities.PartitionCheckConnectedComponents(network,
                                                                                                                          fix_constraints_partitions);


    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_Graph";
    Gedim::Output::CreateFolder(exportFolder);

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(fix_connectedComponents_partitions.size());
      property.assign(fix_connectedComponents_partitions.begin(),
                      fix_connectedComponents_partitions.end());

      exporter.AddSegments(meshDAO.Cell0DsCoordinates(),
                           edges,
                           {
                             {
                               "partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(property.size()),
                               property.data()
                             }
                           });

      exporter.Export(exportFolder + "/Mesh_Cell1Ds.vtu");
    }
  }

  TEST(TestMetisUtilities, TestNetworkPartition_Mesh2D_DualGraph)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    const Gedim::MeshUtilities::MeshGeometricData2D geometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                          meshDAO);

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh2DToDualGraph(meshDAO);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 3;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions,
                                                                                 meshToNetwork.Network);

    const std::vector<unsigned int> fix_constraints_partitions = metisUtilities.PartitionCheckConstraints(meshToNetwork.Network,
                                                                                                          partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions = metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network,
                                                                                                                          fix_constraints_partitions);


    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_DualGraph";
    Gedim::Output::CreateFolder(exportFolder);

    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell2DsCentroids[c];

      const unsigned int graphNumEdges = meshToNetwork.Network.Adjacency.Cols.size();
      const std::vector<std::vector<unsigned int>> graphAdjacency = metisUtilities.MetisAdjacencyToGraphAdjacency(meshToNetwork.Network.Adjacency);
      const Eigen::MatrixXi graphEdges = graphUtilities.GraphAdjacencyToGraphConnectivity(graphNumEdges,
                                                                                          graphAdjacency);

      {
        Gedim::VTKUtilities exporter;

        std::vector<double> partition;
        partition.reserve(partitions.size());
        partition.assign(partitions.begin(), partitions.end());
        std::vector<double> fix_constraints_partition;
        fix_constraints_partition.reserve(fix_constraints_partitions.size());
        fix_constraints_partition.assign(fix_constraints_partitions.begin(), fix_constraints_partitions.end());
        std::vector<double> fix_connectedComponents_partition;
        fix_connectedComponents_partition.reserve(fix_connectedComponents_partitions.size());
        fix_connectedComponents_partition.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());

        exporter.AddPoints(graphVertices,
                           {
                             {
                               "partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(partition.size()),
                               partition.data()
                             },
                             {
                               "fix_constraints_partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(fix_constraints_partition.size()),
                               fix_constraints_partition.data()
                             },
                             {
                               "fix_connectedComponents_partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(fix_connectedComponents_partition.size()),
                               fix_connectedComponents_partition.data()
                             }
                           });

        exporter.Export(exportFolder + "/Graph_Points.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        std::vector<double> property;
        property.reserve(meshToNetwork.Network.EdgesWeight.size());
        property.assign(meshToNetwork.Network.EdgesWeight.begin(),
                        meshToNetwork.Network.EdgesWeight.end());

        exporter.AddSegments(graphVertices,
                             graphEdges,
                             {
                               {
                                 "Weight",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(property.size()),
                                 property.data()
                               }
                             });

        exporter.Export(exportFolder + "/Graph_Edges.vtu");
      }
    }

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(fix_connectedComponents_partitions.size());
      property.assign(fix_connectedComponents_partitions.begin(),
                      fix_connectedComponents_partitions.end());

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

      exporter.Export(exportFolder + "/Mesh_Cell2Ds.vtu");
    }
  }

  TEST(TestMetisUtilities, TestNetworkPartition_Mesh2D_DualGraph_Constraints)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    const Gedim::MeshUtilities::MeshGeometricData2D geometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                          meshDAO);

    std::vector<bool> cell1DsConstrained = std::vector<bool>(meshDAO.Cell1DTotalNumber(), false);
    cell1DsConstrained[2] = true;
    cell1DsConstrained[5] = true;

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh2DToDualGraph(meshDAO,
                                                                                                cell1DsConstrained);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions,
                                                                                 meshToNetwork.Network);

    const std::vector<unsigned int> fix_constraints_partitions = metisUtilities.PartitionCheckConstraints(meshToNetwork.Network,
                                                                                                          partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions = metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network,
                                                                                                                          fix_constraints_partitions);

    for (unsigned int e = 0; e < cell1DsConstrained.size(); e++)
    {
      if (!cell1DsConstrained[e] ||
          meshDAO.Cell1DNumberNeighbourCell2D(e) < 2)
        continue;

      ASSERT_NE(fix_connectedComponents_partitions.at(meshDAO.Cell1DNeighbourCell2D(e,
                                                                                    0)),
                fix_connectedComponents_partitions.at(meshDAO.Cell1DNeighbourCell2D(e,
                                                                                    1)));
    }

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_DualGraph_Constraints";
    Gedim::Output::CreateFolder(exportFolder);

    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell2DsCentroids[c];

      const unsigned int graphNumEdges = meshToNetwork.Network.Adjacency.Cols.size();
      const std::vector<std::vector<unsigned int>> graphAdjacency = metisUtilities.MetisAdjacencyToGraphAdjacency(meshToNetwork.Network.Adjacency);
      const Eigen::MatrixXi graphEdges = graphUtilities.GraphAdjacencyToGraphConnectivity(graphNumEdges,
                                                                                          graphAdjacency);
      {
        Gedim::VTKUtilities exporter;

        std::vector<double> partition;
        partition.reserve(partitions.size());
        partition.assign(partitions.begin(), partitions.end());
        std::vector<double> fix_constraints_partition;
        fix_constraints_partition.reserve(fix_constraints_partitions.size());
        fix_constraints_partition.assign(fix_constraints_partitions.begin(), fix_constraints_partitions.end());
        std::vector<double> fix_connectedComponents_partition;
        fix_connectedComponents_partition.reserve(fix_connectedComponents_partitions.size());
        fix_connectedComponents_partition.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());


        exporter.AddPoints(graphVertices,
                           {
                             {
                               "partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(partition.size()),
                               partition.data()
                             },
                             {
                               "fix_constraints_partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(fix_constraints_partition.size()),
                               fix_constraints_partition.data()
                             },
                             {
                               "fix_connectedComponents_partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(fix_connectedComponents_partition.size()),
                               fix_connectedComponents_partition.data()
                             }
                           });

        exporter.Export(exportFolder + "/Graph_Points.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        std::vector<double> property;
        property.reserve(meshToNetwork.Network.EdgesWeight.size());
        property.assign(meshToNetwork.Network.EdgesWeight.begin(),
                        meshToNetwork.Network.EdgesWeight.end());

        exporter.AddSegments(graphVertices,
                             graphEdges,
                             {
                               {
                                 "Weight",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(property.size()),
                                 property.data()
                               }
                             });

        exporter.Export(exportFolder + "/Graph_Edges.vtu");
      }
    }

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(fix_connectedComponents_partitions.size());
      property.assign(fix_connectedComponents_partitions.begin(),
                      fix_connectedComponents_partitions.end());

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

      exporter.Export(exportFolder + "/Mesh_Cell2Ds.vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> index, constrained, weight;
      index.resize(meshDAO.Cell1DTotalNumber());
      constrained.resize(meshDAO.Cell1DTotalNumber());
      weight.resize(meshDAO.Cell1DTotalNumber(), 0.0);

      for (unsigned int e = 0; e < meshDAO.Cell1DTotalNumber(); e++)
      {
        index[e] = e;
        constrained[e] = cell1DsConstrained[e] ? 1.0 : 0.0;
      }

      for (unsigned int e = 0; e < meshToNetwork.Network.EdgesWeight.size(); e++)
        weight[meshToNetwork.EdgesMeshCellIndex[e]] = meshToNetwork.Network.EdgesWeight[e];

      exporter.AddSegments(meshDAO.Cell0DsCoordinates(),
                           meshDAO.Cell1DsExtremes(),
                           {
                             {
                               "Index",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(index.size()),
                               index.data()
                             },
                             {
                               "Constrained",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(constrained.size()),
                               constrained.data()
                             },
                             {
                               "Weight",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(weight.size()),
                               weight.data()
                             }
                           });

      exporter.Export(exportFolder + "/Mesh_Cell1Ds.vtu");
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

    Gedim::GraphUtilities graphUtilities;
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

    std::vector<bool> cell1DsConstrained = std::vector<bool>(meshDAO.Cell1DTotalNumber(),
                                                             false);
    {
      const unsigned int constrainIndex = meshDAO.Cell1DDoublePropertyIndex("marked");
      for (unsigned int e = 0; e < meshDAO.Cell1DTotalNumber(); e++)
        cell1DsConstrained[e] = meshDAO.Cell1DDoublePropertyValue(e, constrainIndex, 0) == 1.0;
    }

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh2DToDualGraph(meshDAO,
                                                                                                cell1DsConstrained,
                                                                                                weights);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 50;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions,
                                                                                 meshToNetwork.Network);

    const std::vector<unsigned int> fix_constraints_partitions = metisUtilities.PartitionCheckConstraints(meshToNetwork.Network,
                                                                                                          partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions = metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network,
                                                                                                                          fix_constraints_partitions);


    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell2DsCentroids[c];

      const unsigned int graphNumEdges = meshToNetwork.Network.Adjacency.Cols.size();
      const std::vector<std::vector<unsigned int>> graphAdjacency = metisUtilities.MetisAdjacencyToGraphAdjacency(meshToNetwork.Network.Adjacency);
      const Eigen::MatrixXi graphEdges = graphUtilities.GraphAdjacencyToGraphConnectivity(graphNumEdges,
                                                                                          graphAdjacency);

      {
        Gedim::VTKUtilities exporter;

        std::vector<double> partition;
        partition.reserve(partitions.size());
        partition.assign(partitions.begin(), partitions.end());
        std::vector<double> fix_constraints_partition;
        fix_constraints_partition.reserve(fix_constraints_partitions.size());
        fix_constraints_partition.assign(fix_constraints_partitions.begin(), fix_constraints_partitions.end());
        std::vector<double> fix_connectedComponents_partition;
        fix_connectedComponents_partition.reserve(fix_connectedComponents_partitions.size());
        fix_connectedComponents_partition.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());

        exporter.AddPoints(graphVertices,
                           {
                             {
                               "partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(partition.size()),
                               partition.data()
                             },
                             {
                               "fix_constraints_partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(fix_constraints_partition.size()),
                               fix_constraints_partition.data()
                             },
                             {
                               "fix_connectedComponents_partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(fix_connectedComponents_partition.size()),
                               fix_connectedComponents_partition.data()
                             }
                           });

        exporter.Export(exportFolder + "/Graph_Points.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        std::vector<double> property;
        property.reserve(meshToNetwork.Network.EdgesWeight.size());
        property.assign(meshToNetwork.Network.EdgesWeight.begin(),
                        meshToNetwork.Network.EdgesWeight.end());

        exporter.AddSegments(graphVertices,
                             graphEdges,
                             {
                               {
                                 "Weight",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(property.size()),
                                 property.data()
                               }
                             });

        exporter.Export(exportFolder + "/Graph_Edges.vtu");
      }
    }

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> index, constrained, weight;
      index.resize(meshDAO.Cell1DTotalNumber());
      constrained.resize(meshDAO.Cell1DTotalNumber());
      weight.resize(meshDAO.Cell1DTotalNumber(), 0.0);

      for (unsigned int e = 0; e < meshDAO.Cell1DTotalNumber(); e++)
      {
        index[e] = e;
        constrained[e] = cell1DsConstrained[e] ? 1.0 : 0.0;
      }

      for (unsigned int e = 0; e < meshToNetwork.Network.EdgesWeight.size(); e++)
        weight[meshToNetwork.EdgesMeshCellIndex[e]] = meshToNetwork.Network.EdgesWeight[e];

      exporter.AddSegments(meshDAO.Cell0DsCoordinates(),
                           meshDAO.Cell1DsExtremes(),
                           {
                             {
                               "Index",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(index.size()),
                               index.data()
                             },
                             {
                               "Constrained",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(constrained.size()),
                               constrained.data()
                             },
                             {
                               "Weight",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(weight.size()),
                               weight.data()
                             }
                           });

      exporter.Export(exportFolder + "/Mesh_Cell1Ds.vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(fix_connectedComponents_partitions.size());
      property.assign(fix_connectedComponents_partitions.begin(),
                      fix_connectedComponents_partitions.end());

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

      exporter.Export(exportFolder + "/Mesh_Cell2Ds.vtu");
    }

    for (unsigned int e = 0; e < cell1DsConstrained.size(); e++)
    {
      if (!cell1DsConstrained[e] ||
          meshDAO.Cell1DNumberNeighbourCell2D(e) < 2)
        continue;

      ASSERT_NE(fix_connectedComponents_partitions.at(meshDAO.Cell1DNeighbourCell2D(e,
                                                                                    0)),
                fix_connectedComponents_partitions.at(meshDAO.Cell1DNeighbourCell2D(e,
                                                                                    1)));
    }
  }

  TEST(TestMetisUtilities, TestNetworkPartition_Mesh3D_DualGraph)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_3D_22Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    meshUtilities.ComputeCell2DCell3DNeighbours(meshDAO);

    const Gedim::MeshUtilities::MeshGeometricData3D geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                          meshDAO);

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh3DToDualGraph(meshDAO);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 3;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions,
                                                                                 meshToNetwork.Network);

    const std::vector<unsigned int> fix_constraints_partitions = metisUtilities.PartitionCheckConstraints(meshToNetwork.Network,
                                                                                                          partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions = metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network,
                                                                                                                          fix_constraints_partitions);


    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh3D_DualGraph";
    Gedim::Output::CreateFolder(exportFolder);

    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell3DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell3DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell3DsCentroids[c];

      const unsigned int graphNumEdges = meshToNetwork.Network.Adjacency.Cols.size();
      const std::vector<std::vector<unsigned int>> graphAdjacency = metisUtilities.MetisAdjacencyToGraphAdjacency(meshToNetwork.Network.Adjacency);
      const Eigen::MatrixXi graphEdges = graphUtilities.GraphAdjacencyToGraphConnectivity(graphNumEdges,
                                                                                          graphAdjacency);

      {
        Gedim::VTKUtilities exporter;

        std::vector<double> partition;
        partition.reserve(partitions.size());
        partition.assign(partitions.begin(), partitions.end());
        std::vector<double> fix_constraints_partition;
        fix_constraints_partition.reserve(fix_constraints_partitions.size());
        fix_constraints_partition.assign(fix_constraints_partitions.begin(), fix_constraints_partitions.end());
        std::vector<double> fix_connectedComponents_partition;
        fix_connectedComponents_partition.reserve(fix_connectedComponents_partitions.size());
        fix_connectedComponents_partition.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());

        exporter.AddPoints(graphVertices,
                           {
                             {
                               "partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(partition.size()),
                               partition.data()
                             },
                             {
                               "fix_constraints_partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(fix_constraints_partition.size()),
                               fix_constraints_partition.data()
                             },
                             {
                               "fix_connectedComponents_partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(fix_connectedComponents_partition.size()),
                               fix_connectedComponents_partition.data()
                             }
                           });

        exporter.Export(exportFolder + "/Graph_Points.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        std::vector<double> property;
        property.reserve(meshToNetwork.Network.EdgesWeight.size());
        property.assign(meshToNetwork.Network.EdgesWeight.begin(),
                        meshToNetwork.Network.EdgesWeight.end());

        exporter.AddSegments(graphVertices,
                             graphEdges,
                             {
                               {
                                 "Weight",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(property.size()),
                                 property.data()
                               }
                             });

        exporter.Export(exportFolder + "/Graph_Edges.vtu");
      }
    }

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(fix_connectedComponents_partitions.size());
      property.assign(fix_connectedComponents_partitions.begin(),
                      fix_connectedComponents_partitions.end());

      exporter.AddPolyhedrons(meshDAO.Cell0DsCoordinates(),
                              meshDAO.Cell3DsFacesVertices(),
                              {
                                {
                                  "Partition",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(property.size()),
                                  property.data()
                                }
                              });

      exporter.Export(exportFolder + "/Mesh_Cell3Ds.vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> index;
      index.resize(meshDAO.Cell2DTotalNumber());
      for (unsigned int f = 0; f < meshDAO.Cell2DTotalNumber(); f++)
        index[f] = f;

      exporter.AddPolygons(meshDAO.Cell0DsCoordinates(),
                           meshDAO.Cell2DsVertices(),
                           {
                             {
                               "Index",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(index.size()),
                               index.data()
                             }
                           });

      exporter.Export(exportFolder + "/Mesh_Cell2Ds.vtu");
    }
  }

  TEST(TestMetisUtilities, TestNetworkPartition_Mesh3D_DualGraph_Constraints)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_3D_22Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    meshUtilities.ComputeCell2DCell3DNeighbours(meshDAO);

    const Gedim::MeshUtilities::MeshGeometricData3D geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                          meshDAO);

    std::vector<bool> cell2DsConstrained = std::vector<bool>(meshDAO.Cell2DTotalNumber(),
                                                             false);
    cell2DsConstrained[2] = true;
    cell2DsConstrained[23] = true;
    cell2DsConstrained[29] = true;

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh3DToDualGraph(meshDAO,
                                                                                                cell2DsConstrained);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 3;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions,
                                                                                 meshToNetwork.Network);

    const std::vector<unsigned int> fix_constraints_partitions = metisUtilities.PartitionCheckConstraints(meshToNetwork.Network,
                                                                                                          partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions = metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network,
                                                                                                                          fix_constraints_partitions);

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh3D_DualGraph_Constraints";
    Gedim::Output::CreateFolder(exportFolder);

    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell3DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell3DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell3DsCentroids[c];

      const unsigned int graphNumEdges = meshToNetwork.Network.Adjacency.Cols.size();
      const std::vector<std::vector<unsigned int>> graphAdjacency = metisUtilities.MetisAdjacencyToGraphAdjacency(meshToNetwork.Network.Adjacency);
      const Eigen::MatrixXi graphEdges = graphUtilities.GraphAdjacencyToGraphConnectivity(graphNumEdges,
                                                                                          graphAdjacency);

      {
        Gedim::VTKUtilities exporter;

        std::vector<double> partition;
        partition.reserve(partitions.size());
        partition.assign(partitions.begin(), partitions.end());
        std::vector<double> fix_constraints_partition;
        fix_constraints_partition.reserve(fix_constraints_partitions.size());
        fix_constraints_partition.assign(fix_constraints_partitions.begin(), fix_constraints_partitions.end());
        std::vector<double> fix_connectedComponents_partition;
        fix_connectedComponents_partition.reserve(fix_connectedComponents_partitions.size());
        fix_connectedComponents_partition.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());

        exporter.AddPoints(graphVertices,
                           {
                             {
                               "partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(partition.size()),
                               partition.data()
                             },
                             {
                               "fix_constraints_partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(fix_constraints_partition.size()),
                               fix_constraints_partition.data()
                             },
                             {
                               "fix_connectedComponents_partition",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(fix_connectedComponents_partition.size()),
                               fix_connectedComponents_partition.data()
                             }
                           });

        exporter.Export(exportFolder + "/Graph_Points.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        std::vector<double> property;
        property.reserve(meshToNetwork.Network.EdgesWeight.size());
        property.assign(meshToNetwork.Network.EdgesWeight.begin(),
                        meshToNetwork.Network.EdgesWeight.end());

        exporter.AddSegments(graphVertices,
                             graphEdges,
                             {
                               {
                                 "Weight",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(property.size()),
                                 property.data()
                               }
                             });

        exporter.Export(exportFolder + "/Graph_Edges.vtu");
      }
    }

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(fix_connectedComponents_partitions.size());
      property.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());

      exporter.AddPolyhedrons(meshDAO.Cell0DsCoordinates(),
                              meshDAO.Cell3DsFacesVertices(),
                              {
                                {
                                  "Partition",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(property.size()),
                                  property.data()
                                }
                              });

      exporter.Export(exportFolder + "/Mesh_Cell3Ds.vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      std::vector<double> index, constrained, weight;
      index.resize(meshDAO.Cell2DTotalNumber());
      constrained.resize(meshDAO.Cell2DTotalNumber());
      weight.resize(meshDAO.Cell2DTotalNumber(), 0.0);
      for (unsigned int f = 0; f < meshDAO.Cell2DTotalNumber(); f++)
      {
        index[f] = f;
        constrained[f] = cell2DsConstrained[f] ? 1.0 : 0.0;
      }

      for (unsigned int e = 0; e < meshToNetwork.Network.EdgesWeight.size(); e++)
        weight[meshToNetwork.EdgesMeshCellIndex[e]] = meshToNetwork.Network.EdgesWeight[e];

      exporter.AddPolygons(meshDAO.Cell0DsCoordinates(),
                           meshDAO.Cell2DsVertices(),
                           {
                             {
                               "Index",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(index.size()),
                               index.data()
                             },
                             {
                               "Constrained",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(constrained.size()),
                               constrained.data()
                             },
                             {
                               "Weight",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(weight.size()),
                               weight.data()
                             }
                           });

      exporter.Export(exportFolder + "/Mesh_Cell2Ds.vtu");
    }

    for (unsigned int f = 0; f < cell2DsConstrained.size(); f++)
    {
      if (!cell2DsConstrained[f] ||
          meshDAO.Cell2DNumberNeighbourCell3D(f) < 2)
        continue;

      ASSERT_NE(fix_connectedComponents_partitions.at(meshDAO.Cell2DNeighbourCell3D(f,
                                                                                    0)),
                fix_connectedComponents_partitions.at(meshDAO.Cell2DNeighbourCell3D(f,
                                                                                    1)));
    }
  }
}

#endif // __TEST_METISUTILITIES_H
