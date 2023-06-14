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

#include "MeshMatrices_3D_22Cells_Mock.hpp"

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

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh2DToDualGraph(meshDAO);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 3;

    const std::vector<unsigned int> partition = metisUtilities.NetworkPartition(partitionOptions,
                                                                                meshToNetwork.MetisNetwork);


    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_DualGraph";
    Gedim::Output::CreateFolder(exportFolder);

    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell2DsCentroids[c];

      const Eigen::MatrixXi graphEdges = metisUtilities.GraphToConnectivityMatrix(meshToNetwork.MetisNetwork);

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

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh2DToDualGraph(meshDAO,
                                                                                                edgeConstrained);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partition = metisUtilities.NetworkPartition(partitionOptions,
                                                                                meshToNetwork.MetisNetwork);

    for (unsigned int e = 0; e < edgeConstrained.size(); e++)
    {
      if (!edgeConstrained[e] ||
          meshDAO.Cell1DNumberNeighbourCell2D(e) < 2)
        continue;

      ASSERT_NE(partition.at(meshDAO.Cell1DNeighbourCell2D(e,
                                                           0)),
                partition.at(meshDAO.Cell1DNeighbourCell2D(e,
                                                           1)));
    }

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_DualGraph_Constraints";
    Gedim::Output::CreateFolder(exportFolder);

    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell2DsCentroids[c];

      const Eigen::MatrixXi graphEdges = metisUtilities.GraphToConnectivityMatrix(meshToNetwork.MetisNetwork);

      std::vector<double> weights;
      weights.reserve(meshToNetwork.MetisNetwork.EdgeWeights.size());
      weights.assign(meshToNetwork.MetisNetwork.EdgeWeights.begin(),
                     meshToNetwork.MetisNetwork.EdgeWeights.end());

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

      exporter.Export(exportFolder + "/Cell2Ds.vtu");
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
        constrained[e] = edgeConstrained[e] ? 1.0 : 0.0;
      }

      for (unsigned int e = 0; e < meshToNetwork.MetisNetwork.EdgeWeights.size(); e++)
        weight[meshToNetwork.EdgesMeshCellIndex[e]] = meshToNetwork.MetisNetwork.EdgeWeights[e];

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

      exporter.Export(exportFolder + "/Cell1Ds.vtu");
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

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh2DToDualGraph(meshDAO,
                                                                                                edgeIsConstrained,
                                                                                                weights);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 50;

    const std::vector<unsigned int> partition = metisUtilities.NetworkPartition(partitionOptions,
                                                                                meshToNetwork.MetisNetwork);


    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell2DsCentroids[c];

      const Eigen::MatrixXi graphEdges = metisUtilities.GraphToConnectivityMatrix(meshToNetwork.MetisNetwork);

      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(partition.size());
      property.assign(partition.begin(), partition.end());

      std::vector<double> weights;
      weights.reserve(meshToNetwork.MetisNetwork.EdgeWeights.size());
      weights.assign(meshToNetwork.MetisNetwork.EdgeWeights.begin(),
                     meshToNetwork.MetisNetwork.EdgeWeights.end());

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

  TEST(TestMetisUtilities, TestNetworkPartition_Mesh3D_DualGraph)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

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

    const std::vector<unsigned int> partition = metisUtilities.NetworkPartition(partitionOptions,
                                                                                meshToNetwork.MetisNetwork);


    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh3D_DualGraph";
    Gedim::Output::CreateFolder(exportFolder);

    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell3DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell3DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell3DsCentroids[c];

      const Eigen::MatrixXi graphEdges = metisUtilities.GraphToConnectivityMatrix(meshToNetwork.MetisNetwork);

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

      exporter.Export(exportFolder + "/Cell3Ds.vtu");
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

      exporter.Export(exportFolder + "/Cell2Ds.vtu");
    }
  }

  TEST(TestMetisUtilities, TestNetworkPartition_Mesh3D_DualGraph_Constraints)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_3D_22Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    meshUtilities.ComputeCell2DCell3DNeighbours(meshDAO);

    const Gedim::MeshUtilities::MeshGeometricData3D geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                          meshDAO);

    std::vector<bool> facesConstrained = std::vector<bool>(meshDAO.Cell2DTotalNumber(),
                                                           false);
    facesConstrained[2] = true;
    facesConstrained[23] = true;
    facesConstrained[29] = true;

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh3DToDualGraph(meshDAO,
                                                                                                facesConstrained);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 3;

    const std::vector<unsigned int> partition = metisUtilities.NetworkPartition(partitionOptions,
                                                                                meshToNetwork.MetisNetwork);

    const std::vector<unsigned int> fixed_partition = metisUtilities.Mesh3DPartitionCheckConstraints(meshDAO,
                                                                                                     facesConstrained,
                                                                                                     partition);

    for (unsigned int e = 0; e < meshToNetwork.MetisNetwork.EdgeWeights.size(); e++)
    {
      const unsigned int cell2DIndex = meshToNetwork.EdgesMeshCellIndex[e];
      const unsigned int neigh1 = meshDAO.Cell2DNeighbourCell3D(cell2DIndex,
                                                                0);
      const unsigned int neigh2 = meshDAO.Cell2DNeighbourCell3D(cell2DIndex,
                                                                1);

      cerr<< "Face "<< cell2DIndex<< " ";
      cerr<< "constrained "<< facesConstrained[cell2DIndex]<< " ";
      cerr<< "weight "<< meshToNetwork.MetisNetwork.EdgeWeights[e]<< " - ";
      cerr<< "neigh1 "<< neigh1<< " ";
      cerr<< "P "<< fixed_partition.at(neigh1)<< " ";
      cerr<< "neigh2 "<< neigh2<< " ";
      cerr<< "P "<< fixed_partition.at(neigh2)<< endl;
    }

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh3D_DualGraph_Constraints";
    Gedim::Output::CreateFolder(exportFolder);

    {
      Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell3DTotalNumber());
      for (unsigned int c = 0; c < meshDAO.Cell3DTotalNumber(); c++)
        graphVertices.col(c)<< geometricData.Cell3DsCentroids[c];

      const Eigen::MatrixXi graphEdges = metisUtilities.GraphToConnectivityMatrix(meshToNetwork.MetisNetwork);

      Gedim::VTKUtilities exporter;

      std::vector<double> property;
      property.reserve(fixed_partition.size());
      property.assign(fixed_partition.begin(), fixed_partition.end());

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
      property.reserve(fixed_partition.size());
      property.assign(fixed_partition.begin(), fixed_partition.end());

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

      exporter.Export(exportFolder + "/Cell3Ds.vtu");
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
        constrained[f] = facesConstrained[f] ? 1.0 : 0.0;
      }

      for (unsigned int e = 0; e < meshToNetwork.MetisNetwork.EdgeWeights.size(); e++)
        weight[meshToNetwork.EdgesMeshCellIndex[e]] = meshToNetwork.MetisNetwork.EdgeWeights[e];

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

      exporter.Export(exportFolder + "/Cell2Ds.vtu");
    }

    for (unsigned int f = 0; f < facesConstrained.size(); f++)
    {
      if (!facesConstrained[f] ||
          meshDAO.Cell2DNumberNeighbourCell3D(f) < 2)
        continue;

      ASSERT_NE(fixed_partition.at(meshDAO.Cell2DNeighbourCell3D(f,
                                                                 0)),
                fixed_partition.at(meshDAO.Cell2DNeighbourCell3D(f,
                                                                 1)));
    }
  }
}

#endif // __TEST_METISUTILITIES_H
