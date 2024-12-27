#ifndef __TEST_METISUTILITIES_H
#define __TEST_METISUTILITIES_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "FileTextReader.hpp"
#include "Gedim_Macro.hpp"
#include "GeometryUtilities.hpp"
#include "GraphUtilities.hpp"
#include "IOUtilities.hpp"
#include "MeshDAOImporterFromCsv.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "MetisUtilities.hpp"
#include "OpenVolumeMeshInterface.hpp"
#include "VTKUtilities.hpp"

#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "MeshMatrices_2D_4Cells_Mock.hpp"

#include "MeshMatrices_3D_22Cells_Mock.hpp"
#include "MeshMatrices_3D_329Cells_Mock.hpp"
#include "OVM_Mesh_Mock.hpp"

#include "OFF_Mesh_Mock.hpp"
#include "ObjectFileFormatInterface.hpp"

#include "CommonUtilities.hpp"

namespace UnitTesting
{
void ExportMetis3DToVTU(const Gedim::IMeshDAO &meshDAO,
                        const Gedim::MeshUtilities::MeshGeometricData3D &geometricData,
                        const Gedim::MetisUtilities::MeshToNetwork &meshToNetwork,
                        const Gedim::GraphUtilities &graphUtilities,
                        const Gedim::MetisUtilities &metisUtilities,
                        const std::vector<unsigned int> &partitions,
                        const std::vector<unsigned int> &fix_constraints_partitions,
                        const std::vector<unsigned int> &fix_connectedComponents_partitions,
                        const string &exportFolder)
{
    {
        Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell3DTotalNumber());
        for (unsigned int c = 0; c < meshDAO.Cell3DTotalNumber(); c++)
            graphVertices.col(c) << geometricData.Cell3DsCentroids[c];

        const unsigned int graphNumEdges = meshToNetwork.Network.Adjacency.Cols.size();
        const std::vector<std::vector<unsigned int>> graphAdjacency =
            metisUtilities.MetisAdjacencyToGraphAdjacency(meshToNetwork.Network.Adjacency);
        const Eigen::MatrixXi graphEdges = graphUtilities.GraphAdjacencyToGraphConnectivity(graphNumEdges, graphAdjacency);

        {
            Gedim::VTKUtilities exporter;

            std::vector<double> weights;
            if (meshToNetwork.Network.NodesWeight.size() > 0)
            {
                weights.reserve(meshToNetwork.Network.NodesWeight.size());
                weights.assign(meshToNetwork.Network.NodesWeight.begin(), meshToNetwork.Network.NodesWeight.end());
            }
            else
                weights.resize(meshDAO.Cell3DTotalNumber(), 0.0);

            std::vector<double> partition;

            if (partitions.size() > 0)
            {
                partition.reserve(partitions.size());
                partition.assign(partitions.begin(), partitions.end());
            }
            else
                partition.resize(meshDAO.Cell3DTotalNumber(), 0.0);

            std::vector<double> fix_constraints_partition;

            if (fix_constraints_partitions.size() > 0)
            {
                fix_constraints_partition.reserve(fix_constraints_partitions.size());
                fix_constraints_partition.assign(fix_constraints_partitions.begin(), fix_constraints_partitions.end());
            }
            else
                fix_constraints_partition.resize(meshDAO.Cell3DTotalNumber(), 0.0);

            std::vector<double> fix_connectedComponents_partition;
            if (fix_connectedComponents_partitions.size() > 0)
            {
                fix_connectedComponents_partition.reserve(fix_connectedComponents_partitions.size());
                fix_connectedComponents_partition.assign(fix_connectedComponents_partitions.begin(),
                                                         fix_connectedComponents_partitions.end());
            }
            else
                fix_connectedComponents_partition.resize(meshDAO.Cell3DTotalNumber(), 0.0);

            exporter.AddPoints(
                graphVertices,
                {{"weights", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(weights.size()), weights.data()},
                 {"partition",
                  Gedim::VTPProperty::Formats::Cells,
                  static_cast<unsigned int>(partition.size()),
                  partition.data()},
                 {"fix_constraints_partition",
                  Gedim::VTPProperty::Formats::Cells,
                  static_cast<unsigned int>(fix_constraints_partition.size()),
                  fix_constraints_partition.data()},
                 {"fix_connectedComponents_partition",
                  Gedim::VTPProperty::Formats::Cells,
                  static_cast<unsigned int>(fix_connectedComponents_partition.size()),
                  fix_connectedComponents_partition.data()}});

            exporter.Export(exportFolder + "/Graph_Points.vtu");
        }

        {
            Gedim::VTKUtilities exporter;

            std::vector<double> property;
            property.reserve(meshToNetwork.Network.EdgesWeight.size());
            property.assign(meshToNetwork.Network.EdgesWeight.begin(), meshToNetwork.Network.EdgesWeight.end());

            exporter.AddSegments(
                graphVertices,
                graphEdges,
                {{"weights", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(property.size()), property.data()}});

            exporter.Export(exportFolder + "/Graph_Edges.vtu");
        }
    }

    {
        Gedim::VTKUtilities exporter;

        std::vector<double> property;
        property.reserve(fix_connectedComponents_partitions.size());
        property.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());

        exporter.AddPolyhedrons(
            meshDAO.Cell0DsCoordinates(),
            meshDAO.Cell3DsFacesVertices(),
            {{"Partition", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(property.size()), property.data()}});

        exporter.Export(exportFolder + "/Mesh_Cell3Ds.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        std::vector<double> index;
        index.resize(meshDAO.Cell2DTotalNumber());
        for (unsigned int f = 0; f < meshDAO.Cell2DTotalNumber(); f++)
            index[f] = f;

        exporter.AddPolygons(
            meshDAO.Cell0DsCoordinates(),
            meshDAO.Cell2DsVertices(),
            {{"Index", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(index.size()), index.data()}});

        exporter.Export(exportFolder + "/Mesh_Cell2Ds.vtu");
    }
}

TEST(TestMetisUtilities, TestNetworkPartition)
{
    Gedim::MetisUtilities metisUtilities;

    const unsigned int numberVertices = 6;
    Eigen::MatrixXi edges(2, 6);
    edges.col(0) << 0, 1;
    edges.col(1) << 1, 2;
    edges.col(2) << 1, 3;
    edges.col(3) << 1, 5;
    edges.col(4) << 2, 4;
    edges.col(5) << 2, 5;

    const Gedim::MetisUtilities::MetisNetwork network = metisUtilities.MeshToGraph(numberVertices, edges, true);

    ASSERT_EQ(std::vector<unsigned int>({0, 1, 5, 8, 9, 10, 12}), network.Adjacency.Rows);
    ASSERT_EQ(std::vector<unsigned int>({1, 0, 2, 3, 5, 1, 4, 5, 1, 2, 1, 2}), network.Adjacency.Cols);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions, network);

    const std::vector<unsigned int> fix_constraints_partitions = metisUtilities.PartitionCheckConstraints(network, partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions =
        metisUtilities.PartitionCheckConnectedComponents(network, fix_constraints_partitions);

#if ENABLE_METIS == 1
    ASSERT_EQ(std::vector<unsigned int>({1, 1, 0, 1, 0, 0}), partitions);
#else
    ASSERT_EQ(std::vector<unsigned int>({0, 0, 0, 0, 0, 0}), partitions);
#endif
}

TEST(TestMetisUtilities, TestNetworkPartition_Mesh2D_Graph)
{
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    const Eigen::MatrixXi edges = meshDAO.Cell1DsExtremes();

    const Gedim::MetisUtilities::MetisNetwork network = metisUtilities.MeshToGraph(meshDAO.Cell0DTotalNumber(), edges, true);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions, network);

    const std::vector<unsigned int> fix_constraints_partitions = metisUtilities.PartitionCheckConstraints(network, partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions =
        metisUtilities.PartitionCheckConnectedComponents(network, fix_constraints_partitions);

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_Graph";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Gedim::VTKUtilities exporter;

        std::vector<double> property;
        property.reserve(fix_connectedComponents_partitions.size());
        property.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());

        exporter.AddSegments(
            meshDAO.Cell0DsCoordinates(),
            edges,
            {{"partition", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(property.size()), property.data()}});

        exporter.Export(exportFolder + "/Mesh_Cell1Ds.vtu");
    }
}

TEST(TestMetisUtilities, TestNetworkPartition_Mesh2D_Graph_Weights)
{
    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    const Eigen::MatrixXi edges = meshDAO.Cell1DsExtremes();

    std::vector<unsigned int> verticesWeights(meshDAO.Cell0DTotalNumber(), 1.0);
    verticesWeights[3] = 100.0;
    verticesWeights[4] = 100.0;
    verticesWeights[5] = 100.0;

    std::vector<unsigned int> edgesWeights(meshDAO.Cell1DTotalNumber(), 1.0);
    edgesWeights[3] = 100.0;
    edgesWeights[7] = 100.0;
    edgesWeights[9] = 100.0;

    Gedim::MetisUtilities::MetisNetwork network =
        metisUtilities.MeshToGraph(meshDAO.Cell0DTotalNumber(), edges, true, verticesWeights, edgesWeights);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions, network);

    const std::vector<unsigned int> fix_constraints_partitions = metisUtilities.PartitionCheckConstraints(network, partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions =
        metisUtilities.PartitionCheckConnectedComponents(network, fix_constraints_partitions);

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_Graph_Weights";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::MeshUtilities meshUtilities;
    meshUtilities.ExportMeshToVTU(meshDAO, exportFolder, "Mesh");

    {
        Gedim::VTKUtilities exporter;

        std::vector<double> property;
        property.reserve(fix_connectedComponents_partitions.size());
        property.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());
        std::vector<double> weights;
        if (network.NodesWeight.size() > 0)
        {
            weights.reserve(network.NodesWeight.size());
            weights.assign(network.NodesWeight.begin(), network.NodesWeight.end());
        }
        else
            weights.resize(meshDAO.Cell0DTotalNumber(), 1.0);

        exporter.AddPoints(
            meshDAO.Cell0DsCoordinates(),
            {{"weights", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(weights.size()), weights.data()},
             {"partition", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(property.size()), property.data()}});

        exporter.Export(exportFolder + "/Network_Vertices.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        const unsigned int graphNumEdges = network.Adjacency.Cols.size();
        const std::vector<std::vector<unsigned int>> graphAdjacency =
            metisUtilities.MetisAdjacencyToGraphAdjacency(network.Adjacency);
        const Eigen::MatrixXi graphEdges = graphUtilities.GraphAdjacencyToGraphConnectivity(graphNumEdges, graphAdjacency);

        std::vector<double> weights;
        if (network.EdgesWeight.size() > 0)
        {
            weights.reserve(network.EdgesWeight.size());
            weights.assign(network.EdgesWeight.begin(), network.EdgesWeight.end());
        }
        else
            weights.resize(graphNumEdges, 1.0);

        exporter.AddSegments(
            meshDAO.Cell0DsCoordinates(),
            graphEdges,
            {{"weights", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(weights.size()), weights.data()}});

        exporter.Export(exportFolder + "/Network_Edges.vtu");
    }
}

TEST(TestMetisUtilities, TestNetworkPartition_Mesh2D_DualGraph)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    const Gedim::MeshUtilities::MeshGeometricData2D geometricData =
        meshUtilities.FillMesh2DGeometricData(geometryUtilities, meshDAO);

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh2DToDualGraph(meshDAO);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 3;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions, meshToNetwork.Network);

    const std::vector<unsigned int> fix_constraints_partitions =
        metisUtilities.PartitionCheckConstraints(meshToNetwork.Network, partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions =
        metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network, fix_constraints_partitions);

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_DualGraph";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
        for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
            graphVertices.col(c) << geometricData.Cell2DsCentroids[c];

        const unsigned int graphNumEdges = meshToNetwork.Network.Adjacency.Cols.size();
        const std::vector<std::vector<unsigned int>> graphAdjacency =
            metisUtilities.MetisAdjacencyToGraphAdjacency(meshToNetwork.Network.Adjacency);
        const Eigen::MatrixXi graphEdges = graphUtilities.GraphAdjacencyToGraphConnectivity(graphNumEdges, graphAdjacency);

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
            fix_connectedComponents_partition.assign(fix_connectedComponents_partitions.begin(),
                                                     fix_connectedComponents_partitions.end());

            exporter.AddPoints(graphVertices,
                               {{"partition",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(partition.size()),
                                 partition.data()},
                                {"fix_constraints_partition",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(fix_constraints_partition.size()),
                                 fix_constraints_partition.data()},
                                {"fix_connectedComponents_partition",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(fix_connectedComponents_partition.size()),
                                 fix_connectedComponents_partition.data()}});

            exporter.Export(exportFolder + "/Graph_Points.vtu");
        }

        {
            Gedim::VTKUtilities exporter;

            std::vector<double> property;
            property.reserve(meshToNetwork.Network.EdgesWeight.size());
            property.assign(meshToNetwork.Network.EdgesWeight.begin(), meshToNetwork.Network.EdgesWeight.end());

            exporter.AddSegments(
                graphVertices,
                graphEdges,
                {{"Weight", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(property.size()), property.data()}});

            exporter.Export(exportFolder + "/Graph_Edges.vtu");
        }
    }

    {
        Gedim::VTKUtilities exporter;

        std::vector<double> property;
        property.reserve(fix_connectedComponents_partitions.size());
        property.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());

        exporter.AddPolygons(
            meshDAO.Cell0DsCoordinates(),
            meshDAO.Cell2DsVertices(),
            {{"Partition", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(property.size()), property.data()}});

        exporter.Export(exportFolder + "/Mesh_Cell2Ds.vtu");
    }
}

TEST(TestMetisUtilities, TestNetworkPartition_Mesh2D_DualGraph_Constraints)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    const Gedim::MeshUtilities::MeshGeometricData2D geometricData =
        meshUtilities.FillMesh2DGeometricData(geometryUtilities, meshDAO);

    std::vector<bool> cell1DsConstrained = std::vector<bool>(meshDAO.Cell1DTotalNumber(), false);
    cell1DsConstrained[2] = true;
    cell1DsConstrained[5] = true;

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh2DToDualGraph(meshDAO, {}, cell1DsConstrained);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 2;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions, meshToNetwork.Network);

    const std::vector<unsigned int> fix_constraints_partitions =
        metisUtilities.PartitionCheckConstraints(meshToNetwork.Network, partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions =
        metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network, fix_constraints_partitions);

    for (unsigned int e = 0; e < cell1DsConstrained.size(); e++)
    {
        if (!cell1DsConstrained[e] || meshDAO.Cell1DNumberNeighbourCell2D(e) < 2 ||
            !meshDAO.Cell1DHasNeighbourCell2D(e, 0) || !meshDAO.Cell1DHasNeighbourCell2D(e, 1))
            continue;

        ASSERT_NE(fix_connectedComponents_partitions.at(meshDAO.Cell1DNeighbourCell2D(e, 0)),
                  fix_connectedComponents_partitions.at(meshDAO.Cell1DNeighbourCell2D(e, 1)));
    }

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_DualGraph_Constraints";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
        for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
            graphVertices.col(c) << geometricData.Cell2DsCentroids[c];

        const unsigned int graphNumEdges = meshToNetwork.Network.Adjacency.Cols.size();
        const std::vector<std::vector<unsigned int>> graphAdjacency =
            metisUtilities.MetisAdjacencyToGraphAdjacency(meshToNetwork.Network.Adjacency);
        const Eigen::MatrixXi graphEdges = graphUtilities.GraphAdjacencyToGraphConnectivity(graphNumEdges, graphAdjacency);
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
            fix_connectedComponents_partition.assign(fix_connectedComponents_partitions.begin(),
                                                     fix_connectedComponents_partitions.end());

            exporter.AddPoints(graphVertices,
                               {{"partition",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(partition.size()),
                                 partition.data()},
                                {"fix_constraints_partition",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(fix_constraints_partition.size()),
                                 fix_constraints_partition.data()},
                                {"fix_connectedComponents_partition",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(fix_connectedComponents_partition.size()),
                                 fix_connectedComponents_partition.data()}});

            exporter.Export(exportFolder + "/Graph_Points.vtu");
        }

        {
            Gedim::VTKUtilities exporter;

            std::vector<double> property;
            property.reserve(meshToNetwork.Network.EdgesWeight.size());
            property.assign(meshToNetwork.Network.EdgesWeight.begin(), meshToNetwork.Network.EdgesWeight.end());

            exporter.AddSegments(
                graphVertices,
                graphEdges,
                {{"Weight", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(property.size()), property.data()}});

            exporter.Export(exportFolder + "/Graph_Edges.vtu");
        }
    }

    {
        Gedim::VTKUtilities exporter;

        std::vector<double> property;
        property.reserve(fix_connectedComponents_partitions.size());
        property.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());

        exporter.AddPolygons(
            meshDAO.Cell0DsCoordinates(),
            meshDAO.Cell2DsVertices(),
            {{"Partition", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(property.size()), property.data()}});

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

        exporter.AddSegments(
            meshDAO.Cell0DsCoordinates(),
            meshDAO.Cell1DsExtremes(),
            {{"Index", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(index.size()), index.data()},
             {"Constrained",
              Gedim::VTPProperty::Formats::Cells,
              static_cast<unsigned int>(constrained.size()),
              constrained.data()},
             {"Weight", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(weight.size()), weight.data()}});

        exporter.Export(exportFolder + "/Mesh_Cell1Ds.vtu");
    }
}

TEST(TestMetisUtilities, TestNetworkPartition_ImportMesh2D_DualGraph)
{
    GTEST_SKIP();

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_ImportMesh2D_DualGraph";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh);

    Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
    meshImporterConfiguration.Folder = "/home/geoscore/Dropbox/Polito/Articles/GE_MESH/Poj_MeshQuality_ToGe/"
                                       "NumericalTests/METIS/SingleDomain/M1";
    Gedim::MeshFromCsvUtilities meshFromCsvUtilities;
    Gedim::MeshDAOImporterFromCsv importer(meshFromCsvUtilities);
    importer.Import(meshImporterConfiguration, meshDAO);

    meshUtilities.ExportMeshToVTU(meshDAO, exportFolder, "ImportedMesh");

    const Gedim::MeshUtilities::MeshGeometricData2D geometricData =
        meshUtilities.FillMesh2DGeometricData(geometryUtilities, meshDAO);

    Eigen::SparseMatrix<unsigned int> weights(meshDAO.Cell2DTotalNumber(), meshDAO.Cell2DTotalNumber());

    {
        list<Eigen::Triplet<unsigned int>> triplets;

        Gedim::FileReader fileReader("/home/geoscore/Dropbox/Polito/Articles/GE_MESH/Poj_MeshQuality_ToGe/"
                                     "NumericalTests/METIS/SingleDomain/M1/Weights.txt");
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
                    triplets.push_back(
                        Eigen::Triplet<unsigned int>(l, c, round(meshDAO.Cell2DTotalNumber() * meshDAO.Cell2DTotalNumber() / weight)));
                }
            }
        }

        weights.setFromTriplets(triplets.begin(), triplets.end());
        weights.makeCompressed();
    }

    std::vector<bool> cell1DsConstrained = std::vector<bool>(meshDAO.Cell1DTotalNumber(), false);
    {
        const unsigned int constrainIndex = meshDAO.Cell1DDoublePropertyIndex("marked");
        for (unsigned int e = 0; e < meshDAO.Cell1DTotalNumber(); e++)
            cell1DsConstrained[e] = meshDAO.Cell1DDoublePropertyValue(e, constrainIndex, 0) == 1.0;
    }

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork =
        metisUtilities.Mesh2DToDualGraph(meshDAO, {}, cell1DsConstrained, weights);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 50;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions, meshToNetwork.Network);

    const std::vector<unsigned int> fix_constraints_partitions =
        metisUtilities.PartitionCheckConstraints(meshToNetwork.Network, partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions =
        metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network, fix_constraints_partitions);

    for (unsigned int e = 0; e < cell1DsConstrained.size(); e++)
    {
        if (!cell1DsConstrained[e] || meshDAO.Cell1DNumberNeighbourCell2D(e) < 2 ||
            !meshDAO.Cell1DHasNeighbourCell2D(e, 0) || !meshDAO.Cell1DHasNeighbourCell2D(e, 1))
            continue;

        ASSERT_NE(fix_connectedComponents_partitions.at(meshDAO.Cell1DNeighbourCell2D(e, 0)),
                  fix_connectedComponents_partitions.at(meshDAO.Cell1DNeighbourCell2D(e, 1)));
    }
}

TEST(TestMetisUtilities, TestNetworkPartition_Mesh2D_DualGraph_OFFFile)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;
    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh);

    Gedim::ObjectFileFormatInterface offInterface;

    const std::vector<string> lines = GedimUnitTesting::OFF_Mesh_Mock::FileLines();
    const Gedim::ObjectFileFormatInterface::OFFMesh off_mesh = offInterface.StringsToOFFMesh(lines);
    offInterface.OFFMeshToMeshDAO(off_mesh, meshUtilities, meshDAO);

    meshUtilities.ComputeCell1DCell2DNeighbours(meshDAO);

    const Gedim::MeshUtilities::MeshGeometricData2D geometricData =
        meshUtilities.FillMesh2DGeometricData(geometryUtilities, meshDAO);

    std::vector<bool> cell1DsConstrained = {};
    if (false)
    {
        std::vector<bool> cell1DsConstrained = std::vector<bool>(meshDAO.Cell1DTotalNumber(), false);
        cell1DsConstrained[2] = true;
        cell1DsConstrained[5] = true;
    }

    Eigen::SparseMatrix<unsigned int> cell1DsWeight;

    if (false)
    {
        cell1DsWeight.resize(meshDAO.Cell2DTotalNumber(), meshDAO.Cell2DTotalNumber());

        constexpr double toCutEdgeJHorizWeigth = 1.0;
        constexpr double toCutEdgeVertWeigth = 5330.0;
        constexpr double toMantainEdgeWeigth = 26244.0;
        list<Eigen::Triplet<unsigned int>> triplets;

        for (unsigned int e = 0; e < meshDAO.Cell1DTotalNumber(); e++)
        {
            if (meshDAO.Cell1DNumberNeighbourCell2D(e) < 2 || !meshDAO.Cell1DHasNeighbourCell2D(e, 0) ||
                !meshDAO.Cell1DHasNeighbourCell2D(e, 1))
                continue;

            if (geometryUtilities.AreValuesEqual(meshDAO.Cell1DOriginCoordinates(e).x(),
                                                 meshDAO.Cell1DEndCoordinates(e).x(),
                                                 geometryUtilities.Tolerance1D()))
            {
                triplets.push_back(Eigen::Triplet<unsigned int>(meshDAO.Cell1DNeighbourCell2D(e, 0),
                                                                meshDAO.Cell1DNeighbourCell2D(e, 1),
                                                                std::round(toCutEdgeVertWeigth)));
                triplets.push_back(Eigen::Triplet<unsigned int>(meshDAO.Cell1DNeighbourCell2D(e, 1),
                                                                meshDAO.Cell1DNeighbourCell2D(e, 0),
                                                                std::round(toCutEdgeVertWeigth)));
            }
            else if (geometryUtilities.AreValuesEqual(meshDAO.Cell1DOriginCoordinates(e).y(),
                                                      meshDAO.Cell1DEndCoordinates(e).y(),
                                                      geometryUtilities.Tolerance1D()))
            {
                triplets.push_back(Eigen::Triplet<unsigned int>(meshDAO.Cell1DNeighbourCell2D(e, 0),
                                                                meshDAO.Cell1DNeighbourCell2D(e, 1),
                                                                std::round(toCutEdgeJHorizWeigth)));
                triplets.push_back(Eigen::Triplet<unsigned int>(meshDAO.Cell1DNeighbourCell2D(e, 1),
                                                                meshDAO.Cell1DNeighbourCell2D(e, 0),
                                                                std::round(toCutEdgeJHorizWeigth)));
            }
            else
            {
                triplets.push_back(Eigen::Triplet<unsigned int>(meshDAO.Cell1DNeighbourCell2D(e, 0),
                                                                meshDAO.Cell1DNeighbourCell2D(e, 1),
                                                                std::round(toMantainEdgeWeigth)));
                triplets.push_back(Eigen::Triplet<unsigned int>(meshDAO.Cell1DNeighbourCell2D(e, 1),
                                                                meshDAO.Cell1DNeighbourCell2D(e, 0),
                                                                std::round(toMantainEdgeWeigth)));
            }
        }

        cell1DsWeight.setFromTriplets(triplets.begin(), triplets.end());
        cell1DsWeight.makeCompressed();
    }
    else if (false)
    {
        cell1DsWeight.resize(meshDAO.Cell2DTotalNumber(), meshDAO.Cell2DTotalNumber());

        list<Eigen::Triplet<unsigned int>> triplets;

        Gedim::FileReader fileReader("/home/geoscore/Downloads/weights.txt");
        fileReader.Open();

        std::vector<string> lines;
        fileReader.GetAllLines(lines);
        fileReader.Close();

        for (unsigned int l = 0; l < lines.size(); l++)
        {
            istringstream converter(lines[l]);

            unsigned int i, j, weight;
            converter >> i >> j >> weight;

            triplets.push_back(Eigen::Triplet<unsigned int>(i, j, weight));
        }

        cell1DsWeight.setFromTriplets(triplets.begin(), triplets.end());
        cell1DsWeight.makeCompressed();
    }

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork =
        metisUtilities.Mesh2DToDualGraph(meshDAO, {}, cell1DsConstrained, cell1DsWeight);

    std::vector<unsigned int> partitions, fix_constraints_partitions, fix_connectedComponents_partitions;

    std::vector<double> toTests = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.5};
    for (double p : toTests)
    {

        Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
        partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
        partitionOptions.MasterWeight = 100;
        partitionOptions.NumberOfParts = p * meshDAO.Cell2DTotalNumber();
        partitionOptions.ContigousPartitions = true;
        partitionOptions.CompressGraph = true;
        partitionOptions.MinimizeConnectivity = true;
        partitionOptions.DebugLevel = Gedim::MetisUtilities::NetworkPartitionOptions::DebugLevels::None;
        partitionOptions.CoarseningSchema = Gedim::MetisUtilities::NetworkPartitionOptions::CoarseningSchemes::Default;
        partitionOptions.RefinementSchema = Gedim::MetisUtilities::NetworkPartitionOptions::RefinementSchemes::Default;
        partitionOptions.NumberRefinementIterations = 10;
        partitionOptions.InitialPartitioningSchema =
            Gedim::MetisUtilities::NetworkPartitionOptions::InitialPartitioningSchemes::Default;
        partitionOptions.RandomSeed = -1;

        partitions = metisUtilities.NetworkPartition(partitionOptions, meshToNetwork.Network);

        fix_constraints_partitions = metisUtilities.PartitionCheckConstraints(meshToNetwork.Network, partitions);

        fix_connectedComponents_partitions =
            metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network, fix_constraints_partitions);

        if (partitionOptions.DebugLevel != Gedim::MetisUtilities::NetworkPartitionOptions::DebugLevels::None)
        {
            using namespace Gedim;

            std::cout << "Perc " << p << " ";
            std::cout << "partition: " << *std::max_element(partitions.begin(), partitions.end()) + 1 << " / "
                      << partitionOptions.NumberOfParts << " ";
            std::cout
                << "fix_const: " << *std::max_element(fix_constraints_partitions.begin(), fix_constraints_partitions.end()) + 1
                << " / " << partitionOptions.NumberOfParts << " ";
            std::cout
                << "fix_conne: "
                << *std::max_element(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end()) + 1
                << " / " << partitionOptions.NumberOfParts << std::endl;
        }
    }

    for (unsigned int e = 0; e < cell1DsConstrained.size(); e++)
    {
        if (!cell1DsConstrained[e] || meshDAO.Cell1DNumberNeighbourCell2D(e) < 2 ||
            !meshDAO.Cell1DHasNeighbourCell2D(e, 0) || !meshDAO.Cell1DHasNeighbourCell2D(e, 1))
            continue;

        ASSERT_NE(fix_connectedComponents_partitions.at(meshDAO.Cell1DNeighbourCell2D(e, 0)),
                  fix_connectedComponents_partitions.at(meshDAO.Cell1DNeighbourCell2D(e, 1)));
    }

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh2D_DualGraph_OFFFile";
    Gedim::Output::CreateFolder(exportFolder);

    {
        Eigen::MatrixXd graphVertices = Eigen::MatrixXd::Zero(3, meshDAO.Cell2DTotalNumber());
        for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
            graphVertices.col(c) << geometricData.Cell2DsCentroids[c];

        const unsigned int graphNumEdges = meshToNetwork.Network.Adjacency.Cols.size();
        const std::vector<std::vector<unsigned int>> graphAdjacency =
            metisUtilities.MetisAdjacencyToGraphAdjacency(meshToNetwork.Network.Adjacency);
        const Eigen::MatrixXi graphEdges = graphUtilities.GraphAdjacencyToGraphConnectivity(graphNumEdges, graphAdjacency);
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
            fix_connectedComponents_partition.assign(fix_connectedComponents_partitions.begin(),
                                                     fix_connectedComponents_partitions.end());

            exporter.AddPoints(graphVertices,
                               {{"partition",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(partition.size()),
                                 partition.data()},
                                {"fix_constraints_partition",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(fix_constraints_partition.size()),
                                 fix_constraints_partition.data()},
                                {"fix_connectedComponents_partition",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(fix_connectedComponents_partition.size()),
                                 fix_connectedComponents_partition.data()}});

            exporter.Export(exportFolder + "/Graph_Points.vtu");
        }

        {
            Gedim::VTKUtilities exporter;

            std::vector<double> property;
            property.reserve(meshToNetwork.Network.EdgesWeight.size());
            property.assign(meshToNetwork.Network.EdgesWeight.begin(), meshToNetwork.Network.EdgesWeight.end());

            exporter.AddSegments(
                graphVertices,
                graphEdges,
                {{"Weight", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(property.size()), property.data()}});

            exporter.Export(exportFolder + "/Graph_Edges.vtu");
        }
    }

    {
        Gedim::VTKUtilities exporter;

        std::vector<double> property;
        property.reserve(fix_connectedComponents_partitions.size());
        property.assign(fix_connectedComponents_partitions.begin(), fix_connectedComponents_partitions.end());

        exporter.AddPolygons(
            meshDAO.Cell0DsCoordinates(),
            meshDAO.Cell2DsVertices(),
            {{"Partition", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(property.size()), property.data()}});

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
            constrained[e] = cell1DsConstrained.size() > 0 && cell1DsConstrained[e] ? 1.0 : 0.0;
        }

        for (unsigned int e = 0; e < meshToNetwork.Network.EdgesWeight.size(); e++)
            weight[meshToNetwork.EdgesMeshCellIndex[e]] = meshToNetwork.Network.EdgesWeight[e];

        exporter.AddSegments(
            meshDAO.Cell0DsCoordinates(),
            meshDAO.Cell1DsExtremes(),
            {{"Index", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(index.size()), index.data()},
             {"Constrained",
              Gedim::VTPProperty::Formats::Cells,
              static_cast<unsigned int>(constrained.size()),
              constrained.data()},
             {"Weight", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(weight.size()), weight.data()}});

        exporter.Export(exportFolder + "/Mesh_Cell1Ds.vtu");
    }
}

TEST(TestMetisUtilities, TestNetworkPartition_Mesh3D_DualGraph)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_3D_22Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    meshUtilities.ComputeCell2DCell3DNeighbours(meshDAO);

    const Gedim::MeshUtilities::MeshGeometricData3D geometricData =
        meshUtilities.FillMesh3DGeometricData(geometryUtilities, meshDAO);

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh3DToDualGraph(meshDAO);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 3;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions, meshToNetwork.Network);

    const std::vector<unsigned int> fix_constraints_partitions =
        metisUtilities.PartitionCheckConstraints(meshToNetwork.Network, partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions =
        metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network, fix_constraints_partitions);

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh3D_DualGraph";
    Gedim::Output::CreateFolder(exportFolder);

    ExportMetis3DToVTU(meshDAO, geometricData, meshToNetwork, graphUtilities, metisUtilities, partitions, fix_constraints_partitions, fix_connectedComponents_partitions, exportFolder);
}

TEST(TestMetisUtilities, TestNetworkPartition_Mesh3D_DualGraph_Constraints)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    GedimUnitTesting::MeshMatrices_3D_22Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    meshUtilities.ComputeCell2DCell3DNeighbours(meshDAO);

    const Gedim::MeshUtilities::MeshGeometricData3D geometricData =
        meshUtilities.FillMesh3DGeometricData(geometryUtilities, meshDAO);

    std::vector<bool> cell2DsConstrained = std::vector<bool>(meshDAO.Cell2DTotalNumber(), false);
    cell2DsConstrained[2] = true;
    cell2DsConstrained[23] = true;
    cell2DsConstrained[29] = true;

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh3DToDualGraph(meshDAO, {}, cell2DsConstrained);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 3;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions, meshToNetwork.Network);

    const std::vector<unsigned int> fix_constraints_partitions =
        metisUtilities.PartitionCheckConstraints(meshToNetwork.Network, partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions =
        metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network, fix_constraints_partitions);

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh3D_DualGraph_Constraints";
    Gedim::Output::CreateFolder(exportFolder);

    ExportMetis3DToVTU(meshDAO, geometricData, meshToNetwork, graphUtilities, metisUtilities, partitions, fix_constraints_partitions, fix_connectedComponents_partitions, exportFolder);

    for (unsigned int f = 0; f < cell2DsConstrained.size(); f++)
    {
        if (!cell2DsConstrained[f] || meshDAO.Cell2DNumberNeighbourCell3D(f) < 2)
            continue;

        ASSERT_NE(fix_connectedComponents_partitions.at(meshDAO.Cell2DNeighbourCell3D(f, 0)),
                  fix_connectedComponents_partitions.at(meshDAO.Cell2DNeighbourCell3D(f, 1)));
    }
}

TEST(TestMetisUtilities, TestNetworkPartition_Mesh3D_DualGraph_OVM)
{
    GTEST_SKIP();

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    std::vector<std::vector<bool>> meshCell3DsFacesOrientation;
    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh);

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_Mesh3D_DualGraph_OVM";
    Gedim::Output::CreateFolder(exportFolder);

    meshUtilities.ImportOpenVolumeMesh("/home/geoscore/Downloads/mesh.ovm", meshDAO, meshCell3DsFacesOrientation);

    meshUtilities.ComputeCell2DCell3DNeighbours(meshDAO);
    const Gedim::MeshUtilities::MeshGeometricData3D geometricData =
        meshUtilities.FillMesh3DGeometricData(geometryUtilities, meshDAO);

    constexpr unsigned int nodeWeightsActive = true;
    constexpr unsigned int edgeWeightsActive = false;

    std::vector<unsigned int> nodesWeights;
    Eigen::SparseMatrix<unsigned int> edgesWeights;

    if (nodeWeightsActive)
    {
        nodesWeights.resize(meshDAO.Cell3DTotalNumber(), 1);
        std::vector<double> nodesW(meshDAO.Cell3DTotalNumber());
        for (unsigned int c = 0; c < meshDAO.Cell3DTotalNumber(); c++)
            nodesW[c] = c + 1;

        double min = *std::min_element(nodesW.begin(), nodesW.end());
        double max = *std::max_element(nodesW.begin(), nodesW.end());

        ofstream file(exportFolder + "/weights.csv");
        file << "cell"
             << ",";
        file << "nodesW"
             << ",";
        file << "nodesWeights" << endl;

        for (unsigned int c = 0; c < meshDAO.Cell3DTotalNumber(); c++)
        {
            nodesWeights[c] = (max == min) ? std::round(nodesW[c] / 1000) + 1 :
                                           // std::round(nodesW[c] *
                                           //                        meshDAO.Cell3DTotalNumber() / 1000) + 1 :
                                  std::round((nodesW.at(c) - min) / (max - min) * meshDAO.Cell3DTotalNumber() / 1000) + 1;

            file << c << ",";
            file << nodesW.at(c) << ",";
            file << nodesWeights.at(c) << endl;
        }
        file.close();

        std::cout << "nodesW_min " << *std::min_element(nodesW.begin(), nodesW.end()) << std::endl;
        std::cout << "nodesW_min " << *std::max_element(nodesW.begin(), nodesW.end()) << std::endl;
        std::cout << "nodesWeights_min " << *std::min_element(nodesWeights.begin(), nodesWeights.end()) << std::endl;
        std::cout << "nodesWeights_max " << *std::max_element(nodesWeights.begin(), nodesWeights.end()) << std::endl;
    }

    if (edgeWeightsActive)
    {
        edgesWeights.resize(meshDAO.Cell2DTotalNumber(), meshDAO.Cell2DTotalNumber());

        {
            list<Eigen::Triplet<unsigned int>> triplets;

            Gedim::FileReader fileReader("/home/geoscore/Dropbox/Polito/Articles/GE_MESH_3D/NumericalTests/"
                                         "Tet_poisson_5/weights.txt");
            fileReader.Open();

            std::vector<string> lines;
            fileReader.GetAllLines(lines);
            fileReader.Close();

            for (unsigned int l = 0; l < lines.size(); l++)
            {
                istringstream converter(lines[l]);

                unsigned int i, j, weight;
                converter >> i >> j >> weight;

                triplets.push_back(Eigen::Triplet<unsigned int>(i, j, weight));
            }

            edgesWeights.setFromTriplets(triplets.begin(), triplets.end());
            edgesWeights.makeCompressed();
        }
    }

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork =
        metisUtilities.Mesh3DToDualGraph(meshDAO, nodesWeights, {}, edgesWeights);

    ExportMetis3DToVTU(meshDAO, geometricData, meshToNetwork, graphUtilities, metisUtilities, {}, {}, {}, exportFolder);

    unsigned int startTest = 39, numTests = 1;
    std::vector<unsigned int> numPartions(numTests), askPartitions(numTests);
    for (unsigned int p = 1; p <= numTests; p++)
    {
        Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
        partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
        partitionOptions.MasterWeight = 100;
        partitionOptions.ContigousPartitions = true;
        partitionOptions.CompressGraph = true;
        partitionOptions.MinimizeConnectivity = true;
        partitionOptions.NumberOfParts = (meshDAO.Cell3DTotalNumber() * (p + startTest)) / 100;
        partitionOptions.RandomSeed = 10;

        askPartitions[p - 1] = partitionOptions.NumberOfParts;

        std::cout << "partitionOptions.NumberOfParts " << partitionOptions.NumberOfParts << std::endl;

        const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions, meshToNetwork.Network);

        const std::vector<unsigned int> fix_constraints_partitions =
            metisUtilities.PartitionCheckConstraints(meshToNetwork.Network, partitions);

        const std::vector<unsigned int> fix_connectedComponents_partitions =
            metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network, fix_constraints_partitions);

        std::cout << "found "
                  << (*std::max_element(begin(fix_connectedComponents_partitions), end(fix_connectedComponents_partitions)) + 1)
                  << std::endl;

        numPartions[p - 1] =
            (*std::max_element(begin(fix_connectedComponents_partitions), end(fix_connectedComponents_partitions)) + 1);

        ExportMetis3DToVTU(meshDAO, geometricData, meshToNetwork, graphUtilities, metisUtilities, partitions, fix_constraints_partitions, fix_connectedComponents_partitions, exportFolder);
    }

    ofstream file(exportFolder + "/data.csv");
    file << "percentage"
         << ",";
    file << "askPartions"
         << ",";
    file << "numPartions" << endl;

    for (unsigned int p = 1; p <= numTests; p++)
    {
        file << p << ",";
        file << askPartitions[p - 1] << ",";
        file << numPartions[p - 1] << endl;
    }

    file.close();
}

TEST(TestMetisUtilities, TestNetworkPartition_ExportMesh3D)
{
    GTEST_SKIP();

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::GraphUtilities graphUtilities;
    Gedim::MetisUtilities metisUtilities;

    std::vector<std::vector<bool>> meshCell3DsFacesOrientation;
    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh);

    Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
    meshImporterConfiguration.Folder = "/home/geoscore/Dropbox/Polito/Articles/GE_MESH_3D/NumericalTests/conformed";
    Gedim::MeshFromCsvUtilities meshFromCsvUtilities;
    Gedim::MeshDAOImporterFromCsv importer(meshFromCsvUtilities);
    importer.Import(meshImporterConfiguration, meshDAO);

    std::string exportFolder = "./Export/TestMetisUtilities/TestNetworkPartition_ExportMesh3D";
    Gedim::Output::CreateFolder(exportFolder);

    meshUtilities.ExportMeshToVTU(meshDAO, exportFolder, "ImportedMesh");

    meshUtilities.ComputeCell2DCell3DNeighbours(meshDAO);
    const Gedim::MeshUtilities::MeshGeometricData3D geometricData =
        meshUtilities.FillMesh3DGeometricData(geometryUtilities, meshDAO);

    meshUtilities.ExportMeshToOpenVolume(meshDAO, geometricData.Cell3DsFacesNormalDirections, exportFolder + "/mesh3D.ovm");

    const Gedim::MetisUtilities::MeshToNetwork meshToNetwork = metisUtilities.Mesh3DToDualGraph(meshDAO);

    Gedim::MetisUtilities::NetworkPartitionOptions partitionOptions;
    partitionOptions.PartitionType = Gedim::MetisUtilities::NetworkPartitionOptions::PartitionTypes::CutBalancing;
    partitionOptions.MasterWeight = 100;
    partitionOptions.NumberOfParts = 3;

    const std::vector<unsigned int> partitions = metisUtilities.NetworkPartition(partitionOptions, meshToNetwork.Network);

    const std::vector<unsigned int> fix_constraints_partitions =
        metisUtilities.PartitionCheckConstraints(meshToNetwork.Network, partitions);

    const std::vector<unsigned int> fix_connectedComponents_partitions =
        metisUtilities.PartitionCheckConnectedComponents(meshToNetwork.Network, fix_constraints_partitions);

    ExportMetis3DToVTU(meshDAO, geometricData, meshToNetwork, graphUtilities, metisUtilities, partitions, fix_constraints_partitions, fix_connectedComponents_partitions, exportFolder);
}
} // namespace UnitTesting

#endif // __TEST_METISUTILITIES_H
