#ifndef __TEST_GRAPHUTILITIES_H
#define __TEST_GRAPHUTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GraphUtilities.hpp"

namespace UnitTesting
{
  // ***************************************************************************
  TEST(TestGraphUtilities, TestGraphConnectivityAndGraphAdjacency)
  {
    Gedim::GraphUtilities graphUtilities;

    // Create a graph given in the above diagram
    const unsigned int graph_numVertices = 5;
    const unsigned int graph_numEdges = 5;
    Eigen::MatrixXi graph_connectivity_expected(2, graph_numEdges);

    graph_connectivity_expected.col(0)<< 0, 2;
    graph_connectivity_expected.col(1)<< 0, 3;
    graph_connectivity_expected.col(2)<< 1, 0;
    graph_connectivity_expected.col(3)<< 2, 1;
    graph_connectivity_expected.col(4)<< 3, 4;

    const std::vector<std::vector<unsigned int>> graph_adjacency = graphUtilities.GraphConnectivityToGraphAdjacency(graph_numVertices,
                                                                                                                    graph_connectivity_expected);

    ASSERT_EQ(5, graph_adjacency.size());
    ASSERT_EQ(std::vector<unsigned int>({ 2, 3 }), graph_adjacency[0]);
    ASSERT_EQ(std::vector<unsigned int>({ 0 }), graph_adjacency[1]);
    ASSERT_EQ(std::vector<unsigned int>({ 1 }), graph_adjacency[2]);
    ASSERT_EQ(std::vector<unsigned int>({ 4 }), graph_adjacency[3]);
    ASSERT_EQ(std::vector<unsigned int>({ }), graph_adjacency[4]);

    const Eigen::MatrixXi graph_connectivity = graphUtilities.GraphAdjacencyToGraphConnectivity(graph_numEdges,
                                                                                                graph_adjacency);

    ASSERT_EQ(graph_connectivity_expected, graph_connectivity);
  }
  // ***************************************************************************
  TEST(TestGraphUtilities, TestExtractSubGraph)
  {
    Gedim::GraphUtilities graphUtilities;

    // Create a graph given in the above diagram
    const unsigned int graph_numVertices = 5;
    const unsigned int graph_numEdges = 5;
    Eigen::MatrixXi graph_connectivity_expected(2, graph_numEdges);

    graph_connectivity_expected.col(0)<< 0, 2;
    graph_connectivity_expected.col(1)<< 0, 3;
    graph_connectivity_expected.col(2)<< 1, 0;
    graph_connectivity_expected.col(3)<< 2, 1;
    graph_connectivity_expected.col(4)<< 3, 4;

    const std::vector<std::vector<unsigned int>> graph_adjacency = graphUtilities.GraphConnectivityToGraphAdjacency(graph_numVertices,
                                                                                                                    graph_connectivity_expected);
    const std::unordered_map<unsigned int, unsigned int> graph_filter = { {0, 0}, {2, 1}, {3, 2}, {4, 3} };

    const std::vector<std::vector<unsigned int>> subGraph_adjacency = graphUtilities.ExtractSubGraph(graph_adjacency,
                                                                                                     graph_filter);

    ASSERT_EQ(4, subGraph_adjacency.size());
    ASSERT_EQ(std::vector<unsigned int>({ 1, 2 }), subGraph_adjacency[0]);
    ASSERT_EQ(std::vector<unsigned int>({ }), subGraph_adjacency[1]);
    ASSERT_EQ(std::vector<unsigned int>({ 3 }), subGraph_adjacency[2]);
    ASSERT_EQ(std::vector<unsigned int>({ }), subGraph_adjacency[3]);
  }
  // ***************************************************************************
  TEST(TestGraphUtilities, TestStronglyConnectedComponents)
  {
    Gedim::GraphUtilities graphUtilities;

    // Create a graph given in the above diagram
    const unsigned int graph_numVertices = 7;
    const unsigned int graph_numEdges = 7;
    Eigen::MatrixXi graph_connectivity(2, graph_numEdges);

    graph_connectivity.col(0)<< 0, 2;
    graph_connectivity.col(1)<< 0, 3;
    graph_connectivity.col(2)<< 1, 0;
    graph_connectivity.col(3)<< 2, 1;
    graph_connectivity.col(4)<< 3, 4;
    graph_connectivity.col(5)<< 6, 5;
    graph_connectivity.col(6)<< 5, 6;

    const std::vector<std::vector<unsigned int>> graph_adjacency = graphUtilities.GraphConnectivityToGraphAdjacency(graph_numVertices,
                                                                                                                    graph_connectivity);

    const std::vector<std::vector<unsigned int>> stronglyConnectedComponents = graphUtilities.ComputeStronglyConnectedComponents(graph_adjacency);

    ASSERT_EQ(4, stronglyConnectedComponents.size());
    ASSERT_EQ(std::vector<unsigned int>({ 5, 6 }), stronglyConnectedComponents[0]);
    ASSERT_EQ(std::vector<unsigned int>({ 0, 1, 2 }), stronglyConnectedComponents[1]);
    ASSERT_EQ(std::vector<unsigned int>({ 3 }), stronglyConnectedComponents[2]);
    ASSERT_EQ(std::vector<unsigned int>({ 4 }), stronglyConnectedComponents[3]);
  }
  // ***************************************************************************
}

#endif // __TEST_METISUTILITIES_H
