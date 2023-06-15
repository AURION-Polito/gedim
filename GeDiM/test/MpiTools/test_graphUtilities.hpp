#ifndef __TEST_GRAPHUTILITIES_H
#define __TEST_GRAPHUTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GraphUtilities.hpp"

namespace UnitTesting
{
  TEST(TestGraphUtilities, TestConnectedComponents)
  {
    Gedim::GraphUtilities graphUtilities;

    // Create a graph given in the above diagram
    Gedim::GraphUtilities::Graph graph;
    const unsigned int graph_numVertices = 5;
    const unsigned int graph_numEdges = 5;
    Eigen::MatrixXi graph_connectivity(2, graph_numEdges);

    graph_connectivity.col(0)<< 1, 0;
    graph_connectivity.col(1)<< 0, 2;
    graph_connectivity.col(2)<< 2, 1;
    graph_connectivity.col(3)<< 0, 3;
    graph_connectivity.col(4)<< 3, 4;

    const std::vector<std::vector<unsigned int>> graph_adjacency = graphUtilities.ComputeGraphAdjacency(graph_numVertices,
                                                                                                        graph_connectivity);

    const std::vector<std::vector<unsigned int>> stronglyConnectedComponents = graphUtilities.ComputeStronglyConnectedComponents(graph_numVertices,
                                                                                                                                 graph_adjacency);

    ASSERT_EQ(3, stronglyConnectedComponents.size());
    ASSERT_EQ(std::vector<unsigned int>({ 0, 1, 2 }), stronglyConnectedComponents[0]);
    ASSERT_EQ(std::vector<unsigned int>({ 3 }), stronglyConnectedComponents[1]);
    ASSERT_EQ(std::vector<unsigned int>({ 4 }), stronglyConnectedComponents[2]);
  }
}

#endif // __TEST_METISUTILITIES_H
