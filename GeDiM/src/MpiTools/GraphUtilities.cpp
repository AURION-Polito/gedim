#include "GraphUtilities.hpp"
#include <queue>

namespace Gedim
{
// ***************************************************************************
void GraphUtilities::DepthFirstSearch(const unsigned int &vertex,
                                      const std::vector<std::vector<unsigned int>> &graphAdjacency,
                                      std::vector<bool> &visited,
                                      std::list<unsigned int> &visitedVertices) const
{
    // Mark the current node as visited and print it
    visited[vertex] = true;
    visitedVertices.push_back(vertex);

    // Recur for all the vertices adjacent to this vertex
    for (const unsigned int &adj_v : graphAdjacency[vertex])
    {
        if (visited[adj_v])
            continue;

        DepthFirstSearch(adj_v, graphAdjacency, visited, visitedVertices);
    }
}
// ***************************************************************************
std::vector<unsigned int> GraphUtilities::BreadthFirstSearch(const unsigned int &vertex,
                                                             const std::vector<std::vector<unsigned int>> &graphAdjacency) const
{
    unsigned int visitedVertex = vertex;
    std::list<unsigned int> visitedVertices;
    std::vector<bool> markedVertices(graphAdjacency.size(), false);
    std::queue<unsigned int> queue;

    markedVertices[visitedVertex] = true;
    queue.push(visitedVertex);

    while (!queue.empty())
    {
        visitedVertex = queue.front();
        visitedVertices.push_back(visitedVertex);
        queue.pop();

        for (const unsigned int adjacent : graphAdjacency[visitedVertex])
        {
            if (!markedVertices[adjacent])
            {
                markedVertices[adjacent] = true;
                queue.push(adjacent);
            }
        }
    }

    return std::vector<unsigned int>(visitedVertices.begin(), visitedVertices.end());
}
// ***************************************************************************
std::vector<std::vector<unsigned int>> GraphUtilities::ComputeAdjacencyTranspose(const std::vector<std::vector<unsigned int>> &graphAdjacency) const
{
    const unsigned int &graphNumVertices = graphAdjacency.size();
    std::vector<std::vector<unsigned int>> graphT_Adjacency(graphNumVertices);

    std::vector<std::list<unsigned int>> adj(graphNumVertices);
    for (unsigned int v = 0; v < graphNumVertices; v++)
    {
        for (const unsigned int &adj_v : graphAdjacency[v])
            adj[adj_v].push_back(v);
    }

    for (unsigned int v = 0; v < graphNumVertices; v++)
        graphT_Adjacency[v] = std::vector<unsigned int>(adj[v].begin(), adj[v].end());

    return graphT_Adjacency;
}
// ***************************************************************************
void GraphUtilities::FillOrder(const unsigned int &v,
                               const std::vector<std::vector<unsigned int>> &graphAdjacency,
                               std::vector<bool> &visited,
                               std::stack<int> &stack) const
{
    // Mark the current node as visited and print it
    visited[v] = true;

    // Recur for all the vertices adjacent to this vertex
    for (const unsigned int &adj_v : graphAdjacency[v])
    {
        if (visited[adj_v])
            continue;

        FillOrder(adj_v, graphAdjacency, visited, stack);
    }

    // All vertices reachable from v are processed by now, push v
    stack.push(v);
}
// ***************************************************************************
std::vector<std::vector<unsigned int>> GraphUtilities::ExtractSubGraph(const std::vector<std::vector<unsigned int>> &graphAdjacency,
                                                                       const std::unordered_map<unsigned int, unsigned int> &subGraphFilter) const
{
    const unsigned int &numFilteredVertices = subGraphFilter.size();
    std::vector<std::vector<unsigned int>> subGraph_Adjacency(numFilteredVertices);

    std::vector<std::list<unsigned int>> filtered_adj(numFilteredVertices);
    for (const auto &v : subGraphFilter)
    {
        const unsigned int original_v = v.first;
        const unsigned int filtered_v = v.second;

        for (const unsigned int &original_adj_v : graphAdjacency[original_v])
        {
            const std::unordered_map<unsigned int, unsigned int>::const_iterator adj_v = subGraphFilter.find(original_adj_v);
            if (adj_v == subGraphFilter.end())
                continue;

            const unsigned int filtered_adj_v = adj_v->second;
            filtered_adj[filtered_v].push_back(filtered_adj_v);
        }
    }

    for (unsigned int v = 0; v < numFilteredVertices; v++)
        subGraph_Adjacency[v] = std::vector<unsigned int>(filtered_adj[v].begin(), filtered_adj[v].end());

    return subGraph_Adjacency;
}
// ***************************************************************************
GraphUtilities::GraphAdjacencyData GraphUtilities::GraphConnectivityToGraphAdjacency(const unsigned int &graphNumVertices,
                                                                                     const Eigen::MatrixXi &graphConnectivity,
                                                                                     const bool &directEdges) const
{
    struct AdjacencyData
    {
        unsigned int VertexIndex;
        unsigned int EdgeIndex;
    };

    std::vector<std::list<AdjacencyData>> adj(graphNumVertices);
    for (unsigned int e = 0; e < graphConnectivity.cols(); e++)
    {
        adj[graphConnectivity(0, e)].push_back({static_cast<unsigned int>(graphConnectivity(1, e)), e});

        if (!directEdges)
        {
            adj[graphConnectivity(1, e)].push_back({static_cast<unsigned int>(graphConnectivity(0, e)), e});
        }
    }

    GraphAdjacencyData result;
    result.GraphAdjacencyVertices.resize(graphNumVertices);
    result.GraphAdjacencyVerticesMap.resize(graphNumVertices);
    result.GraphAdjacencyEdges.resize(graphNumVertices);
    result.GraphAdjacencyEdgesMap.resize(graphNumVertices);

    for (unsigned int v = 0; v < graphNumVertices; v++)
    {
        result.GraphAdjacencyVertices[v].reserve(adj[v].size());
        result.GraphAdjacencyEdges[v].reserve(adj[v].size());

        for (const AdjacencyData &ad_data : adj[v])
        {
            result.GraphAdjacencyVerticesMap[v].insert(
                std::make_pair(ad_data.VertexIndex, result.GraphAdjacencyVertices[v].size()));
            result.GraphAdjacencyVertices[v].push_back(ad_data.VertexIndex);

            result.GraphAdjacencyEdgesMap[v].insert(std::make_pair(ad_data.EdgeIndex, result.GraphAdjacencyEdges[v].size()));
            result.GraphAdjacencyEdges[v].push_back(ad_data.EdgeIndex);
        }

        result.GraphAdjacencyVertices[v].shrink_to_fit();
        result.GraphAdjacencyEdges[v].shrink_to_fit();
    }

    return result;
}
// ***************************************************************************
Eigen::MatrixXi GraphUtilities::GraphAdjacencyToGraphConnectivity(const unsigned int &graphNumEdges,
                                                                  const std::vector<std::vector<unsigned int>> &graphAdjacency) const
{

    Eigen::MatrixXi graphConnectivity(2, graphNumEdges);

    unsigned int e = 0;
    for (unsigned int v = 0; v < graphAdjacency.size(); v++)
    {
        for (const unsigned int &adj_v : graphAdjacency[v])
            graphConnectivity.col(e++) << v, adj_v;
    }

    return graphConnectivity;
}
// ***************************************************************************
std::vector<std::vector<unsigned int>> GraphUtilities::ComputeStronglyConnectedComponents(
    const std::vector<std::vector<unsigned int>> &graphAdjacency) const
{
    std::stack<int> stack;

    const unsigned int &graphNumVertices = graphAdjacency.size();

    // Mark all the vertices as not visited (For first DFS)
    std::vector<bool> visited(graphNumVertices, false);

    // Fill vertices in stack according to their finishing times
    for (unsigned int i = 0; i < graphNumVertices; i++)
    {
        if (visited[i])
            continue;

        FillOrder(i, graphAdjacency, visited, stack);
    }

    // Create a reversed graph
    const std::vector<std::vector<unsigned int>> graphT_Adjacency = ComputeAdjacencyTranspose(graphAdjacency);

    // Mark all the vertices as not visited (For second DFS)
    for (unsigned int i = 0; i < graphNumVertices; i++)
        visited[i] = false;

    std::list<std::vector<unsigned int>> connectedComponents;

    // Now process all vertices in order defined by Stack
    while (!stack.empty())
    {
        // Pop a vertex from stack
        int v = stack.top();
        stack.pop();

        // Print Strongly connected component of the popped vertex
        if (visited[v])
            continue;

        std::list<unsigned int> connectedComponent;
        DepthFirstSearch(v, graphT_Adjacency, visited, connectedComponent);

        connectedComponents.push_back(std::vector<unsigned int>(connectedComponent.begin(), connectedComponent.end()));
    }

    return std::vector<std::vector<unsigned int>>(connectedComponents.begin(), connectedComponents.end());
}
// ***************************************************************************
} // namespace Gedim
