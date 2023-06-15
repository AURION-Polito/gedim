#ifndef __GraphUtilities_H
#define __GraphUtilities_H

#include "Eigen/Eigen"
#include <iostream>
#include <list>
#include <stack>
#include <unordered_map>
#include <vector>

namespace Gedim
{
  class GraphUtilities final
  {
    private:

      /// Fills Stack with vertices (in increasing order of finishing
      /// times). The top element of stack has the maximum finishing
      /// time
      void FillOrder(const unsigned int& v,
                     const std::vector<std::vector<unsigned int> >& graphAdjacency,
                     std::vector<bool>& visited,
                     std::stack<int>& stack) const;

    public:
      GraphUtilities() { };
      ~GraphUtilities() { };

      /// \param graphAdjacency the original graph adiajency
      /// \param verticesFilter the vertices of the sub-graph, (original_vertex_position, filtered_vertex_position)
      /// \return the resulting filtered graph
      std::vector<std::vector<unsigned int>> ExtractSubGraph(const std::vector<std::vector<unsigned int>>& graphAdjacency,
                                                             const std::unordered_map<unsigned int, unsigned int>& subGraphFilter) const;

      std::vector<std::vector<unsigned int>> GraphConnectivityToGraphAdjacency(const unsigned int& graphNumVertices,
                                                                               const Eigen::MatrixXi& graphConnectivity) const;

      Eigen::MatrixXi GraphAdjacencyToGraphConnectivity(const unsigned int& graphNumEdges,
                                                        const std::vector<std::vector<unsigned int>>& graphAdjacency) const;

      /// \brief Compute the Strongly Connected Components of a direct graph
      /// \param numVertices the numGraphVertices
      /// \param adjacency the graph adiacency matrix, size 2 x numGraphEdges
      /// \return the Strongly Connected Components
      std::vector<std::vector<unsigned int>> ComputeStronglyConnectedComponents(const std::vector<std::vector<unsigned int>>& graphAdjacency) const;

      /// A recursive function to DFS starting from v
      void DepthFirstSearch(const unsigned int& vertex,
                            const std::vector<std::vector<unsigned int>>& graphAdjacency,
                            std::vector<bool>& visited,
                            std::list<unsigned int>& visitedVertices) const;

      /// \return the reverse (or transpose) of a graph
      std::vector<std::vector<unsigned int>> ComputeAdjacencyTranspose(const std::vector<std::vector<unsigned int>>& graphAdjacency) const;
  };

}


#endif // __GraphUtilities_H

