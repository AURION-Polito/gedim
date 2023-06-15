#ifndef __METIS_UTILITIES_H
#define __METIS_UTILITIES_H

#include "Eigen/Eigen"
#include "IMeshDAO.hpp"
#include <vector>

namespace Gedim
{
  class MetisUtilities final
  {
    public:
      struct NetworkPartitionOptions final
      {
          enum struct PartitionTypes
          {
            Unknown = -1,
            CutBalancing = 0,
            VolBalancing = 1
          };

          PartitionTypes PartitionType = PartitionTypes::Unknown;
          unsigned int NumberOfParts = 0;
          unsigned int MasterWeight = 100; ///< 0 de-activated; 100 activated totally
      };

      struct MetisNetwork final
      {
          struct MetisAdjacency final
          {
              std::vector<unsigned int> Rows;
              std::vector<unsigned int> Cols;
          };

          MetisAdjacency Adjacency;
          std::vector<bool> EdgesConstrained;
          std::vector<unsigned int> NodesWeight = {};
          std::vector<unsigned int> EdgesWeight = {};
      };

      struct MeshToNetwork final
      {
          std::vector<unsigned int> EdgesMeshCellIndex;
          MetisNetwork Network;
      };

    public:
      MetisUtilities();
      ~MetisUtilities();

      MetisUtilities::MeshToNetwork Mesh3DToDualGraph(const IMeshDAO& mesh,
                                                      const std::vector<bool>& cell2DsConstrained = {},
                                                      const Eigen::SparseMatrix<unsigned int>& weights = Eigen::SparseMatrix<unsigned int>()) const;
      MetisUtilities::MeshToNetwork Mesh2DToDualGraph(const IMeshDAO& mesh,
                                                      const std::vector<bool>& cell1DsConstrained = {},
                                                      const Eigen::SparseMatrix<unsigned int>& weights = Eigen::SparseMatrix<unsigned int>()) const;
      MetisUtilities::MetisNetwork Mesh2DToGraph(const unsigned int& numVertices,
                                                 const Eigen::MatrixXi& edges,
                                                 const bool& undirectEdges) const;

      MetisNetwork::MetisAdjacency GraphAdjacencyToMetisAdjacency(const std::vector<std::vector<unsigned int>>& graphAdjacency) const;
      std::vector<std::vector<unsigned int>> MetisAdjacencyToGraphAdjacency(const MetisNetwork::MetisAdjacency& metisAdjacency) const;

      std::vector<unsigned int> NetworkPartition(const NetworkPartitionOptions& options,
                                                 const MetisNetwork& network) const;

      std::vector<unsigned int> PartitionCheckConstraints(const MetisNetwork& network,
                                                          const std::vector<unsigned int>& partitions) const;
      std::vector<unsigned int> PartitionCheckConnectedComponents(const MetisNetwork& network,
                                                                  const std::vector<unsigned int>& partitions) const;
  };
}


#endif // __METIS_UTILITIES_H

