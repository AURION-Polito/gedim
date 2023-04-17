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

      struct Network final
      {
          std::vector<unsigned int> AdjacencyRows;
          std::vector<unsigned int> AdjacencyCols;
          std::vector<unsigned int> NodeWeights = {};
          std::vector<unsigned int> EdgeWeights = {};
      };

    public:
      MetisUtilities();
      ~MetisUtilities();

      MetisUtilities::Network Mesh2DToDualGraph(const IMeshDAO& mesh,
                                                const Eigen::MatrixXd& weights) const;

      MetisUtilities::Network Mesh2DToGraph(const unsigned int& numVertices,
                                            const Eigen::MatrixXi& edges,
                                            const bool& undirectEdges) const;
      Eigen::MatrixXi GraphToConnectivityMatrix(const Network& network) const;

      std::vector<unsigned int> NetworkPartition(const NetworkPartitionOptions& options,
                                                 const Network& network) const;
  };
}


#endif // __METIS_UTILITIES_H

