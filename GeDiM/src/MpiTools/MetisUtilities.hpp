#ifndef __METIS_UTILITIES_H
#define __METIS_UTILITIES_H

#include "Eigen/Eigen"
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
          std::vector<unsigned int> NodeWeights = {};
          std::vector<unsigned int> EdgeWeights = {};
      };

      struct NetworkAdjacency final
      {
          std::vector<unsigned int> AdjacencyRows;
          std::vector<unsigned int> AdjacencyCols;
      };

    public:
      MetisUtilities();
      ~MetisUtilities();

      MetisUtilities::NetworkAdjacency GraphToAdjacency(const unsigned int& numVertices,
                                                        const Eigen::MatrixXi& edges,
                                                        const bool& undirectEdges) const;

      std::vector<unsigned int> NetworkPartition(const NetworkPartitionOptions& options,
                                                 const NetworkAdjacency& adjacency) const;
  };
}


#endif // __METIS_UTILITIES_H

