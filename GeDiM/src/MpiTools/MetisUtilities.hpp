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

      struct MeshToNetwork final
      {
          std::vector<unsigned int> EdgesMeshCellIndex;
          Network MetisNetwork;
      };

    public:
      MetisUtilities();
      ~MetisUtilities();

      MetisUtilities::MeshToNetwork Mesh3DToDualGraph(const IMeshDAO& mesh,
                                                      const std::vector<bool>& facesConstrained = {},
                                                      const Eigen::SparseMatrix<unsigned int>& weights = Eigen::SparseMatrix<unsigned int>()) const;


      MetisUtilities::MeshToNetwork Mesh2DToDualGraph(const IMeshDAO& mesh,
                                                      const std::vector<bool>& edgeConstrained = {},
                                                      const Eigen::SparseMatrix<unsigned int>& weights = Eigen::SparseMatrix<unsigned int>()) const;

      MetisUtilities::Network Mesh2DToGraph(const unsigned int& numVertices,
                                            const Eigen::MatrixXi& edges,
                                            const bool& undirectEdges) const;
      Eigen::MatrixXi GraphToConnectivityMatrix(const Network& network) const;

      std::vector<unsigned int> NetworkPartition(const NetworkPartitionOptions& options,
                                                 const Network& network) const;

      std::vector<unsigned int> Mesh3DPartitionCheckConstraints(const IMeshDAO& mesh,
                                                                const std::vector<bool>& facesConstrained,
                                                                const std::vector<unsigned int> partition) const;
  };
}


#endif // __METIS_UTILITIES_H

