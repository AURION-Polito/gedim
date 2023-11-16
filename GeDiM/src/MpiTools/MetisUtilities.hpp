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

          enum struct CoarseningSchemes
          {
            Default = -1, // METIS_CTYPE_SHEM
            METIS_CTYPE_RM = 0,
            METIS_CTYPE_SHEM = 1
          };

          enum struct InitialPartitioningSchemes
          {
            Default = -1, // METIS_IPTYPE_METISRB
            METIS_IPTYPE_GROW = 0,
            METIS_IPTYPE_RANDOM = 1,
            METIS_IPTYPE_EDGE = 2,
            METIS_IPTYPE_NODE = 3,
            METIS_IPTYPE_METISRB = 4
          };

          enum struct RefinementSchemes
          {
            Default = -1, // METIS_RTYPE_GREEDY
            METIS_RTYPE_FM = 0,
            METIS_RTYPE_GREEDY = 1,
            METIS_RTYPE_SEP2SIDED = 2,
            METIS_RTYPE_SEP1SIDED = 3
          };

          enum struct DebugLevels
          {
            None = -1,
            METIS_DBG_INFO       = 1,       /*!< Shows various diagnostic messages */
            METIS_DBG_TIME       = 2,       /*!< Perform timing analysis */
            METIS_DBG_COARSEN    = 4,	  /*!< Show the coarsening progress */
            METIS_DBG_REFINE     = 8,	  /*!< Show the refinement progress */
            METIS_DBG_IPART      = 16, 	  /*!< Show info on initial partitioning */
            METIS_DBG_MOVEINFO   = 32, 	  /*!< Show info on vertex moves during refinement */
            METIS_DBG_SEPINFO    = 64, 	  /*!< Show info on vertex moves during sep refinement */
            METIS_DBG_CONNINFO   = 128,     /*!< Show info on minimization of subdomain connectivity */
            METIS_DBG_CONTIGINFO = 256,     /*!< Show info on elimination of connected components */
            METIS_DBG_MEMORY     = 2048     /*!< Show info related to wspace allocation */
          };

          PartitionTypes PartitionType = PartitionTypes::Unknown;
          CoarseningSchemes CoarseningSchema = CoarseningSchemes::Default;
          RefinementSchemes RefinementSchema = RefinementSchemes::Default;
          InitialPartitioningSchemes InitialPartitioningSchema = InitialPartitioningSchemes::Default;
          DebugLevels DebugLevel = DebugLevels::None;
          unsigned int NumberOfParts = 0;
          unsigned int NumberRefinementIterations = 10;
          bool ContigousPartitions = true;
          bool CompressGraph = false;
          bool MinimizeConnectivity = false;
          unsigned int RandomSeed = -1;
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

