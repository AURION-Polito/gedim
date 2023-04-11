#include "MetisUtilities.hpp"

#include "Macro.hpp"

#include "IOUtilities.hpp"

#if ENABLE_METIS == 1
#include "metis.h"
#endif

namespace Gedim
{
  // ***************************************************************************
  MetisUtilities::MetisUtilities()
  {
  }
  MetisUtilities::~MetisUtilities()
  {
  }
  // ***************************************************************************
  MetisUtilities::NetworkAdjacency MetisUtilities::GraphToAdjacency(const unsigned int& numVertices,
                                                                    const Eigen::MatrixXi& edges,
                                                                    const bool& undirectEdges) const
  {
    const unsigned int& numEdges = edges.cols();

    NetworkAdjacency adjacency;
    adjacency.AdjacencyRows.resize(numVertices + 1, 0);
    std::list<unsigned int> adjacencyCols;

    std::vector<std::list<unsigned int>> verticesConnections(numVertices);
    for (unsigned int e = 0; e < numEdges; e++)
    {
      const Eigen::VectorXi& edge = edges.col(e);
      verticesConnections[edge[0]].push_back(edge[1]);
      if (undirectEdges)
        verticesConnections[edge[1]].push_back(edge[0]);
    }

    for (unsigned int v = 0; v < numVertices; v++)
    {
      const std::list<unsigned int>& vertexConnections = verticesConnections[v];
      const unsigned int& numberConnections = vertexConnections.size();

      adjacency.AdjacencyRows[v + 1] = adjacency.AdjacencyRows[v] +
                                       numberConnections;

      for (const unsigned int& connection : vertexConnections)
        adjacencyCols.push_back(connection);
    }

    adjacency.AdjacencyCols = std::vector<unsigned int>(adjacencyCols.begin(),
                                                        adjacencyCols.end());

    return adjacency;
  }
  // ***************************************************************************
  std::vector<unsigned int> MetisUtilities::NetworkPartition(const NetworkPartitionOptions& options,
                                                             const NetworkAdjacency& adjacency) const
  {
    /// <ul>

    /// <li> Check the number of parts
    const unsigned int& numberElements = adjacency.AdjacencyRows.size() - 1;
    unsigned int numberParts = options.NumberOfParts;
    const unsigned int& masterWeight = options.MasterWeight;

    Output::Assert(numberElements > 0);

    std::vector<unsigned int> partition(numberElements, 0);
    Output::Assert(numberParts > 0);

    if (numberParts == 1)
    {
      if (masterWeight == 0)
        partition.resize(numberElements, 1);
      else
        partition.resize(numberElements, 0);

      return partition;
    }

#if ENABLE_METIS == 1
    unsigned int numberConnections = adjacency.AdjacencyCols.size();

    /// <li> Initialize METIS Partition
    rstatus_et metisResult = METIS_OK;
    idx_t nParts = numberParts;
    idx_t objval;
    idx_t nElements = numberElements;
    std::vector<idx_t> xadj(nElements + 1);
    std::vector<idx_t> adjncy(numberConnections);
    std::vector<idx_t> metisPartition(numberElements, 0);

    const NetworkPartitionOptions::PartitionTypes& metisDivisionType = options.PartitionType;

    idx_t nCon = 1;
    idx_t* vwgt = NULL, *adjwgt = NULL;
    idx_t* v_size = NULL;
    real_t* tpwgts = NULL, *ubvec = NULL;
    idx_t metisOptions[METIS_NOPTIONS];

    /// <li> Build adjncy and xadj for METIS partition
    memcpy(xadj.data(),
           adjacency.AdjacencyRows.data(),
           sizeof(idx_t) * (numberElements + 1));
    memcpy(adjncy.data(),
           adjacency.AdjacencyCols.data(),
           sizeof(idx_t) * numberConnections);

    METIS_SetDefaultOptions(metisOptions);

    switch (metisDivisionType)
    {
      case NetworkPartitionOptions::PartitionTypes::CutBalancing:
      {
        metisOptions[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

        if (options.NodeWeights.size() > 0)
        {
          vwgt = new idx_t[numberElements];
          memcpy(vwgt, options.NodeWeights.data(), sizeof(idx_t) * numberElements);
        }

        if (options.EdgeWeights.size() > 0)
        {
          adjwgt = new idx_t[numberConnections];
          memcpy(adjwgt, options.EdgeWeights.data(), sizeof(idx_t)*numberConnections);
        }
      }
        break;
      case NetworkPartitionOptions::PartitionTypes::VolBalancing:
      {
        metisOptions[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

        if (options.NodeWeights.size() > 0)
        {
          v_size = new idx_t[numberElements];
          memcpy(v_size, options.NodeWeights.data(), sizeof(idx_t) * numberElements);
        }
      }
        break;
      default:
        throw std::runtime_error("MetisDivisionType " +
                                 std::to_string(static_cast<unsigned int>(metisDivisionType)) +
                                 " not supported");
    }

    if (masterWeight > 0)
    {
      tpwgts = new real_t[numberParts];
      tpwgts[0] = (real_t)masterWeight / (real_t)numberParts / 100.0;

      for (unsigned int p = 1; p < numberParts; p++)
        tpwgts[p] = 1.0 / (numberParts - 1.0) * (1.0 - tpwgts[0]);
    }

    /// <li> Partiton
    metisResult = (rstatus_et)METIS_PartGraphKway(&nElements,
                                                  &nCon,
                                                  xadj.data(),
                                                  adjncy.data(),
                                                  vwgt,
                                                  v_size,
                                                  adjwgt,
                                                  &nParts,
                                                  tpwgts,
                                                  ubvec,
                                                  metisOptions,
                                                  &objval,
                                                  metisPartition.data());

    switch(metisResult)
    {
      case METIS_ERROR_INPUT:
        throw std::runtime_error("METIS failed due to erroneous inputs and/or options");
      case METIS_ERROR_MEMORY:
        throw std::runtime_error("METIS failed due to insufficient memory");
      case METIS_ERROR:
        throw std::runtime_error("METIS failed because of an unknown error");
      default:
        break;
    }

    partition.resize(numberElements);

    // Sometimes it happens METIS partition is not correct because of a METIS implementation error
    /// <li> Correct METIS partition
    if(objval == 0)
    {
      if (masterWeight == 0)
        std::fill_n(partition.data(), nElements, 1);
      else
        std::fill_n(partition.data(), nElements, 0);
    }
    else
    {
      if (masterWeight > 0)
        memcpy(partition.data(), metisPartition.data(), numberElements * sizeof(unsigned int));
      else
      {
        for (unsigned int p = 0; p < numberElements; p++)
          partition[p] = metisPartition[p] + 1;
      }
    }
    delete[] tpwgts; tpwgts = NULL;
    delete[] v_size; v_size = NULL;
    delete[] vwgt; vwgt = NULL;
    delete[] adjwgt; adjwgt = NULL;

#endif

    // Reordering partition vector
    unsigned int controlActivation=(masterWeight>0) ? 0 : 1;
    bool control =false;
    unsigned int loop=0;
    while(!control)
    {
      if(find(partition.begin(),partition.end(),controlActivation)==partition.end())
      {
        for (unsigned int p = 0; p < numberElements; p++)
        {
          if(controlActivation<partition[p])
            partition[p]=partition[p]-1;
        }
      }
      else{
        controlActivation++;
      }
      loop++;
      control = (loop==numberParts) ? true : false ;
    }

    return partition;

    /// </ul>
  }
  // ***************************************************************************
}
