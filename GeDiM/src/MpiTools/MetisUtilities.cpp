#include "MetisUtilities.hpp"

#include "GraphUtilities.hpp"
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
  MetisUtilities::MeshToNetwork MetisUtilities::Mesh3DToDualGraph(const IMeshDAO& mesh,
                                                                  const std::vector<bool>& cell2DsConstrained,
                                                                  const Eigen::SparseMatrix<unsigned int>& weights) const
  {
    MeshToNetwork meshToNetwork;

    MetisNetwork& network = meshToNetwork.Network;


    struct Connection final
    {
        bool Constrained;
        unsigned int Cell3DIndex;
        unsigned int Cell2DIndex;
    };


    const unsigned int numVertices = mesh.Cell3DTotalNumber();
    std::set<std::pair<unsigned int, unsigned int>> constraints;

    std::vector<std::list<Connection>> verticesConnections(numVertices);
    const bool checkConstraints = (cell2DsConstrained.size() == mesh.Cell2DTotalNumber());

    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); f++)
    {
      for (unsigned int n1 = 0; n1 < mesh.Cell2DNumberNeighbourCell3D(f); n1++)
      {
        if (!mesh.Cell2DHasNeighbourCell3D(f, n1))
          continue;

        const unsigned int cell3Dn1 = mesh.Cell2DNeighbourCell3D(f, n1);

        for (unsigned int n2 = n1 + 1; n2 < mesh.Cell2DNumberNeighbourCell3D(f); n2++)
        {
          if (!mesh.Cell2DHasNeighbourCell3D(f, n2))
            continue;

          const unsigned int cell3Dn2 = mesh.Cell2DNeighbourCell3D(f, n2);
          const bool cell1DConstrained = checkConstraints ? cell2DsConstrained[f] :
                                                            false;

          verticesConnections[cell3Dn1].push_back(Connection { cell1DConstrained, cell3Dn2, f });
          verticesConnections[cell3Dn2].push_back(Connection { cell1DConstrained, cell3Dn1, f });

          if (cell1DConstrained)
          {
            constraints.insert(std::make_pair(cell3Dn1, cell3Dn2));
            constraints.insert(std::make_pair(cell3Dn2, cell3Dn1));
          }
        }
      }
    }

    std::list<unsigned int> adjacencyCols;
    std::list<unsigned int> adjacencyColsMeshCellIndex;
    std::list<bool> adjacencyColsConstrained;

    network.Adjacency.Rows.resize(numVertices + 1, 0);
    for (unsigned int v = 0; v < numVertices; v++)
    {
      const std::list<Connection>& vertexConnections = verticesConnections[v];
      const unsigned int& numberConnections = vertexConnections.size();

      network.Adjacency.Rows[v + 1] = network.Adjacency.Rows[v] +
                                      numberConnections;

      for (const Connection& connection : vertexConnections)
      {
        adjacencyCols.push_back(connection.Cell3DIndex);
        adjacencyColsMeshCellIndex.push_back(connection.Cell2DIndex);
        adjacencyColsConstrained.push_back(connection.Constrained);
      }
    }

    network.Adjacency.Cols = std::vector<unsigned int>(adjacencyCols.begin(),
                                                       adjacencyCols.end());
    meshToNetwork.EdgesMeshCellIndex = std::vector<unsigned int>(adjacencyColsMeshCellIndex.begin(),
                                                                 adjacencyColsMeshCellIndex.end());
    network.EdgesConstrained = std::vector<bool>(adjacencyColsConstrained.begin(),
                                                 adjacencyColsConstrained.end());

    const unsigned int& numEdges = network.Adjacency.Cols.size();
    network.EdgesWeight.resize(numEdges, 1);

    int counter = 0;
    for (unsigned int v = 0; v < numVertices; v++)
    {
      for (const Connection& connection : verticesConnections[v])
      {
        unsigned int weight = (weights.size() > 0) ? weights.coeff(v, connection.Cell3DIndex) : 1;
        network.EdgesWeight[counter++] =
            constraints.find(std::make_pair(v, connection.Cell3DIndex)) != constraints.end() ?
                                                                             1 :
                                                                             weight + numEdges;
      }
    }

    return meshToNetwork;
  }
  // ***************************************************************************
  MetisUtilities::MetisNetwork MetisUtilities::Mesh2DToGraph(const unsigned int& numVertices,
                                                             const Eigen::MatrixXi& edges,
                                                             const bool& undirectEdges) const
  {
    const unsigned int& numEdges = edges.cols();

    MetisNetwork network;
    std::list<unsigned int> adjacencyCols;

    std::vector<std::list<unsigned int>> verticesConnections(numVertices);
    for (unsigned int e = 0; e < numEdges; e++)
    {
      const Eigen::VectorXi& edge = edges.col(e);
      verticesConnections[edge[0]].push_back(edge[1]);
      if (undirectEdges)
        verticesConnections[edge[1]].push_back(edge[0]);
    }

    network.Adjacency.Rows.resize(numVertices + 1, 0);
    for (unsigned int v = 0; v < numVertices; v++)
    {
      const std::list<unsigned int>& vertexConnections = verticesConnections[v];
      const unsigned int& numberConnections = vertexConnections.size();

      network.Adjacency.Rows[v + 1] = network.Adjacency.Rows[v] +
                                      numberConnections;

      for (const unsigned int& connection : vertexConnections)
        adjacencyCols.push_back(connection);
    }

    network.Adjacency.Cols = std::vector<unsigned int>(adjacencyCols.begin(),
                                                       adjacencyCols.end());

    return network;
  }
  // ***************************************************************************
  MetisUtilities::MetisNetwork::MetisAdjacency MetisUtilities::GraphAdjacencyToMetisAdjacency(const std::vector<std::vector<unsigned int>>& graphAdjacency) const
  {
    MetisNetwork::MetisAdjacency metisAdjacency;

    const unsigned int& numVertices = graphAdjacency.size();

    metisAdjacency.Rows.resize(numVertices + 1, 0);
    std::list<unsigned int> adjacencyCols;

    for (unsigned int v = 0; v < numVertices; v++)
    {
      metisAdjacency.Rows[v + 1] = metisAdjacency.Rows[v] + graphAdjacency[v].size();

      for (const unsigned int& adj_v : graphAdjacency[v])
        adjacencyCols.push_back(adj_v);
    }

    metisAdjacency.Cols = std::vector<unsigned int>(adjacencyCols.begin(),
                                                    adjacencyCols.end());

    return metisAdjacency;
  }
  // ***************************************************************************
  std::vector<std::vector<unsigned int>> MetisUtilities::MetisAdjacencyToGraphAdjacency(const MetisNetwork::MetisAdjacency& metisAdjacency) const
  {
    const unsigned int& numVertices = metisAdjacency.Rows.size() - 1;
    std::vector<std::vector<unsigned int>> graphAdiacency(numVertices);

    for (unsigned int v = 0; v < numVertices; v++)
    {
      const unsigned int numConnections = metisAdjacency.Rows[v + 1] - metisAdjacency.Rows[v];
      graphAdiacency[v].resize(numConnections);
      for (unsigned int c = 0; c < numConnections; c++)
        graphAdiacency[v][c] = metisAdjacency.Cols.at(metisAdjacency.Rows[v] + c);
    }

    return graphAdiacency;
  }
  // ***************************************************************************
  MetisUtilities::MeshToNetwork MetisUtilities::Mesh2DToDualGraph(const IMeshDAO& mesh,
                                                                  const std::vector<bool>& cell1DsConstrained,
                                                                  const Eigen::SparseMatrix<unsigned int>& weights) const
  {
    MeshToNetwork meshToNetwork;

    MetisNetwork& network = meshToNetwork.Network;

    const unsigned int numVertices = mesh.Cell2DTotalNumber();
    std::set<std::pair<unsigned int, unsigned int>> constraints;

    struct Connection final
    {
        bool Constrained;
        unsigned int Cell2DIndex;
        unsigned int Cell1DIndex;
    };

    std::vector<std::list<Connection>> verticesConnections(numVertices);
    const bool checkConstraints = (cell1DsConstrained.size() == mesh.Cell1DTotalNumber());

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
    {
      for (unsigned int n1 = 0; n1 < mesh.Cell1DNumberNeighbourCell2D(e); n1++)
      {
        if (!mesh.Cell1DHasNeighbourCell2D(e, n1))
          continue;

        const unsigned int cell2Dn1 = mesh.Cell1DNeighbourCell2D(e, n1);

        for (unsigned int n2 = n1 + 1; n2 < mesh.Cell1DNumberNeighbourCell2D(e); n2++)
        {
          if (!mesh.Cell1DHasNeighbourCell2D(e, n2))
            continue;

          const unsigned int cell2Dn2 = mesh.Cell1DNeighbourCell2D(e, n2);
          const bool cell1DConstrained = checkConstraints ? cell1DsConstrained[e] :
                                                            false;

          verticesConnections[cell2Dn1].push_back(Connection { cell1DConstrained, cell2Dn2, e });
          verticesConnections[cell2Dn2].push_back(Connection { cell1DConstrained, cell2Dn1, e });

          if (cell1DConstrained)
          {
            constraints.insert(std::make_pair(cell2Dn1, cell2Dn2));
            constraints.insert(std::make_pair(cell2Dn2, cell2Dn1));
          }
        }
      }
    }

    std::list<unsigned int> adjacencyCols;
    std::list<unsigned int> adjacencyColsMeshCellIndex;
    std::list<bool> adjacencyColsConstrained;
    network.Adjacency.Rows.resize(numVertices + 1, 0);

    for (unsigned int v = 0; v < numVertices; v++)
    {
      const std::list<Connection>& vertexConnections = verticesConnections[v];
      const unsigned int& numberConnections = vertexConnections.size();

      network.Adjacency.Rows[v + 1] = network.Adjacency.Rows[v] +
                                      numberConnections;

      for (const Connection& connection : vertexConnections)
      {
        adjacencyCols.push_back(connection.Cell2DIndex);
        adjacencyColsMeshCellIndex.push_back(connection.Cell1DIndex);
        adjacencyColsConstrained.push_back(connection.Constrained);
      }
    }

    network.Adjacency.Cols = std::vector<unsigned int>(adjacencyCols.begin(),
                                                       adjacencyCols.end());
    meshToNetwork.EdgesMeshCellIndex = std::vector<unsigned int>(adjacencyColsMeshCellIndex.begin(),
                                                                 adjacencyColsMeshCellIndex.end());
    network.EdgesConstrained = std::vector<bool>(adjacencyColsConstrained.begin(),
                                                 adjacencyColsConstrained.end());

    const unsigned int& numEdges = network.Adjacency.Cols.size();
    network.EdgesWeight.resize(numEdges, 1);

    int counter = 0;
    for (unsigned int v = 0; v < numVertices; v++)
    {
      for (const Connection& connection : verticesConnections[v])
      {
        unsigned int weight = (weights.size() > 0) ? weights.coeff(v, connection.Cell2DIndex) : 1;
        network.EdgesWeight[counter++] =
            constraints.find(std::make_pair(v, connection.Cell2DIndex)) != constraints.end() ?
                                                                             1 :
                                                                             weight + numEdges;
      }
    }

    return meshToNetwork;
  }
  // ***************************************************************************
  std::vector<unsigned int> MetisUtilities::NetworkPartition(const NetworkPartitionOptions& options,
                                                             const MetisNetwork& network) const
  {
    /// <ul>

    /// <li> Check the number of parts
    const unsigned int& numberElements = network.Adjacency.Rows.size() - 1;
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
    unsigned int numberConnections = network.Adjacency.Cols.size();

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
           network.Adjacency.Rows.data(),
           sizeof(idx_t) * (numberElements + 1));
    memcpy(adjncy.data(),
           network.Adjacency.Cols.data(),
           sizeof(idx_t) * numberConnections);

    METIS_SetDefaultOptions(metisOptions);

    switch (metisDivisionType)
    {
      case NetworkPartitionOptions::PartitionTypes::CutBalancing:
      {
        metisOptions[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

        if (network.NodesWeight.size() > 0)
        {
          vwgt = new idx_t[numberElements];
          memcpy(vwgt, network.NodesWeight.data(), sizeof(idx_t) * numberElements);
        }

        if (network.EdgesWeight.size() > 0)
        {
          adjwgt = new idx_t[numberConnections];
          memcpy(adjwgt, network.EdgesWeight.data(), sizeof(idx_t)*numberConnections);
        }
      }
        break;
      case NetworkPartitionOptions::PartitionTypes::VolBalancing:
      {
        metisOptions[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

        if (network.NodesWeight.size() > 0)
        {
          v_size = new idx_t[numberElements];
          memcpy(v_size, network.NodesWeight.data(), sizeof(idx_t) * numberElements);
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
    while (!control)
    {
      if (find(partition.begin(),partition.end(),controlActivation)==partition.end())
      {
        for (unsigned int p = 0; p < numberElements; p++)
        {
          if(controlActivation<partition[p])
            partition[p] = partition[p] - 1;
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
  std::vector<unsigned int> MetisUtilities::PartitionCheckConstraints(const MetisNetwork& network,
                                                                      const std::vector<unsigned int>& partitions) const
  {
    if (network.EdgesConstrained.size() == 0)
      return partitions;

    std::vector<unsigned int> fixedPartition = partitions;
    unsigned int numMaxPartition = *std::max_element(begin(fixedPartition), end(fixedPartition));

    const unsigned int& numVertices = network.Adjacency.Rows.size() - 1;

    unsigned int edge = 0;
    for (unsigned int v = 0; v < numVertices; v++)
    {
      const unsigned int numConnections = network.Adjacency.Rows[v + 1] - network.Adjacency.Rows[v];

      for (unsigned int c = 0; c < numConnections; c++)
      {
        if (!network.EdgesConstrained[edge++])
          continue;

        const unsigned int& adj_v = network.Adjacency.Cols.at(network.Adjacency.Rows[v] + c);
        if (fixedPartition.at(v) !=
            fixedPartition.at(adj_v))
          continue;

        fixedPartition[adj_v] = ++numMaxPartition;
      }
    }

    return fixedPartition;
  }
  // ***************************************************************************
  std::vector<unsigned int> MetisUtilities::PartitionCheckConnectedComponents(const MetisNetwork& network,
                                                                              const std::vector<unsigned int>& partitions) const
  {
    std::vector<unsigned int> fixedPartition = partitions;
    GraphUtilities graphUtilities;

    const unsigned int& numVertices = network.Adjacency.Rows.size() - 1;

    unsigned int numMaxPartition = *std::max_element(begin(fixedPartition), end(fixedPartition));
    const unsigned int numPartitions = numMaxPartition + 1;

    std::vector<std::list<unsigned int>> subGraphs_vertices(numPartitions);
    std::vector<std::unordered_map<unsigned int, unsigned int>> subGraphs_filter(numPartitions);
    for (unsigned int v = 0; v < numVertices; v++)
    {
      const unsigned int& partion = partitions[v];

      subGraphs_vertices[partion].push_back(v);
      subGraphs_filter[partion].insert(std::make_pair(v,
                                                      subGraphs_filter[partion].size()));
    }

    const std::vector<std::vector<unsigned int>> graph_adjacency = MetisAdjacencyToGraphAdjacency(network.Adjacency);
    for (unsigned int p = 0; p < numPartitions; p++)
    {
      const std::vector<std::vector<unsigned int>> subGraph_adjacency = graphUtilities.ExtractSubGraph(graph_adjacency,
                                                                                                       subGraphs_filter[p]);

      const std::vector<std::vector<unsigned int>> connectedComponents = graphUtilities.ComputeStronglyConnectedComponents(subGraph_adjacency);

      if (connectedComponents.size() < 2)
        continue;

      std::vector<unsigned int> subGraph_vertices(subGraphs_vertices[p].begin(),
                                                  subGraphs_vertices[p].end());

      for (unsigned int cc = 1; cc < connectedComponents.size(); cc++)
      {
        const unsigned int new_partition = ++numMaxPartition;

        for (unsigned int v = 0; v < connectedComponents[cc].size(); v++)
        {
          const unsigned int& original_v = subGraph_vertices.at(connectedComponents[cc][v]);
          fixedPartition[original_v] = new_partition;
        }
      }
    }

    return fixedPartition;
  }
  // ***************************************************************************
}
