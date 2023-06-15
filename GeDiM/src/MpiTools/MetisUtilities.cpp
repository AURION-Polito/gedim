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
  MetisUtilities::MeshToNetwork MetisUtilities::Mesh3DToDualGraph(const IMeshDAO& mesh,
                                                                  const std::vector<bool>& facesConstrained,
                                                                  const Eigen::SparseMatrix<unsigned int>& weights) const
  {
    MeshToNetwork meshToNetwork;

    MetisNetwork& network = meshToNetwork.Network;

    const unsigned int numVertices = mesh.Cell3DTotalNumber();
    network.Adjacency.Rows.resize(numVertices + 1, 0);
    std::list<unsigned int> adjacencyCols;
    std::list<unsigned int> adjacencyColsCellIndex;

    std::set<std::pair<unsigned int, unsigned int>> constraints;

    struct Connection final
    {
        unsigned int Cell3DIndex;
        unsigned int Cell2DIndex;
    };

    std::vector<std::list<Connection>> verticesConnections(numVertices);
    const bool checkConstraints = (facesConstrained.size() == mesh.Cell2DTotalNumber());

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

          verticesConnections[cell3Dn1].push_back(Connection { cell3Dn2, f });
          verticesConnections[cell3Dn2].push_back(Connection { cell3Dn1, f });

          if (checkConstraints &&
              facesConstrained[f])
          {
            constraints.insert(std::make_pair(cell3Dn1, cell3Dn2));
            constraints.insert(std::make_pair(cell3Dn2, cell3Dn1));
          }
        }
      }
    }

    for (unsigned int v = 0; v < numVertices; v++)
    {
      const std::list<Connection>& vertexConnections = verticesConnections[v];
      const unsigned int& numberConnections = vertexConnections.size();

      network.Adjacency.Rows[v + 1] = network.Adjacency.Rows[v] +
                                      numberConnections;

      for (const Connection& connection : vertexConnections)
      {
        adjacencyCols.push_back(connection.Cell3DIndex);
        adjacencyColsCellIndex.push_back(connection.Cell2DIndex);
      }
    }

    network.Adjacency.Cols = std::vector<unsigned int>(adjacencyCols.begin(),
                                                       adjacencyCols.end());
    meshToNetwork.EdgesMeshCellIndex = std::vector<unsigned int>(adjacencyColsCellIndex.begin(),
                                                                 adjacencyColsCellIndex.end());

    const unsigned int& numEdges = network.Adjacency.Cols.size();
    network.EdgeWeights.resize(numEdges, 1);

    int counter = 0;
    for (unsigned int v = 0; v < numVertices; v++)
    {
      const std::list<Connection>& vertexConnections = verticesConnections[v];
      for (const Connection& connection : vertexConnections)
      {
        unsigned int weight = (weights.size() > 0) ? weights.coeff(v, connection.Cell3DIndex) : 1;
        network.EdgeWeights[counter++] =
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
    network.Adjacency.Rows.resize(numVertices + 1, 0);
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
                                                                  const std::vector<bool>& edgesConstrained,
                                                                  const Eigen::SparseMatrix<unsigned int>& weights) const
  {
    MeshToNetwork meshToNetwork;

    MetisNetwork& network = meshToNetwork.Network;

    const unsigned int numVertices = mesh.Cell2DTotalNumber();
    network.Adjacency.Rows.resize(numVertices + 1, 0);
    std::list<unsigned int> adjacencyCols;
    std::list<unsigned int> adjacencyColsCellIndex;

    std::set<std::pair<unsigned int, unsigned int>> constraints;

    struct Connection final
    {
        unsigned int Cell2DIndex;
        unsigned int Cell1DIndex;
    };

    std::vector<std::list<Connection>> verticesConnections(numVertices);
    const bool checkConstraints = (edgesConstrained.size() == mesh.Cell1DTotalNumber());

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

          verticesConnections[cell2Dn1].push_back(Connection { cell2Dn2, e });
          verticesConnections[cell2Dn2].push_back(Connection { cell2Dn1, e });

          if (checkConstraints &&
              edgesConstrained[e])
          {
            constraints.insert(std::make_pair(cell2Dn1, cell2Dn2));
            constraints.insert(std::make_pair(cell2Dn2, cell2Dn1));
          }
        }
      }
    }

    for (unsigned int v = 0; v < numVertices; v++)
    {
      const std::list<Connection>& vertexConnections = verticesConnections[v];
      const unsigned int& numberConnections = vertexConnections.size();

      network.Adjacency.Rows[v + 1] = network.Adjacency.Rows[v] +
                                      numberConnections;

      for (const Connection& connection : vertexConnections)
      {
        adjacencyCols.push_back(connection.Cell2DIndex);
        adjacencyColsCellIndex.push_back(connection.Cell1DIndex);
      }
    }

    network.Adjacency.Cols = std::vector<unsigned int>(adjacencyCols.begin(),
                                                       adjacencyCols.end());
    meshToNetwork.EdgesMeshCellIndex = std::vector<unsigned int>(adjacencyColsCellIndex.begin(),
                                                                 adjacencyColsCellIndex.end());

    const unsigned int& numEdges = network.Adjacency.Cols.size();
    network.EdgeWeights.resize(numEdges, 1);

    int counter = 0;
    for (unsigned int v = 0; v < numVertices; v++)
    {
      const std::list<Connection>& vertexConnections = verticesConnections[v];
      for (const Connection& connection : vertexConnections)
      {
        unsigned int weight = (weights.size() > 0) ? weights.coeff(v, connection.Cell2DIndex) : 1;
        network.EdgeWeights[counter++] =
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

        if (network.NodeWeights.size() > 0)
        {
          vwgt = new idx_t[numberElements];
          memcpy(vwgt, network.NodeWeights.data(), sizeof(idx_t) * numberElements);
        }

        if (network.EdgeWeights.size() > 0)
        {
          adjwgt = new idx_t[numberConnections];
          memcpy(adjwgt, network.EdgeWeights.data(), sizeof(idx_t)*numberConnections);
        }
      }
        break;
      case NetworkPartitionOptions::PartitionTypes::VolBalancing:
      {
        metisOptions[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

        if (network.NodeWeights.size() > 0)
        {
          v_size = new idx_t[numberElements];
          memcpy(v_size, network.NodeWeights.data(), sizeof(idx_t) * numberElements);
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
  std::vector<unsigned int> MetisUtilities::Mesh2DPartitionCheckConstraints(const IMeshDAO& mesh,
                                                                            const std::vector<bool>& edgesConstrained,
                                                                            const std::vector<unsigned int> partition) const
  {
    std::vector<unsigned int> fixedPartition = partition;
    unsigned int numMaxPartition = *std::max_element(begin(fixedPartition), end(fixedPartition));

    for (unsigned int e = 0; e < edgesConstrained.size(); e++)
    {
      if (!edgesConstrained[e])
        continue;

      if (mesh.Cell1DNumberNeighbourCell2D(e) < 2)
        continue;

      if (!mesh.Cell1DHasNeighbourCell2D(e, 0) ||
          !mesh.Cell1DHasNeighbourCell2D(e, 1))
        continue;

      const unsigned int neigh1 = mesh.Cell1DNeighbourCell2D(e,
                                                             0);
      const unsigned int neigh2 = mesh.Cell1DNeighbourCell2D(e,
                                                             1);

      if (fixedPartition.at(neigh1) !=
          fixedPartition.at(neigh2))
        continue;

      // fix partition
      fixedPartition[neigh1] = ++numMaxPartition;
    }

    return fixedPartition;
  }
  // ***************************************************************************
  std::vector<unsigned int> MetisUtilities::Mesh3DPartitionCheckConstraints(const IMeshDAO& mesh,
                                                                            const std::vector<bool>& facesConstrained,
                                                                            const std::vector<unsigned int> partition) const
  {
    std::vector<unsigned int> fixedPartition = partition;
    unsigned int numMaxPartition = *std::max_element(begin(fixedPartition), end(fixedPartition));

    for (unsigned int f = 0; f < facesConstrained.size(); f++)
    {
      if (!facesConstrained[f])
        continue;

      if (mesh.Cell2DNumberNeighbourCell3D(f) < 2)
        continue;

      const unsigned int neigh1 = mesh.Cell2DNeighbourCell3D(f,
                                                             0);
      const unsigned int neigh2 = mesh.Cell2DNeighbourCell3D(f,
                                                             1);

      if (fixedPartition.at(neigh1) !=
          fixedPartition.at(neigh2))
        continue;

      // fix partition
      fixedPartition[neigh1] = ++numMaxPartition;
    }

    return fixedPartition;
  }
  // ***************************************************************************
}
