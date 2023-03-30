#include "RefinementUtilities.hpp"


namespace Gedim
{
  // ***************************************************************************
  RefinementUtilities::RefinementUtilities(const GeometryUtilities& geometryUtilities,
                                           const MeshUtilities& meshUtilities) :
    geometryUtilities(geometryUtilities),
    meshUtilities(meshUtilities)
  {
  }
  RefinementUtilities::~RefinementUtilities()
  {
  }
  // ***************************************************************************
  RefinementUtilities::SplitCell1D_Result RefinementUtilities::SplitCell1D_MiddlePoint(const unsigned int& cell1DIndex,
                                                                                       IMeshDAO& mesh) const
  {
    SplitCell1D_Result result;

    const unsigned int cell1DOriginIndex = mesh.Cell1DOrigin(cell1DIndex);
    const unsigned int cell1DEndIndex = mesh.Cell1DEnd(cell1DIndex);

    result.NewCell0DIndex = mesh.Cell0DAppend(1);
    const unsigned int& newCell0DIndex = result.NewCell0DIndex;
    mesh.Cell0DInsertCoordinates(newCell0DIndex,
                                 0.5 * (mesh.Cell0DCoordinates(cell1DOriginIndex) +
                                        mesh.Cell0DCoordinates(cell1DEndIndex)));
    mesh.Cell0DSetMarker(newCell0DIndex,
                         mesh.Cell1DMarker(cell1DIndex));
    mesh.Cell0DSetState(newCell0DIndex,
                        true);

    // Split max edge into sub-edges
    Eigen::MatrixXi newCell1DsExtreme(2, 2);
    newCell1DsExtreme.col(0)<< cell1DOriginIndex, newCell0DIndex;
    newCell1DsExtreme.col(1)<< newCell0DIndex, cell1DEndIndex;
    result.NewCell1DsIndex = meshUtilities.SplitCell1D(cell1DIndex,
                                                       newCell1DsExtreme,
                                                       mesh);

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::SplitPolygon_Result RefinementUtilities::SplitPolygon_NoNewVertices(const unsigned int& cell2DIndex,
                                                                                           const unsigned int cell2DNumVertices,
                                                                                           const unsigned int& fromVertex,
                                                                                           const unsigned int& toVertex,
                                                                                           IMeshDAO& mesh) const
  {
    SplitPolygon_Result result;

    // Create new cell1D from vertex to vertex
    result.NewCell1DIndex = mesh.Cell1DAppend(1);
    const unsigned int& newCell1DIndex = result.NewCell1DIndex;
    mesh.Cell1DInsertExtremes(newCell1DIndex,
                              mesh.Cell2DVertex(cell2DIndex, toVertex),
                              mesh.Cell2DVertex(cell2DIndex, fromVertex));
    mesh.Cell1DSetMarker(newCell1DIndex, 0);
    mesh.Cell1DSetState(newCell1DIndex, true);
    mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex, 2);

    std::list<unsigned int> firstPolygonIndices;
    std::list<unsigned int> secondPolygonIndices;

    for (unsigned int v = 0; ((fromVertex + v) % cell2DNumVertices) != toVertex; v++)
      firstPolygonIndices.push_back((fromVertex + v) % cell2DNumVertices);
    for (unsigned int v = 0; ((toVertex + v) % cell2DNumVertices) != fromVertex; v++)
      secondPolygonIndices.push_back((toVertex + v) % cell2DNumVertices);

    // Split cell2D into sub-cells
    std::vector<Eigen::MatrixXi> subCells(2);

    subCells[0].resize(2, firstPolygonIndices.size() + 1);
    unsigned int v = 0;
    for (const unsigned int& index : firstPolygonIndices)
    {
      subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
      subCells[0](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
      v++;
    }
    subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, toVertex);
    subCells[0](1, v) = newCell1DIndex;

    subCells[1].resize(2, secondPolygonIndices.size() + 1);
    v = 0;
    for (const unsigned int& index : secondPolygonIndices)
    {
      subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
      subCells[1](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
      v++;
    }
    subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, fromVertex);
    subCells[1](1, v) = newCell1DIndex;

    result.NewCell2DsIndex = meshUtilities.SplitCell2D(cell2DIndex,
                                                       subCells,
                                                       mesh);

    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     0,
                                     result.NewCell2DsIndex[1]); // right
    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     1,
                                     result.NewCell2DsIndex[0]); // left

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::SplitPolygon_Result RefinementUtilities::SplitPolygon_NewVertexFrom(const unsigned int& cell2DIndex,
                                                                                           const unsigned int cell2DNumVertices,
                                                                                           const unsigned int& fromEdge,
                                                                                           const unsigned int& toVertex,
                                                                                           const unsigned int& fromNewCell0DIndex,
                                                                                           const std::vector<unsigned int>& fromSplitCell1DsIndex,
                                                                                           const bool& fromEdgeDirection,
                                                                                           IMeshDAO& mesh) const
  {
    SplitPolygon_Result result;

    // Create new cell1D from new vertex to vertex
    result.NewCell1DIndex = mesh.Cell1DAppend(1);
    const unsigned int& newCell1DIndex = result.NewCell1DIndex;
    mesh.Cell1DInsertExtremes(newCell1DIndex,
                              mesh.Cell2DVertex(cell2DIndex, toVertex),
                              fromNewCell0DIndex);
    mesh.Cell1DSetMarker(newCell1DIndex, 0);
    mesh.Cell1DSetState(newCell1DIndex, true);
    mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex, 2);

    std::list<unsigned int> firstPolygonIndices;
    std::list<unsigned int> secondPolygonIndices;

    for (unsigned int v = 0; ((fromEdge + 1 + v) % cell2DNumVertices) != toVertex; v++)
      firstPolygonIndices.push_back((fromEdge + 1 + v) % cell2DNumVertices);
    for (unsigned int v = 0; ((toVertex + v) % cell2DNumVertices) != fromEdge; v++)
      secondPolygonIndices.push_back((toVertex + v) % cell2DNumVertices);

    // Split cell2D into sub-cells
    std::vector<Eigen::MatrixXi> subCells(2);

    subCells[0].resize(2, firstPolygonIndices.size() + 2);
    unsigned int v = 0;
    subCells[0](0, v) = fromNewCell0DIndex;
    subCells[0](1, v) = fromEdgeDirection ? fromSplitCell1DsIndex[1] : fromSplitCell1DsIndex[0];
    v++;
    for (const unsigned int& index : firstPolygonIndices)
    {
      subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
      subCells[0](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
      v++;
    }
    subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, toVertex);
    subCells[0](1, v) = newCell1DIndex;

    subCells[1].resize(2, secondPolygonIndices.size() + 2);
    v = 0;
    for (const unsigned int& index : secondPolygonIndices)
    {
      subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
      subCells[1](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
      v++;
    }
    subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, fromEdge);
    subCells[1](1, v) = fromEdgeDirection ? fromSplitCell1DsIndex[0] : fromSplitCell1DsIndex[1];
    v++;
    subCells[1](0, v) = fromNewCell0DIndex;
    subCells[1](1, v) = newCell1DIndex;

    result.NewCell2DsIndex = meshUtilities.SplitCell2D(cell2DIndex,
                                                       subCells,
                                                       mesh);

    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     0,
                                     result.NewCell2DsIndex[1]); // right
    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     1,
                                     result.NewCell2DsIndex[0]); // left

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::SplitPolygon_Result RefinementUtilities::SplitPolygon_NewVertexTo(const unsigned int& cell2DIndex,
                                                                                         const unsigned int cell2DNumVertices,
                                                                                         const unsigned int& fromVertex,
                                                                                         const unsigned int& toEdge,
                                                                                         const unsigned int& toNewCell0DIndex,
                                                                                         const std::vector<unsigned int>& toSplitCell1DsIndex,
                                                                                         const bool& toEdgeDirection,
                                                                                         IMeshDAO& mesh) const
  {
    SplitPolygon_Result result;

    // Create new cell1D from vertex to new vertex
    result.NewCell1DIndex = mesh.Cell1DAppend(1);
    const unsigned int& newCell1DIndex = result.NewCell1DIndex;
    mesh.Cell1DInsertExtremes(newCell1DIndex,
                              toNewCell0DIndex,
                              mesh.Cell2DVertex(cell2DIndex, fromVertex));
    mesh.Cell1DSetMarker(newCell1DIndex, 0);
    mesh.Cell1DSetState(newCell1DIndex, true);
    mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex, 2);

    std::list<unsigned int> firstPolygonIndices;
    std::list<unsigned int> secondPolygonIndices;

    for (unsigned int v = 0; ((fromVertex + v) % cell2DNumVertices) != toEdge; v++)
      firstPolygonIndices.push_back((fromVertex + v) % cell2DNumVertices);
    for (unsigned int v = 0; ((toEdge + 1 + v) % cell2DNumVertices) != fromVertex; v++)
      secondPolygonIndices.push_back((toEdge + 1 + v) % cell2DNumVertices);

    // Split cell2D into sub-cells
    std::vector<Eigen::MatrixXi> subCells(2);

    subCells[0].resize(2, firstPolygonIndices.size() + 2);
    unsigned int v = 0;
    for (const unsigned int& index : firstPolygonIndices)
    {
      subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
      subCells[0](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
      v++;
    }
    subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, toEdge);
    subCells[0](1, v) = toEdgeDirection ? toSplitCell1DsIndex[0] : toSplitCell1DsIndex[1];
    v++;
    subCells[0](0, v) = toNewCell0DIndex;
    subCells[0](1, v) = newCell1DIndex;

    subCells[1].resize(2, secondPolygonIndices.size() + 2);
    v = 0;
    subCells[1](0, v) = toNewCell0DIndex;
    subCells[1](1, v) = toEdgeDirection ? toSplitCell1DsIndex[1] : toSplitCell1DsIndex[0];
    v++;
    for (const unsigned int& index : secondPolygonIndices)
    {
      subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
      subCells[1](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
      v++;
    }
    subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, fromVertex);
    subCells[1](1, v) = newCell1DIndex;

    result.NewCell2DsIndex = meshUtilities.SplitCell2D(cell2DIndex,
                                                       subCells,
                                                       mesh);

    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     0,
                                     result.NewCell2DsIndex[1]); // right
    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     1,
                                     result.NewCell2DsIndex[0]); // left

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::SplitPolygon_Result RefinementUtilities::SplitPolygon_NewVertices(const unsigned int& cell2DIndex,
                                                                                         const unsigned int cell2DNumVertices,
                                                                                         const unsigned int& fromEdge,
                                                                                         const unsigned int& toEdge,
                                                                                         const unsigned int& fromNewCell0DIndex,
                                                                                         const unsigned int& toNewCell0DIndex,
                                                                                         const std::vector<unsigned int>& fromSplitCell1DsIndex,
                                                                                         const std::vector<unsigned int>& toSplitCell1DsIndex,
                                                                                         const bool& fromEdgeDirection,
                                                                                         const bool& toEdgeDirection,
                                                                                         IMeshDAO& mesh) const
  {
    SplitPolygon_Result result;

    // Create new cell1D from new vertex to new vertex
    result.NewCell1DIndex = mesh.Cell1DAppend(1);
    const unsigned int& newCell1DIndex = result.NewCell1DIndex;
    mesh.Cell1DInsertExtremes(newCell1DIndex,
                              toNewCell0DIndex,
                              fromNewCell0DIndex);
    mesh.Cell1DSetMarker(newCell1DIndex, 0);
    mesh.Cell1DSetState(newCell1DIndex, true);
    mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex, 2);

    std::list<unsigned int> firstPolygonIndices;
    std::list<unsigned int> secondPolygonIndices;

    for (unsigned int v = 0; ((fromEdge + 1 + v) % cell2DNumVertices) != toEdge; v++)
      firstPolygonIndices.push_back((fromEdge + 1 + v) % cell2DNumVertices);
    for (unsigned int v = 0; ((toEdge + 1 + v) % cell2DNumVertices) != fromEdge; v++)
      secondPolygonIndices.push_back((toEdge + 1 + v) % cell2DNumVertices);

    // Split cell2D into sub-cells
    std::vector<Eigen::MatrixXi> subCells(2);

    subCells[0].resize(2, firstPolygonIndices.size() + 3);
    unsigned int v = 0;
    subCells[0](0, v) = fromNewCell0DIndex;
    subCells[0](1, v) = fromEdgeDirection ? fromSplitCell1DsIndex[1] : fromSplitCell1DsIndex[0];
    v++;
    for (const unsigned int& index : firstPolygonIndices)
    {
      subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
      subCells[0](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
      v++;
    }
    subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, toEdge);
    subCells[0](1, v) = toEdgeDirection ? toSplitCell1DsIndex[0] : toSplitCell1DsIndex[1];
    v++;
    subCells[0](0, v) = toNewCell0DIndex;
    subCells[0](1, v) = newCell1DIndex;

    subCells[1].resize(2, secondPolygonIndices.size() + 3);
    v = 0;
    subCells[1](0, v) = toNewCell0DIndex;
    subCells[1](1, v) = toEdgeDirection ? toSplitCell1DsIndex[1] : toSplitCell1DsIndex[0];
    v++;
    for (const unsigned int& index : secondPolygonIndices)
    {
      subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
      subCells[1](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
      v++;
    }
    subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, fromEdge);
    subCells[1](1, v) = fromEdgeDirection ? fromSplitCell1DsIndex[0] : fromSplitCell1DsIndex[1];
    v++;
    subCells[1](0, v) = fromNewCell0DIndex;
    subCells[1](1, v) = newCell1DIndex;

    result.NewCell2DsIndex = meshUtilities.SplitCell2D(cell2DIndex,
                                                       subCells,
                                                       mesh);

    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     0,
                                     result.NewCell2DsIndex[1]); // right
    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     1,
                                     result.NewCell2DsIndex[0]); // left

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::MaxEdgeDirection RefinementUtilities::ComputeTriangleMaxEdgeDirection(const Eigen::VectorXd& edgesLength) const
  {
    MaxEdgeDirection result;

    Eigen::VectorXd::Index maxEdgeLocalIndex;
    edgesLength.maxCoeff(&maxEdgeLocalIndex);
    result.MaxEdgeIndex = maxEdgeLocalIndex;
    result.OppositeVertexIndex = (maxEdgeLocalIndex + 2) % 3;

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::Cell2Ds_GeometricData RefinementUtilities::RefinePolygonCell_InitializeGeometricData(const IMeshDAO& mesh) const
  {
    Cell2Ds_GeometricData geometricData;

    std::vector<unsigned int> cell2DsIndex(mesh.Cell2DTotalNumber());
    std::iota(std::begin(cell2DsIndex), std::end(cell2DsIndex), 0);

    RefinePolygonCell_UpdateGeometricData(mesh,
                                          cell2DsIndex,
                                          geometricData);

    return geometricData;
  }
  // ***************************************************************************
  void RefinementUtilities::RefinePolygonCell_UpdateGeometricData(const IMeshDAO& mesh,
                                                                  const std::vector<unsigned int>& cell2DsIndex,
                                                                  Cell2Ds_GeometricData& geometricData) const
  {
    geometricData.Cell1Ds.Quality.resize(mesh.Cell1DTotalNumber());

    geometricData.Cell2Ds.Vertices.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.Area.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.Centroid.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.EdgesDirection.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.EdgesNormal.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.EdgesLength.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.Triangulations.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.Inertia.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.InRadius.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.Quality.resize(mesh.Cell2DTotalNumber());

    for (unsigned int c = 0; c < cell2DsIndex.size(); c++)
    {
      const unsigned int cell2DIndex = cell2DsIndex.at(c);

      Eigen::MatrixXd& convexCell2DVertices = geometricData.Cell2Ds.Vertices.at(cell2DIndex);
      double& convexCell2DArea = geometricData.Cell2Ds.Area.at(cell2DIndex);
      std::vector<Eigen::Matrix3d>& convexCell2DTriangulationPoints = geometricData.Cell2Ds.Triangulations.at(cell2DIndex);
      Eigen::Vector3d& convexCell2DCentroid = geometricData.Cell2Ds.Centroid.at(cell2DIndex);

      convexCell2DVertices = mesh.Cell2DVerticesCoordinates(cell2DIndex);

      // compute original cell2D triangulation
      const std::vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(convexCell2DVertices);
      const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(convexCell2DVertices,
                                                                                            convexCell2DUnalignedVerticesFilter);

      const std::vector<unsigned int> convexCell2DTriangulationFiltered = geometryUtilities.PolygonTriangulationByFirstVertex(convexCell2DUnalignedVertices);
      std::vector<unsigned int> convexCell2DTriangulation(convexCell2DTriangulationFiltered.size());
      for (unsigned int ocf = 0; ocf < convexCell2DTriangulationFiltered.size(); ocf++)
        convexCell2DTriangulation[ocf] = convexCell2DUnalignedVerticesFilter[convexCell2DTriangulationFiltered[ocf]];

      convexCell2DTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(convexCell2DVertices,
                                                                                     convexCell2DTriangulation);

      const unsigned int& numTriangles = convexCell2DTriangulationPoints.size();

      // compute original cell2D area and centroids
      Eigen::VectorXd convexCell2DTriangulationAreas(numTriangles);
      Eigen::MatrixXd convexCell2DTriangulationCentroids(3, numTriangles);
      for (unsigned int cct = 0; cct < numTriangles; cct++)
      {
        convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cct]);
        convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonBarycenter(convexCell2DTriangulationPoints[cct]);
      }

      convexCell2DArea = convexCell2DTriangulationAreas.sum();
      convexCell2DCentroid = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                               convexCell2DTriangulationAreas,
                                                               convexCell2DArea);

      geometricData.Cell2Ds.EdgesDirection[cell2DIndex].resize(mesh.Cell2DNumberEdges(cell2DIndex));
      for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(cell2DIndex); e++)
      {
        const unsigned int origin = mesh.Cell2DVertex(cell2DIndex, e);
        const unsigned int end = mesh.Cell2DVertex(cell2DIndex,
                                                   (e + 1) % mesh.Cell2DNumberEdges(cell2DIndex));

        geometricData.Cell2Ds.EdgesDirection[cell2DIndex][e] = mesh.Cell1DExists(origin,
                                                                                 end);
      }

      geometricData.Cell2Ds.EdgesLength[cell2DIndex] = geometryUtilities.PolygonEdgeLengths(convexCell2DVertices);

      geometricData.Cell2Ds.EdgesNormal[cell2DIndex] = geometryUtilities.PolygonEdgeNormals(convexCell2DVertices);
      geometricData.Cell2Ds.InRadius[cell2DIndex] = geometryUtilities.PolygonInRadius(convexCell2DVertices,
                                                                                      convexCell2DCentroid,
                                                                                      geometricData.Cell2Ds.EdgesNormal[cell2DIndex]);
      geometricData.Cell2Ds.Inertia[cell2DIndex] = geometryUtilities.PolygonInertia(convexCell2DCentroid,
                                                                                    convexCell2DTriangulationPoints);

      geometricData.Cell2Ds.Quality[cell2DIndex] = std::min(geometricData.Cell2Ds.EdgesLength[cell2DIndex].minCoeff(),
                                                            geometricData.Cell2Ds.InRadius[cell2DIndex]);

      for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(cell2DIndex); e++)
      {
        const unsigned int cell1DIndex = mesh.Cell2DEdge(cell2DIndex, e);
        geometricData.Cell1Ds.Quality[cell1DIndex] = 0.0;

        for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell2D(cell1DIndex); n++)
        {
          if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, n))
            continue;

          const unsigned int neighCell2DIndex = mesh.Cell1DNeighbourCell2D(cell1DIndex, n);
          if (geometricData.Cell1Ds.Quality[cell1DIndex] < geometricData.Cell2Ds.Quality[neighCell2DIndex])
            geometricData.Cell1Ds.Quality[cell1DIndex] = geometricData.Cell2Ds.Quality[neighCell2DIndex];
        }
      }
    }
  }
  // ***************************************************************************
  bool RefinementUtilities::RefinePolygonCell_IsCell1DToSplit(const unsigned int& cell1DIndex,
                                                              const unsigned int& cell2DIndex,
                                                              const unsigned int& cell2DNumVertices,
                                                              const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersection,
                                                              const Eigen::VectorXd& cell2DEdgesLength,
                                                              const double& cell1DsQualityWeight,
                                                              const double& cell1DQuality,
                                                              const IMeshDAO& mesh) const
  {
    // check if the edge intersection is inside edge
    const bool isInsideIntersection = (edgeIntersection.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);

    if (!isInsideIntersection)
      return false;

    // check if the edge quality is enough
    const bool isQualityEnough = geometryUtilities.IsValue1DGreater(0.5 * cell2DEdgesLength[edgeIntersection.Index],
        cell1DsQualityWeight * cell1DQuality);

    if (!isQualityEnough)
      return false;

    // check if edge is part of aligned edges
    if (mesh.Cell1DHasOriginalCell1D(cell1DIndex))
    {
      const unsigned int originalCell1DIndex = mesh.Cell1DOriginalCell1D(cell1DIndex);

      // check previous edge
      const unsigned int previousCell1DIndex = edgeIntersection.Index == 0 ? mesh.Cell2DEdge(cell2DIndex,
                                                                                             cell2DNumVertices - 1) :
                                                                             mesh.Cell2DEdge(cell2DIndex,
                                                                                             edgeIntersection.Index - 1);
      if (mesh.Cell1DHasOriginalCell1D(previousCell1DIndex) &&
          originalCell1DIndex == mesh.Cell1DOriginalCell1D(previousCell1DIndex))
        return false;

      // check next edge
      const unsigned int nextCell1DIndex = mesh.Cell2DEdge(cell2DIndex,
                                                           (edgeIntersection.Index + 1) % cell2DNumVertices);
      if (mesh.Cell1DHasOriginalCell1D(nextCell1DIndex) &&
          originalCell1DIndex == mesh.Cell1DOriginalCell1D(nextCell1DIndex))
        return false;
    }

    return true;
  }
  // ***************************************************************************
  RefinementUtilities::PolygonDirection RefinementUtilities::ComputePolygonMaxDiameterDirection(const Eigen::MatrixXd& vertices,
                                                                                                const Eigen::Vector3d& centroid) const
  {
    PolygonDirection result;

    const std::vector<unsigned int> unalignedPoints = geometryUtilities.UnalignedPoints(vertices);
    const Eigen::MatrixXd unalignedVertices = geometryUtilities.ExtractPoints(vertices,
                                                                              unalignedPoints);

    const Eigen::MatrixXd distances = geometryUtilities.PointsDistance(unalignedVertices);
    Eigen::MatrixXd::Index row, col;
    distances.maxCoeff(&row, &col); // get vertex indices of maximum distance in the polygon
    Gedim::Output::Assert(col != row);

    result.LineTangent = geometryUtilities.SegmentNormal(unalignedVertices.col(row),
                                                         unalignedVertices.col(col));
    result.LineOrigin = centroid;

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::PolygonDirection RefinementUtilities::ComputePolygonMaxInertiaDirection(const Eigen::Vector3d& centroid,
                                                                                               const Eigen::Matrix3d& inertia) const
  {
    PolygonDirection result;

    // get maximux direction of inertia tensor
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(inertia.block(0, 0, 2, 2));
    if (eigensolver.info() != Eigen::Success)
      throw std::runtime_error("Inertia not correct");

    result.LineTangent.setZero();
    result.LineTangent.segment(0, 2) = eigensolver.eigenvectors().col(1);
    result.LineOrigin = centroid;

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolygon_Result RefinementUtilities::RefineTriangleCell_ByEdge(const unsigned int& cell2DIndex,
                                                                                           const unsigned int& edgeIndex,
                                                                                           const unsigned int& oppositeVertexIndex,
                                                                                           const std::vector<bool>& cell2DEdgesDirection,
                                                                                           const double& cell2DArea,
                                                                                           const Eigen::VectorXd& cell2DEdgesLength,
                                                                                           IMeshDAO& mesh) const
  {
    RefinePolygon_Result result;

    if (mesh.Cell2DHasUpdatedCell2Ds(cell2DIndex))
      return result;

    // Check if the cell refinement creates null cell2Ds or cell1Ds
    if (geometryUtilities.IsValue2DZero(cell2DArea / 2.0) ||
        geometryUtilities.IsValue1DZero(cell2DEdgesLength[edgeIndex] / 2.0))
      return result;

    // Add new mesh vertex, middle point of edge
    const unsigned int cell1DIndex = mesh.Cell2DEdge(cell2DIndex, edgeIndex);
    const bool cell1DDirection = cell2DEdgesDirection[edgeIndex];

    const SplitCell1D_Result splitCell1D = SplitCell1D_MiddlePoint(cell1DIndex,
                                                                   mesh);

    const SplitPolygon_Result splitResult = SplitPolygon_NewVertexTo(cell2DIndex,
                                                                     3,
                                                                     oppositeVertexIndex,
                                                                     edgeIndex,
                                                                     splitCell1D.NewCell0DIndex,
                                                                     splitCell1D.NewCell1DsIndex,
                                                                     cell1DDirection,
                                                                     mesh);

    result.NewCell0DsIndex = { splitCell1D.NewCell0DIndex };
    result.NewCell1DsIndex.resize(2);
    result.NewCell1DsIndex[0].Type = RefinePolygon_Result::RefinedCell1D::Types::Updated;
    result.NewCell1DsIndex[0].NewCell1DsIndex = splitCell1D.NewCell1DsIndex;
    result.NewCell1DsIndex[0].OriginalCell1DIndex = cell1DIndex;
    result.NewCell1DsIndex[0].OriginalCell2DEdgeIndex = edgeIndex;
    result.NewCell1DsIndex[0].NewCell0DIndex = splitCell1D.NewCell0DIndex;
    result.NewCell1DsIndex[1].Type = RefinePolygon_Result::RefinedCell1D::Types::New;
    result.NewCell1DsIndex[1].NewCell1DsIndex = { splitResult.NewCell1DIndex };

    result.NewCell2DsIndex = { splitResult.NewCell2DsIndex };

    return result;
  }
  // ***************************************************************************
  void RefinementUtilities::RefineTriangleCell_UpdateNeighbours(const unsigned int& cell2DIndex,
                                                                const unsigned int& cell1DIndex,
                                                                const unsigned int& newCell0DIndex,
                                                                const std::vector<unsigned int>& splitCell1DsIndex,
                                                                const bool& cell2DEdgeDirection,
                                                                IMeshDAO& mesh) const
  {
    // update neighbour cells
    for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell2D(cell1DIndex); n++)
    {
      if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, n))
        continue;

      const unsigned int neighCell2DIndex = mesh.Cell1DNeighbourCell2D(cell1DIndex, n);
      if (neighCell2DIndex == cell2DIndex)
        continue;

      if (mesh.Cell2DHasUpdatedCell2Ds(neighCell2DIndex))
        continue;

      const unsigned int neighEdgeIndex = mesh.Cell2DFindEdge(neighCell2DIndex,
                                                              cell1DIndex);
      const unsigned int neighOppositeVertexIndex = (neighEdgeIndex + 2) % 3;

      SplitPolygon_NewVertexTo(neighCell2DIndex,
                               3,
                               neighOppositeVertexIndex,
                               neighEdgeIndex,
                               newCell0DIndex,
                               splitCell1DsIndex,
                               !cell2DEdgeDirection,
                               mesh);
    }
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolygon_Result RefinementUtilities::RefinePolygonCell_ByDirection(const unsigned int& cell2DIndex,
                                                                                               const Eigen::MatrixXd& cell2DVertices,
                                                                                               const Eigen::Vector3d& lineTangent,
                                                                                               const Eigen::Vector3d& lineOrigin,
                                                                                               const std::vector<double>& cell1DsQuality,
                                                                                               const double& cell1DsQualityWeight,
                                                                                               const Eigen::VectorXd& cell2DEdgesLength,
                                                                                               const std::vector<bool>& cell2DEdgesDirection,
                                                                                               IMeshDAO& mesh) const
  {
    RefinePolygon_Result result;

    if (mesh.Cell2DHasUpdatedCell2Ds(cell2DIndex))
      return result;

    const unsigned int cell2DNumVertices = cell2DVertices.cols();

    GeometryUtilities::LinePolygonPositionResult linePolygonPosition = geometryUtilities.LinePolygonPosition(lineTangent,
                                                                                                             lineOrigin,
                                                                                                             cell2DVertices);
    Output::Assert(linePolygonPosition.Type != GeometryUtilities::LinePolygonPositionResult::Types::Unknown);

    if (linePolygonPosition.Type != GeometryUtilities::LinePolygonPositionResult::Types::Intersecting)
      return result;

    if (linePolygonPosition.EdgeIntersections.size() < 2)
      throw std::runtime_error("Not enough intersections found");

    const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersectionOne = linePolygonPosition.EdgeIntersections[0];
    Output::Assert(edgeIntersectionOne.Type != GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::Unknown);
    if (edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::Parallel)
      return result;

    const unsigned int cell1DIndexOne = mesh.Cell2DEdge(cell2DIndex,
                                                        edgeIntersectionOne.Index);

    GeometryUtilities::LinePolygonPositionResult::EdgeIntersection* findEdgeIntersectionTwo = nullptr;

    for (unsigned int i = 1; i < linePolygonPosition.EdgeIntersections.size(); i++)
    {
      const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersectionTwo = linePolygonPosition.EdgeIntersections[i];

      switch (edgeIntersectionTwo.Type)
      {
        case GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::Parallel:
          return result;
        case GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge:
          findEdgeIntersectionTwo = &linePolygonPosition.EdgeIntersections[i];
          break;
        case GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeOrigin:
        {
          const unsigned int cell1DIndexTwo = mesh.Cell2DEdge(cell2DIndex, edgeIntersectionTwo.Index);
          const unsigned int cell1DOriginIndex = cell2DEdgesDirection[edgeIntersectionTwo.Index] ? mesh.Cell1DOrigin(cell1DIndexTwo) :
                                                                                                   mesh.Cell1DEnd(cell1DIndexTwo);
          if (mesh.Cell1DOrigin(cell1DIndexOne) != cell1DOriginIndex &&
              mesh.Cell1DEnd(cell1DIndexOne) != cell1DOriginIndex)
            findEdgeIntersectionTwo = &linePolygonPosition.EdgeIntersections[i];
        }
          break;
        case GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeEnd:
        {
          const unsigned int cell1DIndexTwo = mesh.Cell2DEdge(cell2DIndex, edgeIntersectionTwo.Index);
          const unsigned int cell1DEndIndex = cell2DEdgesDirection[edgeIntersectionTwo.Index] ? mesh.Cell1DOrigin(cell1DIndexTwo) :
                                                                                                mesh.Cell1DEnd(cell1DIndexTwo);
          if (mesh.Cell1DOrigin(cell1DIndexOne) != cell1DEndIndex &&
              mesh.Cell1DEnd(cell1DIndexOne) != cell1DEndIndex)
            findEdgeIntersectionTwo = &linePolygonPosition.EdgeIntersections[i];
        }
          break;
        default:
          throw std::runtime_error("Unmanaged edgeIntersectionTwo.Type");
      }

      if (findEdgeIntersectionTwo != nullptr)
        break;
    }

    if (findEdgeIntersectionTwo == nullptr)
      throw std::runtime_error("Second edge intersection not found");

    const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersectionTwo = *findEdgeIntersectionTwo;

    const unsigned int cell1DIndexTwo = mesh.Cell2DEdge(cell2DIndex,
                                                        edgeIntersectionTwo.Index);

    const bool createNewVertexOne = RefinePolygonCell_IsCell1DToSplit(cell1DIndexOne,
                                                                      cell2DIndex,
                                                                      cell2DNumVertices,
                                                                      edgeIntersectionOne,
                                                                      cell2DEdgesLength,
                                                                      cell1DsQualityWeight,
                                                                      cell1DsQuality[cell1DIndexOne],
                                                                      mesh);
    const bool createNewVertexTwo = RefinePolygonCell_IsCell1DToSplit(cell1DIndexTwo,
                                                                      cell2DIndex,
                                                                      cell2DNumVertices,
                                                                      edgeIntersectionTwo,
                                                                      cell2DEdgesLength,
                                                                      cell1DsQualityWeight,
                                                                      cell1DsQuality[cell1DIndexTwo],
                                                                      mesh);

    if (!createNewVertexOne && !createNewVertexTwo)
    {
      // no new vertices
      const unsigned int fromVertex = edgeIntersectionOne.CurvilinearCoordinate <= 0.5 ?
                                        edgeIntersectionOne.Index :
                                        (edgeIntersectionOne.Index + 1) % cell2DNumVertices;
      const unsigned int toVertex = edgeIntersectionTwo.CurvilinearCoordinate <= 0.5 ?
                                      edgeIntersectionTwo.Index :
                                      (edgeIntersectionTwo.Index + 1) % cell2DNumVertices;

      if ((fromVertex + 1) % cell2DNumVertices == toVertex ||
          (toVertex + 1) % cell2DNumVertices == fromVertex)
        return result;

      if (geometryUtilities.PointIsAligned(cell2DVertices.col(fromVertex),
                                           cell2DVertices.col(toVertex),
                                           cell2DVertices.col((fromVertex + 1) % cell2DNumVertices)) ||
          geometryUtilities.PointIsAligned(cell2DVertices.col(fromVertex),
                                           cell2DVertices.col(toVertex),
                                           cell2DVertices.col((toVertex + 1) % cell2DNumVertices)))
        return result;

      const SplitPolygon_Result splitResult = SplitPolygon_NoNewVertices(cell2DIndex,
                                                                         cell2DNumVertices,
                                                                         fromVertex,
                                                                         toVertex,
                                                                         mesh);

      result.NewCell1DsIndex.resize(1);
      result.NewCell1DsIndex[0].Type = RefinePolygon_Result::RefinedCell1D::Types::New;
      result.NewCell1DsIndex[0].NewCell1DsIndex = { splitResult.NewCell1DIndex };

      result.NewCell2DsIndex = splitResult.NewCell2DsIndex;
    }
    else if (createNewVertexOne && !createNewVertexTwo)
    {
      // new vertex one
      Gedim::Output::Assert(edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);

      const unsigned int toVertex = edgeIntersectionTwo.CurvilinearCoordinate <= 0.5 ?
                                      edgeIntersectionTwo.Index :
                                      (edgeIntersectionTwo.Index + 1) % cell2DNumVertices;

      const Eigen::Vector3d middleCoordinate = 0.5 * (cell2DVertices.col(edgeIntersectionOne.Index) +
                                                      cell2DVertices.col((edgeIntersectionOne.Index + 1) % cell2DNumVertices));

      if (geometryUtilities.PointIsAligned(cell2DVertices.col(edgeIntersectionOne.Index),
                                           cell2DVertices.col(toVertex),
                                           middleCoordinate) ||
          geometryUtilities.PointIsAligned(cell2DVertices.col((edgeIntersectionOne.Index + 1) % cell2DNumVertices),
                                           cell2DVertices.col(toVertex),
                                           middleCoordinate))
        return result;

      const SplitCell1D_Result splitCell1DOne = SplitCell1D_MiddlePoint(cell1DIndexOne,
                                                                        mesh);

      const SplitPolygon_Result splitResult = SplitPolygon_NewVertexFrom(cell2DIndex,
                                                                         cell2DNumVertices,
                                                                         edgeIntersectionOne.Index,
                                                                         toVertex,
                                                                         splitCell1DOne.NewCell0DIndex,
                                                                         splitCell1DOne.NewCell1DsIndex,
                                                                         cell2DEdgesDirection.at(edgeIntersectionOne.Index),
                                                                         mesh);

      result.NewCell0DsIndex = { splitCell1DOne.NewCell0DIndex };

      result.NewCell1DsIndex.resize(2);
      result.NewCell1DsIndex[0].Type = RefinePolygon_Result::RefinedCell1D::Types::Updated;
      result.NewCell1DsIndex[0].NewCell1DsIndex = splitCell1DOne.NewCell1DsIndex;
      result.NewCell1DsIndex[0].OriginalCell1DIndex = cell1DIndexOne;
      result.NewCell1DsIndex[0].OriginalCell2DEdgeIndex = edgeIntersectionOne.Index;
      result.NewCell1DsIndex[0].NewCell0DIndex = splitCell1DOne.NewCell0DIndex;
      result.NewCell1DsIndex[1].Type = RefinePolygon_Result::RefinedCell1D::Types::New;
      result.NewCell1DsIndex[1].NewCell1DsIndex = { splitResult.NewCell1DIndex };

      result.NewCell2DsIndex = splitResult.NewCell2DsIndex;
    }
    else if (!createNewVertexOne && createNewVertexTwo)
    {
      // new vertex two
      Gedim::Output::Assert(edgeIntersectionTwo.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);

      const unsigned int fromVertex = edgeIntersectionOne.CurvilinearCoordinate <= 0.5 ?
                                        edgeIntersectionOne.Index :
                                        (edgeIntersectionOne.Index + 1) % cell2DNumVertices;

      const Eigen::Vector3d middleCoordinate = 0.5 * (cell2DVertices.col(edgeIntersectionTwo.Index) +
                                                      cell2DVertices.col((edgeIntersectionTwo.Index + 1) % cell2DNumVertices));

      if (geometryUtilities.PointIsAligned(cell2DVertices.col(fromVertex),
                                           cell2DVertices.col(edgeIntersectionTwo.Index),
                                           middleCoordinate) ||
          geometryUtilities.PointIsAligned(cell2DVertices.col(fromVertex),
                                           cell2DVertices.col((edgeIntersectionTwo.Index + 1) % cell2DNumVertices),
                                           middleCoordinate))
        return result;

      const SplitCell1D_Result splitCell1DTwo = SplitCell1D_MiddlePoint(cell1DIndexTwo,
                                                                        mesh);

      const SplitPolygon_Result splitResult = SplitPolygon_NewVertexTo(cell2DIndex,
                                                                       cell2DNumVertices,
                                                                       fromVertex,
                                                                       edgeIntersectionTwo.Index,
                                                                       splitCell1DTwo.NewCell0DIndex,
                                                                       splitCell1DTwo.NewCell1DsIndex,
                                                                       cell2DEdgesDirection.at(edgeIntersectionTwo.Index),
                                                                       mesh);

      result.NewCell0DsIndex = { splitCell1DTwo.NewCell0DIndex };

      result.NewCell1DsIndex.resize(2);
      result.NewCell1DsIndex[0].Type = RefinePolygon_Result::RefinedCell1D::Types::Updated;
      result.NewCell1DsIndex[0].NewCell1DsIndex = splitCell1DTwo.NewCell1DsIndex;
      result.NewCell1DsIndex[0].OriginalCell1DIndex = cell1DIndexTwo;
      result.NewCell1DsIndex[0].OriginalCell2DEdgeIndex = edgeIntersectionTwo.Index;
      result.NewCell1DsIndex[0].NewCell0DIndex = splitCell1DTwo.NewCell0DIndex;
      result.NewCell1DsIndex[1].Type = RefinePolygon_Result::RefinedCell1D::Types::New;
      result.NewCell1DsIndex[1].NewCell1DsIndex = { splitResult.NewCell1DIndex };

      result.NewCell2DsIndex = splitResult.NewCell2DsIndex;
    }
    else
    {
      // two new vertices
      Gedim::Output::Assert(edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);
      Gedim::Output::Assert(edgeIntersectionTwo.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);

      const SplitCell1D_Result splitCell1DOne = SplitCell1D_MiddlePoint(cell1DIndexOne,
                                                                        mesh);
      const SplitCell1D_Result splitCell1DTwo = SplitCell1D_MiddlePoint(cell1DIndexTwo,
                                                                        mesh);


      const SplitPolygon_Result splitResult = SplitPolygon_NewVertices(cell2DIndex,
                                                                       cell2DNumVertices,
                                                                       edgeIntersectionOne.Index,
                                                                       edgeIntersectionTwo.Index,
                                                                       splitCell1DOne.NewCell0DIndex,
                                                                       splitCell1DTwo.NewCell0DIndex,
                                                                       splitCell1DOne.NewCell1DsIndex,
                                                                       splitCell1DTwo.NewCell1DsIndex,
                                                                       cell2DEdgesDirection.at(edgeIntersectionOne.Index),
                                                                       cell2DEdgesDirection.at(edgeIntersectionTwo.Index),
                                                                       mesh);


      result.NewCell0DsIndex = { splitCell1DOne.NewCell0DIndex,
                                 splitCell1DTwo.NewCell0DIndex };

      result.NewCell1DsIndex.resize(3);
      result.NewCell1DsIndex[0].Type = RefinePolygon_Result::RefinedCell1D::Types::Updated;
      result.NewCell1DsIndex[0].NewCell1DsIndex = splitCell1DOne.NewCell1DsIndex;
      result.NewCell1DsIndex[0].OriginalCell1DIndex = cell1DIndexOne;
      result.NewCell1DsIndex[0].OriginalCell2DEdgeIndex = edgeIntersectionOne.Index;
      result.NewCell1DsIndex[0].NewCell0DIndex = splitCell1DOne.NewCell0DIndex;
      result.NewCell1DsIndex[1].Type = RefinePolygon_Result::RefinedCell1D::Types::Updated;
      result.NewCell1DsIndex[1].NewCell1DsIndex = splitCell1DTwo.NewCell1DsIndex;
      result.NewCell1DsIndex[1].OriginalCell1DIndex = cell1DIndexTwo;
      result.NewCell1DsIndex[1].OriginalCell2DEdgeIndex = edgeIntersectionTwo.Index;
      result.NewCell1DsIndex[1].NewCell0DIndex = splitCell1DTwo.NewCell0DIndex;
      result.NewCell1DsIndex[2].Type = RefinePolygon_Result::RefinedCell1D::Types::New;
      result.NewCell1DsIndex[2].NewCell1DsIndex = { splitResult.NewCell1DIndex };

      result.NewCell2DsIndex = splitResult.NewCell2DsIndex;
    }

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolygon_UpdateNeighbour_Result RefinementUtilities::RefinePolygonCell_UpdateNeighbours(const unsigned int& cell2DIndex,
                                                                                                                    const unsigned int& cell1DIndex,
                                                                                                                    const unsigned int& newCell0DIndex,
                                                                                                                    const std::vector<unsigned int>& splitCell1DsIndex,
                                                                                                                    const std::vector<std::vector<bool>>& cell2DsEdgesDirection,
                                                                                                                    IMeshDAO& mesh) const
  {
    RefinePolygon_UpdateNeighbour_Result result;

    std::list<RefinePolygon_UpdateNeighbour_Result::UpdatedCell2D> newCell2DsIndex;

    // update neighbour cells
    for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell2D(cell1DIndex); n++)
    {
      if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, n))
        continue;

      const unsigned int neighCell2DIndex = mesh.Cell1DNeighbourCell2D(cell1DIndex, n);
      if (neighCell2DIndex == cell2DIndex)
        continue;

      if (mesh.Cell2DHasUpdatedCell2Ds(neighCell2DIndex))
        continue;

      const unsigned int neighCell2DNumVertices = mesh.Cell2DNumberVertices(neighCell2DIndex);
      const unsigned int neighEdgeIndex = mesh.Cell2DFindEdge(neighCell2DIndex,
                                                              cell1DIndex);


      const bool& cell2DEdgeDirection = cell2DsEdgesDirection[neighCell2DIndex][neighEdgeIndex];

      std::vector<Eigen::MatrixXi> subCells(1);

      subCells[0].resize(2, neighCell2DNumVertices + 1);
      unsigned int ne = 0;
      for (unsigned int e = 0; e < neighEdgeIndex; e++)
      {
        subCells[0](0, ne) = mesh.Cell2DVertex(neighCell2DIndex, e);
        subCells[0](1, ne) = mesh.Cell2DEdge(neighCell2DIndex, e);
        ne++;
      }
      subCells[0](0, ne) = mesh.Cell2DVertex(neighCell2DIndex, neighEdgeIndex);
      subCells[0](1, ne) = cell2DEdgeDirection ? splitCell1DsIndex[0] : splitCell1DsIndex[1];
      ne++;
      subCells[0](0, ne) = newCell0DIndex;
      subCells[0](1, ne) = cell2DEdgeDirection ? splitCell1DsIndex[1] : splitCell1DsIndex[0];
      ne++;
      for (unsigned int e = neighEdgeIndex + 1; e < neighCell2DNumVertices; e++)
      {
        subCells[0](0, ne) = mesh.Cell2DVertex(neighCell2DIndex, e);
        subCells[0](1, ne) = mesh.Cell2DEdge(neighCell2DIndex, e);
        ne++;
      }

      const std::vector<unsigned int> newCell2Ds = meshUtilities.SplitCell2D(neighCell2DIndex,
                                                                             subCells,
                                                                             mesh);

      for (const unsigned int& newCell2D : newCell2Ds)
      {
        RefinePolygon_UpdateNeighbour_Result::UpdatedCell2D updatedCell2D;
        updatedCell2D.OriginalCell2DIndex = neighCell2DIndex;
        updatedCell2D.NewCell2DIndex = newCell2D;
        newCell2DsIndex.push_back(updatedCell2D);
      }
    }

    result.UpdatedCell2Ds = std::vector<RefinePolygon_UpdateNeighbour_Result::UpdatedCell2D>(newCell2DsIndex.begin(),
                                                                                             newCell2DsIndex.end());

    return result;
  }
  // ***************************************************************************
}
