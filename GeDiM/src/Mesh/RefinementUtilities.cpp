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
  RefinementUtilities::MaxEdgeDirection RefinementUtilities::ComputeTriangleMaxEdgeDirection(const Eigen::VectorXd& edgesLength)
  {
    MaxEdgeDirection result;

    Eigen::VectorXd::Index maxEdgeLocalIndex;
    edgesLength.maxCoeff(&maxEdgeLocalIndex);
    result.MaxEdgeIndex = maxEdgeLocalIndex;
    result.OppositeVertexIndex = (maxEdgeLocalIndex + 2) % 3;

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolygon_Result RefinementUtilities::RefineTriangleCellByEdge(const unsigned int& cell2DIndex,
                                                                                          const unsigned int& edgeIndex,
                                                                                          const unsigned int& oppositeVertexIndex,
                                                                                          const std::vector<bool>& cell2DEdgesDirection,
                                                                                          IMeshDAO& mesh) const
  {
    RefinePolygon_Result result;

    if (mesh.Cell2DHasUpdatedCell2Ds(cell2DIndex))
      return result;

    // Add new mesh vertex, middle point of edge
    const unsigned int cell1DIndex = mesh.Cell2DEdge(cell2DIndex, edgeIndex);
    const bool cell1DDirection = cell2DEdgesDirection[edgeIndex];
    const unsigned int cell1DOriginIndex = mesh.Cell1DOrigin(cell1DIndex);
    const unsigned int cell1DEndIndex = mesh.Cell1DEnd(cell1DIndex);

    const unsigned int newCell0DIndex = mesh.Cell0DAppend(1);
    mesh.Cell0DInsertCoordinates(newCell0DIndex,
                                 0.5 * (mesh.Cell0DCoordinates(cell1DOriginIndex) + mesh.Cell0DCoordinates(cell1DEndIndex)));
    mesh.Cell0DSetMarker(newCell0DIndex, mesh.Cell1DMarker(cell1DIndex));
    mesh.Cell0DSetState(newCell0DIndex, true);

    // Split max edge into sub-edges
    Eigen::MatrixXi newCell1DsExtreme(2, 2);
    newCell1DsExtreme.col(0)<< cell1DOriginIndex, newCell0DIndex;
    newCell1DsExtreme.col(1)<< newCell0DIndex, cell1DEndIndex;
    const std::vector<unsigned int> splitCell1DsIndex = meshUtilities.SplitCell1D(cell1DIndex,
                                                                                  newCell1DsExtreme,
                                                                                  mesh);
    SplitPolygon_NewVertexTo(cell2DIndex,
                             3,
                             oppositeVertexIndex,
                             edgeIndex,
                             newCell0DIndex,
                             splitCell1DsIndex,
                             cell1DDirection,
                             mesh);

    // update neighbour cell
    const unsigned int neighIndex = cell1DDirection ? 0 : 1;
    if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, neighIndex))
      return result;

    const unsigned int neighCell2DIndex = mesh.Cell1DNeighbourCell2D(cell1DIndex, neighIndex);
    const unsigned int neighEdgeIndex = mesh.Cell2DFindEdge(neighCell2DIndex,
                                                            cell1DIndex);
    const unsigned int neighOppositeVertexIndex = (neighEdgeIndex + 2) % 3;

    SplitPolygon_NewVertexTo(neighCell2DIndex,
                             3,
                             neighOppositeVertexIndex,
                             neighEdgeIndex,
                             newCell0DIndex,
                             splitCell1DsIndex,
                             !cell1DDirection,
                             mesh);

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolygon_Result RefinementUtilities::RefinePolygonalCellByDirection(const unsigned int& cell2DIndex,
                                                                                                const Eigen::MatrixXd& cell2DVertices,
                                                                                                const Eigen::Vector3d& lineTangent,
                                                                                                const Eigen::Vector3d& lineOrigin,
                                                                                                const std::vector<double>& cell1DsQualityParameter,
                                                                                                const Eigen::VectorXd& cell2DEdgesLength,
                                                                                                const std::vector<bool>& cell2DEdgesDirection,
                                                                                                IMeshDAO& mesh) const
  {
    if (mesh.Cell2DHasUpdatedCell2Ds(cell2DIndex))
      return;

    const unsigned int cell2DNumVertices = cell2DVertices.cols();

    GeometryUtilities::LinePolygonPositionResult result = geometryUtilities.LinePolygonPosition(lineTangent,
                                                                                                lineOrigin,
                                                                                                cell2DVertices);
    Output::Assert(result.Type != GeometryUtilities::LinePolygonPositionResult::Types::Unknown);

    if (result.Type != GeometryUtilities::LinePolygonPositionResult::Types::Intersecting)
      return;

    if (result.EdgeIntersections.size() < 2)
      throw std::runtime_error("Not enough intersections found");

    const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersectionOne = result.EdgeIntersections[0];
    Output::Assert(edgeIntersectionOne.Type != GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::Unknown);
    if (edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::Parallel)
      return;

    const unsigned int cell1DIndexOne = mesh.Cell2DEdge(cell2DIndex, edgeIntersectionOne.Index);

    GeometryUtilities::LinePolygonPositionResult::EdgeIntersection* findEdgeIntersectionTwo = nullptr;

    for (unsigned int i = 1; i < result.EdgeIntersections.size(); i++)
    {
      const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersectionTwo = result.EdgeIntersections[i];

      switch (edgeIntersectionTwo.Type)
      {
        case GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::Parallel:
          return;
        case GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge:
          findEdgeIntersectionTwo = &result.EdgeIntersections[i];
          break;
        case GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeOrigin:
        {
          const unsigned int cell1DIndexTwo = mesh.Cell2DEdge(cell2DIndex, edgeIntersectionTwo.Index);
          if (mesh.Cell1DOrigin(cell1DIndexOne) != mesh.Cell1DOrigin(cell1DIndexTwo) &&
              mesh.Cell1DEnd(cell1DIndexOne) != mesh.Cell1DOrigin(cell1DIndexTwo))
            findEdgeIntersectionTwo = &result.EdgeIntersections[i];
        }
          break;
        case GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeEnd:
        {
          const unsigned int cell1DIndexTwo = mesh.Cell2DEdge(cell2DIndex, edgeIntersectionTwo.Index);
          if (mesh.Cell1DOrigin(cell1DIndexOne) != mesh.Cell1DEnd(cell1DIndexTwo) &&
              mesh.Cell1DEnd(cell1DIndexOne) != mesh.Cell1DEnd(cell1DIndexTwo))
            findEdgeIntersectionTwo = &result.EdgeIntersections[i];
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

    const unsigned int cell1DIndexTwo = mesh.Cell2DEdge(cell2DIndex, edgeIntersectionTwo.Index);

    const bool createNewVertexOne =
        (edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge) &&
        geometryUtilities.IsValue1DGreater(cell2DEdgesLength[edgeIntersectionOne.Index], cell1DsQualityParameter[cell1DIndexOne]);
    const bool createNewVertexTwo =
        (edgeIntersectionTwo.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge) &&
        geometryUtilities.IsValue1DGreater(cell2DEdgesLength[edgeIntersectionTwo.Index], cell1DsQualityParameter[cell1DIndexTwo]);

    unsigned int newCell0DIndexOne = 0;
    unsigned int newCell0DIndexTwo = 0;
    std::vector<unsigned int> splitCell1DsIndexOne;
    std::vector<unsigned int> splitCell1DsIndexTwo;

    if (createNewVertexOne)
    {
      const unsigned int cell1DOriginIndex = mesh.Cell1DOrigin(cell1DIndexOne);
      const unsigned int cell1DEndIndex = mesh.Cell1DEnd(cell1DIndexOne);

      // Add new mesh vertex, middle point of edge one
      newCell0DIndexOne = mesh.Cell0DAppend(1);
      mesh.Cell0DInsertCoordinates(newCell0DIndexOne,
                                   0.5 * (mesh.Cell0DCoordinates(cell1DOriginIndex) +
                                          mesh.Cell0DCoordinates(cell1DEndIndex)));
      mesh.Cell0DSetMarker(newCell0DIndexOne, mesh.Cell1DMarker(cell1DIndexOne));
      mesh.Cell0DSetState(newCell0DIndexOne, true);

      // Split max edge into sub-edges
      Eigen::MatrixXi newCell1DsExtreme(2, 2);
      newCell1DsExtreme.col(0)<< cell1DOriginIndex, newCell0DIndexOne;
      newCell1DsExtreme.col(1)<< newCell0DIndexOne, cell1DEndIndex;
      splitCell1DsIndexOne = meshUtilities.SplitCell1D(cell1DIndexOne,
                                                       newCell1DsExtreme,
                                                       mesh);
    }

    if (createNewVertexTwo)
    {
      const unsigned int cell1DOriginIndex = mesh.Cell1DOrigin(cell1DIndexTwo);
      const unsigned int cell1DEndIndex = mesh.Cell1DEnd(cell1DIndexTwo);

      // Add new mesh vertex, middle point of edge two
      newCell0DIndexTwo = mesh.Cell0DAppend(1);
      mesh.Cell0DInsertCoordinates(newCell0DIndexTwo,
                                   0.5 * (mesh.Cell0DCoordinates(cell1DOriginIndex) +
                                          mesh.Cell0DCoordinates(cell1DEndIndex)));
      mesh.Cell0DSetMarker(newCell0DIndexTwo, mesh.Cell1DMarker(cell1DIndexTwo));
      mesh.Cell0DSetState(newCell0DIndexTwo, true);

      // Split max edge into sub-edges
      Eigen::MatrixXi newCell1DsExtreme(2, 2);
      newCell1DsExtreme.col(0)<< cell1DOriginIndex, newCell0DIndexTwo;
      newCell1DsExtreme.col(1)<< newCell0DIndexTwo, cell1DEndIndex;
      splitCell1DsIndexTwo = meshUtilities.SplitCell1D(cell1DIndexTwo,
                                                       newCell1DsExtreme,
                                                       mesh);
    }


    if (!createNewVertexOne && !createNewVertexTwo)
    {
      // no new vertices
      Gedim::Output::Assert(edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeOrigin ||
                            edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeEnd);
      Gedim::Output::Assert(edgeIntersectionTwo.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeOrigin ||
                            edgeIntersectionTwo.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeEnd);

      const unsigned int fromVertex = edgeIntersectionOne.Type ==
                                      GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeOrigin ?
                                        edgeIntersectionOne.Index :
                                        (edgeIntersectionOne.Index + 1) % cell2DNumVertices;
      const unsigned int toVertex = edgeIntersectionTwo.Type ==
                                    GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeOrigin ?
                                      edgeIntersectionTwo.Index :
                                      (edgeIntersectionTwo.Index + 1) % cell2DNumVertices;

      SplitPolygon_NoNewVertices(cell2DIndex,
                                 cell2DNumVertices,
                                 fromVertex,
                                 toVertex,
                                 mesh);
    }
    else if (createNewVertexOne && !createNewVertexTwo)
    {
      // new vertex one
      Gedim::Output::Assert(edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);
      Gedim::Output::Assert(edgeIntersectionTwo.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeOrigin ||
                            edgeIntersectionTwo.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeEnd);

      const unsigned int toVertex = edgeIntersectionTwo.Type ==
                                    GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeOrigin ?
                                      edgeIntersectionTwo.Index :
                                      (edgeIntersectionTwo.Index + 1) % cell2DNumVertices;

      SplitPolygon_NewVertexFrom(cell2DIndex,
                                 cell2DNumVertices,
                                 edgeIntersectionOne.Index,
                                 toVertex,
                                 newCell0DIndexOne,
                                 splitCell1DsIndexOne,
                                 cell2DEdgesDirection.at(edgeIntersectionOne.Index),
                                 mesh);
    }
    else if (!createNewVertexOne && createNewVertexTwo)
    {
      // new vertex two
      Gedim::Output::Assert(edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeOrigin ||
                            edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeEnd);
      Gedim::Output::Assert(edgeIntersectionTwo.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);

      const unsigned int fromVertex = edgeIntersectionOne.Type ==
                                      GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::OnEdgeOrigin ?
                                        edgeIntersectionOne.Index :
                                        (edgeIntersectionOne.Index + 1) % cell2DNumVertices;

      SplitPolygon_NewVertexTo(cell2DIndex,
                               cell2DNumVertices,
                               fromVertex,
                               edgeIntersectionTwo.Index,
                               newCell0DIndexTwo,
                               splitCell1DsIndexTwo,
                               cell2DEdgesDirection.at(edgeIntersectionTwo.Index),
                               mesh);
    }
    else
    {
      // two new vertices
      Gedim::Output::Assert(edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);
      Gedim::Output::Assert(edgeIntersectionTwo.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);

      SplitPolygon_NewVertices(cell2DIndex,
                               cell2DNumVertices,
                               edgeIntersectionOne.Index,
                               edgeIntersectionTwo.Index,
                               newCell0DIndexOne,
                               newCell0DIndexTwo,
                               splitCell1DsIndexOne,
                               splitCell1DsIndexTwo,
                               cell2DEdgesDirection.at(edgeIntersectionOne.Index),
                               cell2DEdgesDirection.at(edgeIntersectionTwo.Index),
                               mesh);
    }


  }
  // ***************************************************************************
}
