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
  void RefinementUtilities::SplitTriangle_FromFourVertices(const unsigned int& cell2DIndex,
                                                           const std::vector<unsigned int>& splitCell1DsIndex,
                                                           const bool& cell1DDirection,
                                                           const unsigned int& oppositeVertexIndex,
                                                           const unsigned int& cell0DOppositeIndex,
                                                           const unsigned int& newCell0DIndex,
                                                           IMeshDAO& mesh) const
  {
    // Create new cell1D from new vertex to opposite vertex
    const unsigned int newCell1DIndex = mesh.Cell1DAppend(1);
    mesh.Cell1DInsertExtremes(newCell1DIndex, newCell0DIndex, cell0DOppositeIndex);
    mesh.Cell1DSetMarker(newCell1DIndex, 0);
    mesh.Cell1DSetState(newCell1DIndex, true);
    mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex, 2);

    // Split cell2D into sub-cells
    std::vector<Eigen::MatrixXi> subCells(2);

    subCells[0].resize(2, 3);
    subCells[0].row(0)<< cell0DOppositeIndex,
        mesh.Cell2DVertex(cell2DIndex, (oppositeVertexIndex + 1) % 3),
        newCell0DIndex;
    subCells[0].row(1)<< mesh.Cell2DEdge(cell2DIndex, oppositeVertexIndex),
        cell1DDirection ? splitCell1DsIndex.at(0) : splitCell1DsIndex.at(1),
        newCell1DIndex;

    subCells[1].resize(2, 3);
    subCells[1].row(0)<< cell0DOppositeIndex,
        newCell0DIndex,
        mesh.Cell2DVertex(cell2DIndex, (oppositeVertexIndex + 2) % 3);
    subCells[1].row(1)<< newCell1DIndex,
        cell1DDirection ? splitCell1DsIndex.at(1) : splitCell1DsIndex.at(0),
        mesh.Cell2DEdge(cell2DIndex, (oppositeVertexIndex + 2) % 3);

    const std::vector<unsigned int> newCell2DsIndex = meshUtilities.SplitCell2D(cell2DIndex,
                                                                                subCells,
                                                                                mesh);

    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     0,
                                     newCell2DsIndex[1]); // right
    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     1,
                                     newCell2DsIndex[0]); // left
  }
  // ***************************************************************************
  void RefinementUtilities::RefineTriangleCellByMaxEdge(const unsigned int& cell2DIndex,
                                                        const Eigen::VectorXd& cell2DEdgesLength,
                                                        const std::vector<bool>& cell2DEdgesDirection,
                                                        IMeshDAO& mesh) const
  {
    Gedim::Output::Assert(cell2DEdgesLength.size() == 3);

    if (mesh.Cell2DHasUpdatedCell2Ds(cell2DIndex))
      return;

    Eigen::VectorXd::Index maxEdgeLocalIndex;
    cell2DEdgesLength.maxCoeff(&maxEdgeLocalIndex);
    const unsigned int oppositeVertexIndex = (maxEdgeLocalIndex + 2) % 3;

    // Add new mesh vertex, middle point of edge
    const unsigned int cell0DOppositeIndex = mesh.Cell2DVertex(cell2DIndex, oppositeVertexIndex);
    const unsigned int cell1DMaxIndex = mesh.Cell2DEdge(cell2DIndex, maxEdgeLocalIndex);
    const bool cell1DDirection = cell2DEdgesDirection[maxEdgeLocalIndex];
    const unsigned int cell1DOriginIndex = mesh.Cell1DOrigin(cell1DMaxIndex);
    const unsigned int cell1DEndIndex = mesh.Cell1DEnd(cell1DMaxIndex);

    const unsigned int newCell0DIndex = mesh.Cell0DAppend(1);
    mesh.Cell0DInsertCoordinates(newCell0DIndex,
                                 0.5 * (mesh.Cell0DCoordinates(cell1DOriginIndex) + mesh.Cell0DCoordinates(cell1DEndIndex)));
    mesh.Cell0DSetMarker(newCell0DIndex, mesh.Cell1DMarker(cell1DMaxIndex));
    mesh.Cell0DSetState(newCell0DIndex, true);

    // Split max edge into sub-edges
    Eigen::MatrixXi newCell1DsExtreme(2, 3);
    newCell1DsExtreme.col(0)<< cell1DOriginIndex, newCell0DIndex;
    newCell1DsExtreme.col(1)<< newCell0DIndex, cell1DEndIndex;
    const std::vector<unsigned int> splitCell1DsIndex = meshUtilities.SplitCell1D(cell1DMaxIndex,
                                                                                  newCell1DsExtreme,
                                                                                  mesh);

    SplitTriangle_FromFourVertices(cell2DIndex,
                                   splitCell1DsIndex,
                                   cell1DDirection,
                                   oppositeVertexIndex,
                                   cell0DOppositeIndex,
                                   newCell0DIndex,
                                   mesh);

    // update neighbour cell
    const unsigned int neighIndex = cell1DDirection ? 0 : 1;
    if (!mesh.Cell1DHasNeighbourCell2D(cell1DMaxIndex, neighIndex))
      return;

    const unsigned int neighCell2DIndex = mesh.Cell1DNeighbourCell2D(cell1DMaxIndex, neighIndex);
    unsigned int neighMaxEdgeIndex = 3;
    for (unsigned int e = 0; e < 3; e++)
    {
      if (mesh.Cell2DEdge(neighCell2DIndex, e) == cell1DMaxIndex)
      {
        neighMaxEdgeIndex = e;
        break;
      }
    }

    Gedim::Output::Assert(neighMaxEdgeIndex < 3);
    const unsigned int neighOppositeVertexIndex = (neighMaxEdgeIndex + 2) % 3;
    const unsigned int neighCell0DOppositeIndex = mesh.Cell2DVertex(neighCell2DIndex, neighOppositeVertexIndex);

    SplitTriangle_FromFourVertices(neighCell2DIndex,
                                   splitCell1DsIndex,
                                   !cell1DDirection,
                                   neighOppositeVertexIndex,
                                   neighCell0DOppositeIndex,
                                   newCell0DIndex,
                                   mesh);
  }
  // ***************************************************************************
}
