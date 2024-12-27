#include "RefinementUtilities.hpp"
#include "VTKUtilities.hpp"

namespace Gedim
{
// ***************************************************************************
RefinementUtilities::RefinementUtilities(const GeometryUtilities &geometryUtilities, const MeshUtilities &meshUtilities)
    : geometryUtilities(geometryUtilities), meshUtilities(meshUtilities)
{
}
RefinementUtilities::~RefinementUtilities()
{
}
// ***************************************************************************
RefinementUtilities::SplitCell1D_Result RefinementUtilities::SplitCell1D(const unsigned int &cell1DIndex,
                                                                         const Eigen::Vector3d &newVertexCoordinate,
                                                                         IMeshDAO &mesh) const
{
    SplitCell1D_Result result;

    const unsigned int cell1DOriginIndex = mesh.Cell1DOrigin(cell1DIndex);
    const unsigned int cell1DEndIndex = mesh.Cell1DEnd(cell1DIndex);

    result.NewCell0DIndex = mesh.Cell0DAppend(1);
    const unsigned int &newCell0DIndex = result.NewCell0DIndex;
    mesh.Cell0DInsertCoordinates(newCell0DIndex, newVertexCoordinate);
    mesh.Cell0DSetMarker(newCell0DIndex, mesh.Cell1DMarker(cell1DIndex));
    mesh.Cell0DSetState(newCell0DIndex, true);

    // Split max edge into sub-edges
    Eigen::MatrixXi newCell1DsExtreme(2, 2);
    newCell1DsExtreme.col(0) << cell1DOriginIndex, newCell0DIndex;
    newCell1DsExtreme.col(1) << newCell0DIndex, cell1DEndIndex;
    result.NewCell1DsIndex = meshUtilities.SplitCell1D(cell1DIndex, newCell1DsExtreme, mesh);

    return result;
}
// ***************************************************************************
unsigned int RefinementUtilities::UpdateCell2D_NewVertex(const unsigned int cell2DIndex,
                                                         const bool cell2DEdgeDirection,
                                                         const unsigned int cell2DEdgePosition,
                                                         const std::vector<unsigned int> &newCell1DsIndex,
                                                         const unsigned int newCell0DIndex,
                                                         IMeshDAO &mesh) const
{
    const unsigned int numVertices = mesh.Cell2DNumberVertices(cell2DIndex);
    const std::vector<unsigned int> originalVertices = mesh.Cell2DVertices(cell2DIndex);
    const std::vector<unsigned int> originalEdges = mesh.Cell2DEdges(cell2DIndex);

    // create new face 2D
    Eigen::MatrixXi newFace(2, numVertices + 1);

    for (unsigned int e = 0; e < cell2DEdgePosition; e++)
        newFace.col(e) << originalVertices.at(e), originalEdges.at(e);

    if (cell2DEdgeDirection)
    {
        newFace.col(cell2DEdgePosition) << originalVertices.at(cell2DEdgePosition), newCell1DsIndex.at(0);
        newFace.col(cell2DEdgePosition + 1) << newCell0DIndex, newCell1DsIndex.at(1);
    }
    else
    {
        newFace.col(cell2DEdgePosition) << originalVertices.at(cell2DEdgePosition), newCell1DsIndex.at(1);
        newFace.col(cell2DEdgePosition + 1) << newCell0DIndex, newCell1DsIndex.at(0);
    }

    for (unsigned int e = cell2DEdgePosition + 1; e < numVertices; e++)
        newFace.col(e + 1) << originalVertices.at(e), originalEdges.at(e);

    const std::vector<unsigned int> newCell2DIndices = meshUtilities.SplitCell2D(cell2DIndex, {newFace}, mesh);

    Gedim::Output::Assert(newCell2DIndices.size() == 1);
    return newCell2DIndices.at(0);
}
// ***************************************************************************
bool RefinementUtilities::AreVerticesAligned(const Eigen::MatrixXd &cell2DVertices, const unsigned int fromVertex, const unsigned int toVertex) const
{
    const unsigned int cell2DNumVertices = cell2DVertices.cols();

    // check contigous vertices
    if ((fromVertex + 1) % cell2DNumVertices == toVertex || (toVertex + 1) % cell2DNumVertices == fromVertex)
    {
        return true;
    }

    // check aligned vertices
    if (geometryUtilities.PointIsAligned(cell2DVertices.col(fromVertex),
                                         cell2DVertices.col(toVertex),
                                         cell2DVertices.col((fromVertex + 1) % cell2DNumVertices)) ||
        geometryUtilities.PointIsAligned(cell2DVertices.col(fromVertex),
                                         cell2DVertices.col(toVertex),
                                         cell2DVertices.col((toVertex + 1) % cell2DNumVertices)))
    {
        return true;
    }

    return false;
}
// ***************************************************************************
} // namespace Gedim
