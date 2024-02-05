#include "RefinementUtilities.hpp"
#include "VTKUtilities.hpp"


namespace Gedim
{
  // ***************************************************************************
  RefinementUtilities::TetrahedronMaxEdgeDirection RefinementUtilities::ComputeTetrahedronMaxEdgeDirection(const Eigen::MatrixXi& polyhedronEdges,
                                                                                                           const Eigen::VectorXd& edgesLength) const
  {
    Gedim::Output::Assert(polyhedronEdges.cols() == 6);
    Gedim::Output::Assert(edgesLength.size() == 6);

    TetrahedronMaxEdgeDirection result;

    Eigen::VectorXd::Index maxEdgeLocalIndex;
    edgesLength.maxCoeff(&maxEdgeLocalIndex);

    result.MaxEdgeIndex = maxEdgeLocalIndex;

    const unsigned int edgeOrigin = polyhedronEdges(0, result.MaxEdgeIndex);
    const unsigned int edgeEnd = polyhedronEdges(1, result.MaxEdgeIndex);

    unsigned int vertexNumber = 0, vertexFound = 0;
    while (vertexFound < 2 &&
           vertexNumber < 4)
    {
      if (vertexNumber != edgeOrigin &&
          vertexNumber != edgeEnd)
      {
        result.OppositeVerticesIndex[vertexFound++] = vertexNumber;
      }

      vertexNumber++;
    }

    Gedim::Output::Assert(vertexFound == 2);

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolyhedron_Result RefinementUtilities::RefinePolyhedronCell_ByPlane(const unsigned int& cell3DIndex,
                                                                                                 const Eigen::MatrixXd& cell3DVertices,
                                                                                                 const Eigen::MatrixXi& cell3DEdges,
                                                                                                 const Eigen::VectorXd& cell3DEdgesLength,
                                                                                                 const std::vector<Eigen::MatrixXi>& cell3DFaces,
                                                                                                 const std::vector<Eigen::MatrixXd>& cell3DFaces3DVertices,
                                                                                                 const std::vector<Eigen::MatrixXd>& cell3DFacesEdges3DTangent,
                                                                                                 const std::vector<Eigen::Vector3d>& cell3DFacesTranslation,
                                                                                                 const std::vector<Eigen::Matrix3d>& cell3DFacesRotationMatrix,
                                                                                                 const double& cell3DVolume,
                                                                                                 const Eigen::Vector3d& planeNormal,
                                                                                                 const Eigen::Vector3d& planeOrigin,
                                                                                                 const Eigen::Matrix3d& planeRotationMatrix,
                                                                                                 const Eigen::Vector3d& planeTranslation,
                                                                                                 IMeshDAO& mesh) const
  {
    RefinePolyhedron_Result result;

    if (mesh.Cell3DHasUpdatedCell3Ds(cell3DIndex))
    {
      result.ResultType = RefinePolyhedron_Result::ResultTypes::Cell3DAlreadySplitted;
      return result;
    }

    // Check if the cell refinement creates null cell3Ds
    if (geometryUtilities.IsValueZero(0.5 * cell3DVolume,
                                      geometryUtilities.Tolerance3D()))
    {
      result.ResultType = RefinePolyhedron_Result::ResultTypes::Cell3DSplitUnderTolerance;
      return result;
    }

    const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult split_cell = geometryUtilities.SplitPolyhedronWithPlane(cell3DVertices,
                                                                                                                           cell3DEdges,
                                                                                                                           cell3DFaces,
                                                                                                                           cell3DFaces3DVertices,
                                                                                                                           cell3DFacesEdges3DTangent,
                                                                                                                           cell3DFacesTranslation,
                                                                                                                           cell3DFacesRotationMatrix,
                                                                                                                           planeNormal,
                                                                                                                           planeOrigin,
                                                                                                                           planeRotationMatrix,
                                                                                                                           planeTranslation);


    switch (split_cell.Type)
    {
      case Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::None:
        result.ResultType = RefinePolyhedron_Result::ResultTypes::Cell3DSplitNone;
        return result;
      case Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::Split:
      {
        // check new edges length
        for (unsigned int nv = 0; nv < split_cell.Vertices.NewVerticesOriginalEdge.size(); nv++)
        {
          const unsigned int originalEdgeIndex = split_cell.Vertices.NewVerticesOriginalEdge.at(nv);

          if (geometryUtilities.IsValueZero(0.5 * cell3DEdgesLength[originalEdgeIndex],
                                            geometryUtilities.Tolerance1D()))
          {
            result.ResultType = RefinePolyhedron_Result::ResultTypes::Cell3DSplitUnderTolerance;
            return result;
          }
        }
      }
        break;
      default:
        throw std::runtime_error("Split polyhedron unknown");
    }

    // Create new mesh elements
    const unsigned int numOriginalVertices = cell3DVertices.cols();
    const unsigned int numNewVertices = split_cell.Vertices.NewVerticesOriginalEdge.size();

    std::vector<unsigned int> splitCell0DsIndex(split_cell.Vertices.Vertices.cols(),
                                                std::numeric_limits<unsigned int>::max());
    std::vector<unsigned int> splitCell1DsIndex(split_cell.Edges.NewEdgesOriginalEdges.size(),
                                                std::numeric_limits<unsigned int>::max());
    std::vector<unsigned int> splitCell2DsIndex(split_cell.Faces.NewFacesOriginalFaces.size(),
                                                std::numeric_limits<unsigned int>::max());

    std::list<RefinePolyhedron_Result::RefinedCell1D> newCell1Ds;
    std::list<RefinePolyhedron_Result::RefinedCell2D> newCell2Ds;
    std::unordered_map<unsigned int, unsigned int> edgeIndexToNewCell1Ds;

    for (unsigned int v = 0; v < numOriginalVertices; v++)
      splitCell0DsIndex[v] = mesh.Cell3DVertex(cell3DIndex, v);

    // Create new Cell0Ds
    result.NewCell0DsIndex.resize(numNewVertices);
    for (unsigned int nv = 0; nv < numNewVertices; nv++)
    {
      const unsigned int originalEdgeIndex = split_cell.Vertices.NewVerticesOriginalEdge.at(nv);
      const unsigned int originalCell1DIndex = mesh.Cell3DEdge(cell3DIndex,
                                                               originalEdgeIndex);

      const SplitCell1D_Result splitResult = SplitCell1D(originalCell1DIndex,
                                                         split_cell.Vertices.Vertices.col(numOriginalVertices + nv),
                                                         mesh);

      result.NewCell0DsIndex[nv] = splitResult.NewCell0DIndex;
      splitCell0DsIndex[numOriginalVertices + nv] = splitResult.NewCell0DIndex;

      RefinePolyhedron_Result::RefinedCell1D refinedCell1D;
      refinedCell1D.Type = RefinePolyhedron_Result::RefinedCell1D::Types::Updated;
      refinedCell1D.NewCell1DsIndex = splitResult.NewCell1DsIndex;
      refinedCell1D.OriginalCell1DIndex = originalCell1DIndex;
      refinedCell1D.NewCell0DIndex = splitResult.NewCell0DIndex;
      refinedCell1D.OriginalCell3DEdgeIndex = originalEdgeIndex;
      const std::vector<unsigned int>& newEdges = split_cell.OriginalEdgesNewEdges[originalEdgeIndex];

      const unsigned int newEdgeOneCell0DOrigin = splitCell0DsIndex.at(split_cell.Edges.Edges(0, newEdges[0]));
      const unsigned int newEdgeOneCell0DEnd = splitCell0DsIndex.at(split_cell.Edges.Edges(1, newEdges[0]));

      const unsigned int newCell1DOneOrigin = mesh.Cell1DOrigin(refinedCell1D.NewCell1DsIndex[0]);
      const unsigned int newCell1DOneEnd = mesh.Cell1DEnd(refinedCell1D.NewCell1DsIndex[0]);
      const unsigned int newCell1DTwoOrigin = mesh.Cell1DOrigin(refinedCell1D.NewCell1DsIndex[1]);
      const unsigned int newCell1DTwoEnd = mesh.Cell1DEnd(refinedCell1D.NewCell1DsIndex[1]);

      if ((newEdgeOneCell0DOrigin == newCell1DOneOrigin &&
           newEdgeOneCell0DEnd == newCell1DOneEnd) ||
          (newEdgeOneCell0DEnd == newCell1DOneOrigin &&
           newEdgeOneCell0DOrigin == newCell1DOneEnd))
      {
        splitCell1DsIndex[newEdges[0]] = refinedCell1D.NewCell1DsIndex[0];
        splitCell1DsIndex[newEdges[1]] = refinedCell1D.NewCell1DsIndex[1];
      }
      else if ((newEdgeOneCell0DOrigin == newCell1DTwoOrigin &&
                newEdgeOneCell0DEnd == newCell1DTwoEnd) ||
               (newEdgeOneCell0DEnd == newCell1DTwoOrigin &&
                newEdgeOneCell0DOrigin == newCell1DTwoEnd))
      {
        splitCell1DsIndex[newEdges[0]] = refinedCell1D.NewCell1DsIndex[1];
        splitCell1DsIndex[newEdges[1]] = refinedCell1D.NewCell1DsIndex[0];
      }
      else
        throw std::runtime_error("Error on edge identification");


      edgeIndexToNewCell1Ds.insert(std::make_pair(newEdges[0],
                                   newCell1Ds.size()));
      edgeIndexToNewCell1Ds.insert(std::make_pair(newEdges[1],
                                   newCell1Ds.size()));
      newCell1Ds.push_back(refinedCell1D);
    }

    // Create new Cell1Ds
    for (unsigned int ne = 0; ne < split_cell.Edges.NewEdgesOriginalEdges.size(); ne++)
    {
      if (split_cell.Edges.NewEdgesOriginalEdges[ne] >= 0)
      {
        if (splitCell1DsIndex[ne] != std::numeric_limits<unsigned int>::max())
          continue;

        splitCell1DsIndex[ne] = mesh.Cell3DEdge(cell3DIndex,
                                                split_cell.Edges.NewEdgesOriginalEdges[ne]);
        continue;
      }

      const unsigned int originalFaceIndex = split_cell.Edges.NewEdgesOriginalFace.at(ne);
      const unsigned int originalCell2DIndex = mesh.Cell3DFace(cell3DIndex,
                                                               originalFaceIndex);

      const unsigned int newEdgeOrigin = split_cell.Edges.Edges(0, ne);
      const unsigned int newEdgeEnd = split_cell.Edges.Edges(1, ne);

      const unsigned int cell0DIndexOrigin = splitCell0DsIndex.at(newEdgeOrigin);
      const unsigned int cell0DIndexEnd = splitCell0DsIndex.at(newEdgeEnd);

      RefinePolyhedron_Result::RefinedCell1D newCell1D;
      newCell1D.Type = RefinePolyhedron_Result::RefinedCell1D::Types::New;
      newCell1D.NewCell1DsIndex.resize(1);

      newCell1D.NewCell1DsIndex[0] = mesh.Cell1DAppend(1);
      const unsigned int newCell1DIndex = newCell1D.NewCell1DsIndex[0];
      splitCell1DsIndex[ne] = newCell1DIndex;

      mesh.Cell1DInsertExtremes(newCell1DIndex,
                                cell0DIndexEnd,
                                cell0DIndexOrigin);
      mesh.Cell1DSetMarker(newCell1DIndex,
                           mesh.Cell2DMarker(originalCell2DIndex));
      mesh.Cell1DSetState(newCell1DIndex, true);
      mesh.Cell1DInitializeNeighbourCell3Ds(newCell1DIndex,
                                            3);
      for (unsigned int nf = 0; nf < mesh.Cell2DNumberNeighbourCell3D(originalCell2DIndex); nf++)
      {
        if (!mesh.Cell2DHasNeighbourCell3D(originalCell2DIndex, nf))
          continue;

        const unsigned int neighCell3DIndex = mesh.Cell2DNeighbourCell3D(originalCell2DIndex, nf);

        if (neighCell3DIndex == cell3DIndex)
          continue;

        mesh.Cell1DInsertNeighbourCell3D(newCell1DIndex,
                                         0,
                                         neighCell3DIndex);
      }

      edgeIndexToNewCell1Ds.insert(std::make_pair(ne,
                                                  newCell1Ds.size()));
      newCell1Ds.push_back(newCell1D);
    }

    // Create new Cell2Ds
    for (unsigned int ne = 0; ne < split_cell.Edges.NewEdgesOriginalEdges.size(); ne++)
    {
      if (split_cell.Edges.NewEdgesOriginalEdges[ne] >= 0)
        continue;

      const unsigned int originalFaceIndex = split_cell.Edges.NewEdgesOriginalFace.at(ne);
      const unsigned int originalCell2DIndex = mesh.Cell3DFace(cell3DIndex,
                                                               originalFaceIndex);

      const std::vector<unsigned int>& newFaces = split_cell.OriginalFacesNewFaces[originalFaceIndex];

      const Eigen::MatrixXi& newFaceOne = split_cell.Faces.Faces[newFaces[0]];
      const Eigen::MatrixXi& newFaceTwo = split_cell.Faces.Faces[newFaces[1]];

      std::vector<Eigen::MatrixXi> subCells(2);
      subCells[0].resize(2, newFaceOne.cols());
      subCells[1].resize(2, newFaceTwo.cols());

      std::set<unsigned int> newCell1DsPosition;

      for (unsigned int v = 0; v < newFaceOne.cols(); v++)
      {
        subCells[0](0, v) = splitCell0DsIndex.at(newFaceOne(0, v));
        subCells[0](1, v) = splitCell1DsIndex.at(newFaceOne(1, v));

        const auto& edgeIndexToNewCell1D = edgeIndexToNewCell1Ds.find(newFaceOne(1, v));
        if (edgeIndexToNewCell1D !=
            edgeIndexToNewCell1Ds.end())
        {
          if (newCell1DsPosition.find(edgeIndexToNewCell1D->second) ==
              newCell1DsPosition.end())
            newCell1DsPosition.insert(edgeIndexToNewCell1D->second);
        }
      }

      std::list<unsigned int> newFaceTwoNewEdges;
      for (unsigned int v = 0; v < newFaceTwo.cols(); v++)
      {
        subCells[1](0, v) = splitCell0DsIndex.at(newFaceTwo(0, v));
        subCells[1](1, v) = splitCell1DsIndex.at(newFaceTwo(1, v));

        const auto& edgeIndexToNewCell1D = edgeIndexToNewCell1Ds.find(newFaceTwo(1, v));
        if (edgeIndexToNewCell1D !=
            edgeIndexToNewCell1Ds.end())
        {
          if (newCell1DsPosition.find(edgeIndexToNewCell1D->second) ==
              newCell1DsPosition.end())
            newCell1DsPosition.insert(edgeIndexToNewCell1D->second);
        }
      }

      const std::vector<unsigned int> newCell2DIndices = meshUtilities.SplitCell2D(originalCell2DIndex,
                                                                                   subCells,
                                                                                   mesh);

      RefinePolyhedron_Result::RefinedCell2D refinedCell2D;
      refinedCell2D.Type = RefinePolyhedron_Result::RefinedCell2D::Types::Updated;
      refinedCell2D.NewCell2DsIndex = newCell2DIndices;
      refinedCell2D.OriginalCell2DIndex = originalCell2DIndex;
      refinedCell2D.NewCell1DIndex = splitCell1DsIndex[ne];
      refinedCell2D.OriginalCell3DFaceIndex = originalFaceIndex;
      refinedCell2D.NewCell1DsPosition = std::vector<unsigned int>(newCell1DsPosition.begin(),
                                                                   newCell1DsPosition.end());

      splitCell2DsIndex[newFaces[0]] = newCell2DIndices[0];
      splitCell2DsIndex[newFaces[1]] = newCell2DIndices[1];

      newCell2Ds.push_back(refinedCell2D);
    }

    for (unsigned int nf = 0; nf < split_cell.Faces.NewFacesOriginalFaces.size(); nf++)
    {
      if (split_cell.Faces.NewFacesOriginalFaces[nf] >= 0)
      {
        if (splitCell2DsIndex[nf] != std::numeric_limits<unsigned int>::max())
          continue;

        splitCell2DsIndex[nf] = mesh.Cell3DFace(cell3DIndex,
                                                split_cell.Faces.NewFacesOriginalFaces[nf]);
        continue;
      }

      const Eigen::MatrixXi newFace = split_cell.Faces.Faces[nf];

      Eigen::MatrixXi cell2DExtremes(2, newFace.cols());

      for (unsigned int nfv = 0; nfv < newFace.cols(); nfv++)
      {
        cell2DExtremes(0, nfv) = splitCell0DsIndex.at(newFace(0, nfv));
        cell2DExtremes(1, nfv) = splitCell1DsIndex.at(newFace(1, nfv));

        if (split_cell.Edges.NewEdgesOriginalEdges.at(newFace(1, nfv)) >= 0)
        {
          const unsigned int cell1DIndex = cell2DExtremes(1, nfv);
          const unsigned int numCell1DNeigh3Ds = mesh.Cell1DNumberNeighbourCell3D(cell1DIndex);
          mesh.Cell1DInitializeNeighbourCell3Ds(cell1DIndex,
                                                numCell1DNeigh3Ds + 1);
        }
      }

      RefinePolyhedron_Result::RefinedCell2D newCell2D;
      newCell2D.Type = RefinePolyhedron_Result::RefinedCell2D::Types::New;
      newCell2D.NewCell2DsIndex.resize(1);
      newCell2D.NewCell2DsIndex[0] = mesh.Cell2DAppend(1);
      const unsigned int newCell2DIndex = newCell2D.NewCell2DsIndex[0];
      splitCell2DsIndex[nf] = newCell2DIndex;
      mesh.Cell2DAddVerticesAndEdges(newCell2DIndex,
                                     cell2DExtremes);

      mesh.Cell2DSetMarker(newCell2DIndex, 0);
      mesh.Cell2DSetState(newCell2DIndex, true);
      mesh.Cell2DInitializeNeighbourCell3Ds(newCell2DIndex,
                                            2);

      newCell2Ds.push_back(newCell2D);
    }

    // Create new cell3Ds
    std::vector<std::vector<unsigned int>> newCell3DsVertices(2);
    std::vector<std::vector<unsigned int>> newCell3DsEdges(2);
    std::vector<std::vector<unsigned int>> newCell3DsFaces(2);

    for (unsigned int nc = 0; nc < 2; nc++)
    {
      const Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::NewPolyhedron newPolyhedron =
          (nc == 0) ?
            split_cell.PositivePolyhedron :
            split_cell.NegativePolyhedron;

      newCell3DsVertices[nc].resize(newPolyhedron.Vertices.size());
      newCell3DsEdges[nc].resize(newPolyhedron.Edges.size());
      newCell3DsFaces[nc].resize(newPolyhedron.Faces.size());

      for (unsigned int v = 0; v < newPolyhedron.Vertices.size(); v++)
        newCell3DsVertices[nc][v] = splitCell0DsIndex.at(newPolyhedron.Vertices.at(v));
      for (unsigned int e = 0; e < newPolyhedron.Edges.size(); e++)
        newCell3DsEdges[nc][e] = splitCell1DsIndex.at(newPolyhedron.Edges.at(e));
      for (unsigned int f = 0; f < newPolyhedron.Faces.size(); f++)
        newCell3DsFaces[nc][f] = splitCell2DsIndex.at(newPolyhedron.Faces.at(f));
    }

    const std::vector<unsigned int> newCell3DIndices = meshUtilities.SplitCell3D(cell3DIndex,
                                                                                 newCell3DsVertices,
                                                                                 newCell3DsEdges,
                                                                                 newCell3DsFaces,
                                                                                 mesh);

    for (unsigned int nf = 0; nf < split_cell.Faces.NewFacesOriginalFaces.size(); nf++)
    {
      if (split_cell.Faces.NewFacesOriginalFaces[nf] >= 0)
        continue;

      const unsigned newCell2DIndex = splitCell2DsIndex.at(nf);
      mesh.Cell2DInsertNeighbourCell3D(newCell2DIndex,
                                       0,
                                       newCell3DIndices[0]);
      mesh.Cell2DInsertNeighbourCell3D(newCell2DIndex,
                                       1,
                                       newCell3DIndices[1]);

      const Eigen::MatrixXi newFace = split_cell.Faces.Faces[nf];

      for (unsigned int nfe = 0; nfe < newFace.cols(); nfe++)
      {
        if (split_cell.Edges.NewEdgesOriginalEdges.at(newFace(1, nfe)) >= 0)
        {
          const unsigned int cell1DIndex = splitCell1DsIndex.at(newFace(1, nfe));

          const unsigned int numCell1DNeigh3Ds = mesh.Cell1DNumberNeighbourCell3D(cell1DIndex);
          bool hasPositiveCell3D = false;

          for (unsigned int nfen = 0; nfen < numCell1DNeigh3Ds - 1; nfen++)
          {
            if (!mesh.Cell1DHasNeighbourCell3D(cell1DIndex, nfen))
              continue;

            if (mesh.Cell1DNeighbourCell3D(cell1DIndex, nfen) == newCell3DIndices.at(0))
            {
              hasPositiveCell3D = true;
              break;
            }
          }

          mesh.Cell1DInsertNeighbourCell3D(cell1DIndex,
                                           numCell1DNeigh3Ds - 1,
                                           hasPositiveCell3D ?
                                             newCell3DIndices.at(1) :
                                             newCell3DIndices.at(0));

        }
      }

    }

    for (const RefinePolyhedron_Result::RefinedCell1D& newCell1D : newCell1Ds)
    {
      if (newCell1D.Type != RefinePolyhedron_Result::RefinedCell1D::Types::New)
        continue;

      const unsigned int cell1DIndex = newCell1D.NewCell1DsIndex[0];

      mesh.Cell1DInsertNeighbourCell3D(cell1DIndex,
                                       1,
                                       newCell3DIndices[0]);
      mesh.Cell1DInsertNeighbourCell3D(cell1DIndex,
                                       2,
                                       newCell3DIndices[1]);
    }

    result.NewCell1DsIndex = std::vector<RefinePolyhedron_Result::RefinedCell1D>(newCell1Ds.begin(),
                                                                                 newCell1Ds.end());
    result.NewCell2DsIndex = std::vector<RefinePolyhedron_Result::RefinedCell2D>(newCell2Ds.begin(),
                                                                                 newCell2Ds.end());

    result.NewCell3DsIndex = newCell3DIndices;
    result.ResultType = RefinePolyhedron_Result::ResultTypes::Successfull;

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolyhedron_UpdateNeighbour_Result RefinementUtilities::RefinePolyhedronCell_UpdateFaceNeighbours(const unsigned int& cell3DIndex,
                                                                                                                              const unsigned int& cell2DIndex,
                                                                                                                              const unsigned int& newCell1DIndex,
                                                                                                                              const std::vector<unsigned int>& splitCell1DsOriginalIndex,
                                                                                                                              const std::vector<unsigned int>& splitCell1DsNewCell0DIndex,
                                                                                                                              const std::vector<std::vector<unsigned int>>& splitCell1DsUpdatedIndices,
                                                                                                                              const std::vector<unsigned int>& splitCell2DsIndex,
                                                                                                                              const std::vector<std::vector<std::vector<bool>>>& cell3DsFacesEdgesDirection,
                                                                                                                              std::map<unsigned int, unsigned int>& updatedCell2Ds,
                                                                                                                              IMeshDAO& mesh) const
  {
    RefinePolyhedron_UpdateNeighbour_Result result;

    std::list<RefinePolyhedron_UpdateNeighbour_Result::UpdatedCell3D> newCell3DsIndex;


    // update neighbour cells
    for (unsigned int n = 0; n < mesh.Cell2DNumberNeighbourCell3D(cell2DIndex); n++)
    {
      if (!mesh.Cell2DHasNeighbourCell3D(cell2DIndex, n))
        continue;

      const unsigned int neighCell3DIndex = mesh.Cell2DNeighbourCell3D(cell2DIndex, n);

      if (neighCell3DIndex == cell3DIndex)
        continue;

      if (mesh.Cell3DHasUpdatedCell3Ds(neighCell3DIndex))
        continue;

      // update new cell3D
      std::vector<std::vector<unsigned int>> newCell3DsVertices(1);
      std::vector<std::vector<unsigned int>> newCell3DsEdges(1);
      std::vector<std::vector<unsigned int>> newCell3DsFaces(1);

      const std::vector<unsigned int> originalVertices = mesh.Cell3DVertices(neighCell3DIndex);
      const std::vector<unsigned int> originalEdges = mesh.Cell3DEdges(neighCell3DIndex);
      const std::vector<unsigned int> originalFaces = mesh.Cell3DFaces(neighCell3DIndex);

      const unsigned int neighFaceIndex = mesh.Cell3DFindFace(neighCell3DIndex,
                                                              cell2DIndex);
      Gedim::Output::Assert(neighFaceIndex < originalFaces.size());

      newCell3DsVertices[0].resize(originalVertices.size() + splitCell1DsOriginalIndex.size());
      newCell3DsEdges[0].resize(originalEdges.size() + splitCell1DsOriginalIndex.size() + 1);
      newCell3DsFaces[0].resize(originalFaces.size() + 1);

      std::copy(originalVertices.begin(), originalVertices.end(), newCell3DsVertices[0].begin());
      std::copy(originalEdges.begin(), originalEdges.end(), newCell3DsEdges[0].begin());
      std::copy(originalFaces.begin(), originalFaces.end(), newCell3DsFaces[0].begin());

      for (unsigned int ne = 0; ne < splitCell1DsOriginalIndex.size(); ne++)
      {
        const unsigned int cell1DIndex = splitCell1DsOriginalIndex[ne];
        const unsigned int neighEdgeIndex = mesh.Cell3DFindEdge(neighCell3DIndex,
                                                                cell1DIndex);

        Gedim::Output::Assert(neighEdgeIndex < originalEdges.size());

        newCell3DsEdges[0][neighEdgeIndex] = splitCell1DsUpdatedIndices[ne][0];
        newCell3DsEdges[0][originalEdges.size() + ne] = splitCell1DsUpdatedIndices[ne][1];
        newCell3DsVertices[0][originalVertices.size() + ne] = splitCell1DsNewCell0DIndex[ne];

        // update faces with no new edges
        for (unsigned int nef = 0; nef < originalFaces.size(); nef++)
        {
          if (nef == neighFaceIndex)
            continue;

          const unsigned int cell1DCell2DIndex = originalFaces.at(nef);

          const unsigned int localFaceEdgeIndex = mesh.Cell2DFindEdge(cell1DCell2DIndex,
                                                                      cell1DIndex);
          if (localFaceEdgeIndex == mesh.Cell2DNumberVertices(cell1DCell2DIndex))
            continue;

          const unsigned int localFaceIndex = mesh.Cell3DFindFace(neighCell3DIndex,
                                                                  cell1DCell2DIndex);
          Gedim::Output::Assert(localFaceIndex < originalFaces.size());

          if (updatedCell2Ds.find(cell1DCell2DIndex) ==
              updatedCell2Ds.end())
          {
            const unsigned int updatedCell2DIndex = UpdateCell2D_NewVertex(cell1DCell2DIndex,
                                                                           cell3DsFacesEdgesDirection.at(neighCell3DIndex).at(localFaceIndex).at(localFaceEdgeIndex),
                                                                           localFaceEdgeIndex,
                                                                           splitCell1DsUpdatedIndices[ne],
                                                                           splitCell1DsNewCell0DIndex[ne],
                                                                           mesh);
            newCell3DsFaces[0][localFaceIndex] = updatedCell2DIndex;
            updatedCell2Ds.insert(std::make_pair(cell1DCell2DIndex,
                                                 updatedCell2DIndex));
          }
          else
            newCell3DsFaces[0][localFaceIndex] = updatedCell2Ds.at(cell1DCell2DIndex);
        }
      }

      newCell3DsEdges[0][originalEdges.size() + splitCell1DsOriginalIndex.size()] = newCell1DIndex;

      newCell3DsFaces[0][neighFaceIndex] = splitCell2DsIndex[0];
      newCell3DsFaces[0][originalFaces.size()] = splitCell2DsIndex[1];

      const std::vector<unsigned int> newCell3DIndices = meshUtilities.SplitCell3D(neighCell3DIndex,
                                                                                   newCell3DsVertices,
                                                                                   newCell3DsEdges,
                                                                                   newCell3DsFaces,
                                                                                   mesh);

      for (const unsigned int& newCell3D : newCell3DIndices)
      {
        RefinePolyhedron_UpdateNeighbour_Result::UpdatedCell3D updatedCell3D;
        updatedCell3D.OriginalCell3DIndex = neighCell3DIndex;
        updatedCell3D.NewCell3DIndex = newCell3D;
        newCell3DsIndex.push_back(updatedCell3D);
      }
    }

    result.UpdatedCell3Ds = std::vector<RefinePolyhedron_UpdateNeighbour_Result::UpdatedCell3D>(newCell3DsIndex.begin(),
                                                                                                newCell3DsIndex.end());
    return result;
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolyhedron_UpdateNeighbour_Result RefinementUtilities::RefinePolyhedronCell_UpdateEdgeNeighbours(const unsigned int& cell3DIndex,
                                                                                                                              const unsigned int& cell1DIndex,
                                                                                                                              const std::vector<unsigned int>& newCell1DsIndex,
                                                                                                                              const unsigned int& newCell0DIndex,
                                                                                                                              const std::vector<std::vector<std::vector<bool>>>& cell3DsFacesEdgesDirection,
                                                                                                                              std::map<unsigned int, unsigned int>& updatedCell2Ds,
                                                                                                                              IMeshDAO& mesh) const
  {
    RefinePolyhedron_UpdateNeighbour_Result result;

    std::list<RefinePolyhedron_UpdateNeighbour_Result::UpdatedCell3D> newCell3DsIndex;

    // update neighbour cells
    for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell3D(cell1DIndex); n++)
    {
      if (!mesh.Cell1DHasNeighbourCell3D(cell1DIndex, n))
        continue;

      const unsigned int neighCell3DIndex = mesh.Cell1DNeighbourCell3D(cell1DIndex, n);
      if (neighCell3DIndex == cell3DIndex)
        continue;

      if (mesh.Cell3DHasUpdatedCell3Ds(neighCell3DIndex))
        continue;

      // update new cell3D
      std::vector<std::vector<unsigned int>> newCell3DsVertices(1);
      std::vector<std::vector<unsigned int>> newCell3DsEdges(1);
      std::vector<std::vector<unsigned int>> newCell3DsFaces(1);

      const std::vector<unsigned int> originalVertices = mesh.Cell3DVertices(neighCell3DIndex);
      const std::vector<unsigned int> originalEdges = mesh.Cell3DEdges(neighCell3DIndex);
      const std::vector<unsigned int> originalFaces = mesh.Cell3DFaces(neighCell3DIndex);

      const unsigned int neighEdgeIndex = mesh.Cell3DFindEdge(neighCell3DIndex,
                                                              cell1DIndex);
      Gedim::Output::Assert(neighEdgeIndex < originalEdges.size());

      newCell3DsVertices[0].resize(originalVertices.size() + 1);
      newCell3DsEdges[0].resize(originalEdges.size() + 1);
      newCell3DsFaces[0].resize(originalFaces.size());

      std::copy(originalVertices.begin(), originalVertices.end(), newCell3DsVertices[0].begin());
      std::copy(originalEdges.begin(), originalEdges.end(), newCell3DsEdges[0].begin());
      std::copy(originalFaces.begin(), originalFaces.end(), newCell3DsFaces[0].begin());

      newCell3DsVertices[0][originalVertices.size()] = newCell0DIndex;
      newCell3DsEdges[0][neighEdgeIndex] = newCell1DsIndex[0];
      newCell3DsEdges[0][originalEdges.size()] = newCell1DsIndex[1];

      // update faces with new edges
      for (unsigned int nef = 0; nef < originalFaces.size(); nef++)
      {
        const unsigned int cell2DIndex = originalFaces.at(nef);

        const unsigned int localFaceEdgeIndex = mesh.Cell2DFindEdge(cell2DIndex,
                                                                    cell1DIndex);
        if (localFaceEdgeIndex == mesh.Cell2DNumberVertices(cell2DIndex))
          continue;

        const unsigned int localFaceIndex = mesh.Cell3DFindFace(neighCell3DIndex,
                                                                cell2DIndex);
        Gedim::Output::Assert(localFaceIndex < originalFaces.size());

        if (updatedCell2Ds.find(cell2DIndex) ==
            updatedCell2Ds.end())
        {
          const unsigned int updatedCell2DIndex = UpdateCell2D_NewVertex(cell2DIndex,
                                                                         cell3DsFacesEdgesDirection.at(neighCell3DIndex).at(localFaceIndex).at(localFaceEdgeIndex),
                                                                         localFaceEdgeIndex,
                                                                         newCell1DsIndex,
                                                                         newCell0DIndex,
                                                                         mesh);
          newCell3DsFaces[0][localFaceIndex] = updatedCell2DIndex;
          updatedCell2Ds.insert(std::make_pair(cell2DIndex,
                                               updatedCell2DIndex));
        }
        else
          newCell3DsFaces[0][localFaceIndex] = updatedCell2Ds.at(cell2DIndex);
      }

      const std::vector<unsigned int> newCell3DIndices = meshUtilities.SplitCell3D(neighCell3DIndex,
                                                                                   newCell3DsVertices,
                                                                                   newCell3DsEdges,
                                                                                   newCell3DsFaces,
                                                                                   mesh);

      for (const unsigned int& newCell3D : newCell3DIndices)
      {
        RefinePolyhedron_UpdateNeighbour_Result::UpdatedCell3D updatedCell3D;
        updatedCell3D.OriginalCell3DIndex = neighCell3DIndex;
        updatedCell3D.NewCell3DIndex = newCell3D;
        newCell3DsIndex.push_back(updatedCell3D);
      }
    }

    result.UpdatedCell3Ds = std::vector<RefinePolyhedron_UpdateNeighbour_Result::UpdatedCell3D>(newCell3DsIndex.begin(),
                                                                                                newCell3DsIndex.end());


    return result;
  }
  // ***************************************************************************
}
