#include "RefinementUtilities.hpp"
#include "VTKUtilities.hpp"


namespace Gedim
{
  // ***************************************************************************
  RefinementUtilities::CheckSplitType_Result RefinementUtilities::SplitPolygon_CheckSplitType(const Gedim::GeometryUtilities::PolygonTypes& cell2DPolygonType,
                                                                                              const Gedim::GeometryUtilities::PolygonTypes& cell2DUnalignedPolygonType,
                                                                                              const Eigen::MatrixXd& cell2DVertices,
                                                                                              const RefinePolygon_CheckResult& cell2DCheckToRefine) const
  {
    CheckSplitType_Result result;

    const bool isTriangle = (cell2DPolygonType == Gedim::GeometryUtilities::PolygonTypes::Triangle);
    const bool isUnalignedTriangle = (cell2DUnalignedPolygonType == Gedim::GeometryUtilities::PolygonTypes::Triangle);

    const unsigned int cell2DNumVertices = cell2DVertices.cols();
    Gedim::Output::Assert(cell2DCheckToRefine.Cell1DsIndex.size() == 2);

    Gedim::Output::Assert(cell2DCheckToRefine.Cell1DsIntersection.size() == 2);
    const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersectionOne = cell2DCheckToRefine.Cell1DsIntersection[0];
    const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersectionTwo = cell2DCheckToRefine.Cell1DsIntersection[1];

    Gedim::Output::Assert(cell2DCheckToRefine.Cell1DsToSplit.size() == 2);
    const RefinePolygon_CheckResult::Cell1DToSplit& cell1DToSplitOne = cell2DCheckToRefine.Cell1DsToSplit[0];
    const RefinePolygon_CheckResult::Cell1DToSplit& cell1DToSplitTwo = cell2DCheckToRefine.Cell1DsToSplit[1];

    if (isTriangle &&
        cell1DToSplitTwo.Type ==
        RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside
        )
    {
      result.Type = CheckSplitType_Result::SplitTypes::NewVertexFrom;
      return result;
    }

    if (isTriangle &&
        cell1DToSplitOne.Type ==
        RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside)
    {
      result.Type = CheckSplitType_Result::SplitTypes::NewVertexTo;
      return result;
    }

    if (!isTriangle &&
        !cell1DToSplitOne.IsToSplit &&
        !cell1DToSplitTwo.IsToSplit)
    {
      unsigned int fromVertex = geometryUtilities.IsValueGreaterOrEqual(0.5, edgeIntersectionOne.CurvilinearCoordinate,
                                                                        geometryUtilities.Tolerance1D()) ?
                                  edgeIntersectionOne.Index :
                                  (edgeIntersectionOne.Index + 1) % cell2DNumVertices;
      unsigned int toVertex = geometryUtilities.IsValueGreaterOrEqual(0.5, edgeIntersectionTwo.CurvilinearCoordinate,
                                                                      geometryUtilities.Tolerance1D()) ?
                                edgeIntersectionTwo.Index :
                                (edgeIntersectionTwo.Index + 1) % cell2DNumVertices;

      // check aligned vertices
      if (AreVerticesAligned(cell2DVertices, fromVertex, toVertex))
      {
        // try to flip vertices
        if (cell1DToSplitOne.Type !=
            RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside &&
            cell1DToSplitTwo.Type !=
            RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside)
        {
          // no fixed vertices, flip toVertex
          toVertex = (toVertex == edgeIntersectionTwo.Index) ?
                       (edgeIntersectionTwo.Index + 1) % cell2DNumVertices :
                       edgeIntersectionTwo.Index;
        }
        else if (cell1DToSplitOne.Type ==
                 RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside &&
                 cell1DToSplitTwo.Type !=
                 RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside)
        {
          // fromVertex is fixed
          toVertex = (toVertex == edgeIntersectionTwo.Index) ?
                       (edgeIntersectionTwo.Index + 1) % cell2DNumVertices :
                       edgeIntersectionTwo.Index;
        }
        else if (cell1DToSplitOne.Type !=
                 RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside &&
                 cell1DToSplitTwo.Type ==
                 RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside)
        {
          // toVertex is fixed
          fromVertex = (fromVertex == edgeIntersectionOne.Index) ?
                         (edgeIntersectionOne.Index + 1) % cell2DNumVertices :
                         edgeIntersectionOne.Index;
        }


        if (AreVerticesAligned(cell2DVertices, fromVertex, toVertex))
        {
          if (isUnalignedTriangle)
          {
            if (cell1DToSplitOne.Type !=
                RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside &&
                cell1DToSplitTwo.Type ==
                RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside)
            {
              result.Type = CheckSplitType_Result::SplitTypes::NewVertexFrom;
              return result;
            }

            if (cell1DToSplitOne.Type ==
                RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside &&
                cell1DToSplitTwo.Type !=
                RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside)
            {
              result.Type = CheckSplitType_Result::SplitTypes::NewVertexTo;
              return result;
            }

            return result;
          }
        }
      }

      result.NoNewVerticesIndex = { fromVertex, toVertex };
      result.Type = CheckSplitType_Result::SplitTypes::NoNewVertices;
      return result;
    }

    if (cell1DToSplitOne.IsToSplit &&
        !cell1DToSplitTwo.IsToSplit)
    {
      result.Type = CheckSplitType_Result::SplitTypes::NewVertexFrom;
      return result;
    }

    if (!cell1DToSplitOne.IsToSplit &&
        cell1DToSplitTwo.IsToSplit)
    {
      result.Type = CheckSplitType_Result::SplitTypes::NewVertexTo;
      return result;
    }

    if (!isTriangle)
    {
      result.Type = CheckSplitType_Result::SplitTypes::NewVertices;
      return result;
    }

    // not managed case
    return result;
  }
  // ***************************************************************************
  bool RefinementUtilities::SplitPolygon_CheckIsNotToExtend(const RefinePolygon_CheckResult::Cell1DToSplit& cell1DSplitOne,
                                                            const RefinePolygon_CheckResult::Cell1DToSplit& cell1DSplitTwo) const
  {
    if ((cell1DSplitOne.Type == RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyNeighQualityNotEnough ||
         cell1DSplitOne.Type == RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyNeighAlignedNotRespect) &&
        (cell1DSplitTwo.Type != RefinePolygon_CheckResult::Cell1DToSplit::Types::BothQualityNotEnough &&
         cell1DSplitTwo.Type != RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyLocalQualityNotEnough &&
         cell1DSplitTwo.Type != RefinePolygon_CheckResult::Cell1DToSplit::Types::BothAlignedNotRespect &&
         cell1DSplitTwo.Type != RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyLocalAlignedNotRespect)
        )
      return false;

    if ((cell1DSplitTwo.Type == RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyNeighQualityNotEnough ||
         cell1DSplitTwo.Type == RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyNeighAlignedNotRespect) &&
        (cell1DSplitOne.Type != RefinePolygon_CheckResult::Cell1DToSplit::Types::BothQualityNotEnough &&
         cell1DSplitOne.Type != RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyLocalQualityNotEnough &&
         cell1DSplitOne.Type != RefinePolygon_CheckResult::Cell1DToSplit::Types::BothAlignedNotRespect &&
         cell1DSplitOne.Type != RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyLocalAlignedNotRespect)
        )
      return false;

    return true;
  }
  // ***************************************************************************
  bool RefinementUtilities::SplitPolygon_CheckIsToSplit_Relaxed(const RefinePolygon_CheckResult::Cell1DToSplit& cell1DSplitOne,
                                                                const RefinePolygon_CheckResult::Cell1DToSplit& cell1DSplitTwo) const
  {
    if ((cell1DSplitOne.Type == RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyNeighQualityNotEnough ||
         cell1DSplitOne.Type == RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyNeighAlignedNotRespect))
      return false;

    if ((cell1DSplitTwo.Type == RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyNeighQualityNotEnough ||
         cell1DSplitTwo.Type == RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyNeighAlignedNotRespect))
      return false;

    return true;
  }
  // ***************************************************************************
  bool RefinementUtilities::SplitPolygon_IsAreaPositive(const Eigen::VectorXi& newCell2D_Indices,
                                                        const Eigen::Matrix3d& cell2DRotation,
                                                        const Eigen::Vector3d& cell2DTranslation,
                                                        IMeshDAO& mesh) const
  {
    double area_2D = 0.0;

    {
      Eigen::MatrixXd polygon_3D_vertices(3, newCell2D_Indices.size());
      for (unsigned int v = 0; v < newCell2D_Indices.size(); v++)
        polygon_3D_vertices.col(v)<< mesh.Cell0DCoordinates(newCell2D_Indices[v]);

      const Eigen::MatrixXd newCell2D_Vertices = geometryUtilities.RotatePointsFrom3DTo2D(polygon_3D_vertices,
                                                                                          cell2DRotation.transpose(),
                                                                                          cell2DTranslation);
      const std::vector<unsigned int> unaligned = geometryUtilities.UnalignedPoints(newCell2D_Vertices);
      const Eigen::MatrixXd unalignedVertices = geometryUtilities.ExtractPoints(newCell2D_Vertices,
                                                                                unaligned);
      const std::vector<unsigned int> triangles = geometryUtilities.PolygonTriangulationByFirstVertex(unalignedVertices);

      const std::vector<Eigen::Matrix3d> trianglesVertices = geometryUtilities.ExtractTriangulationPoints(unalignedVertices,
                                                                                                          triangles);

      for (unsigned int cct = 0; cct < trianglesVertices.size(); cct++)
        area_2D += geometryUtilities.PolygonArea(trianglesVertices[cct]);
    }

    return !geometryUtilities.IsValueZero(area_2D, geometryUtilities.Tolerance2D());
  }
  // ***************************************************************************
  RefinementUtilities::SplitPolygon_Result RefinementUtilities::SplitPolygon_NoNewVertices(const unsigned int cell2DIndex,
                                                                                           const unsigned int cell2DNumVertices,
                                                                                           const unsigned int fromVertex,
                                                                                           const unsigned int toVertex,
                                                                                           const Eigen::Matrix3d& cell2DRotation,
                                                                                           const Eigen::Vector3d& cell2DTranslation,
                                                                                           IMeshDAO& mesh) const
  {
    SplitPolygon_Result result;

    // Create new polygons list
    std::list<unsigned int> firstPolygonIndices;
    std::list<unsigned int> secondPolygonIndices;

    for (unsigned int v = 0; ((fromVertex + v) % cell2DNumVertices) != toVertex; v++)
      firstPolygonIndices.push_back((fromVertex + v) % cell2DNumVertices);
    for (unsigned int v = 0; ((toVertex + v) % cell2DNumVertices) != fromVertex; v++)
      secondPolygonIndices.push_back((toVertex + v) % cell2DNumVertices);

    // Create sub-cells vertices
    std::vector<Eigen::MatrixXi> subCells(2);
    subCells[0].resize(2, firstPolygonIndices.size() + 1);
    subCells[1].resize(2, secondPolygonIndices.size() + 1);

    {
      unsigned int v = 0;
      for (const unsigned int& index : firstPolygonIndices)
      {
        subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
        v++;
      }
      subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, toVertex);

      // Check split areas
      if (!SplitPolygon_IsAreaPositive(subCells[0].row(0).transpose(),
                                       cell2DRotation,
                                       cell2DTranslation,
                                       mesh))
      {
        result.Type = SplitPolygon_Result::Types::NoSplit;
        return result;
      }

      v = 0;
      for (const unsigned int& index : secondPolygonIndices)
      {
        subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
        v++;
      }
      subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, fromVertex);


      if (!SplitPolygon_IsAreaPositive(subCells[1].row(0).transpose(),
                                       cell2DRotation,
                                       cell2DTranslation,
                                       mesh))
      {
        result.Type = SplitPolygon_Result::Types::NoSplit;
        return result;
      }
    }

    // Create new cell1D from vertex to vertex
    {
      // check new edge length
      if (geometryUtilities.IsValueZero(geometryUtilities.SegmentLength(mesh.Cell2DVertexCoordinates(cell2DIndex, toVertex),
                                                                        mesh.Cell2DVertexCoordinates(cell2DIndex, fromVertex)),
                                        geometryUtilities.Tolerance1D()))
      {
        result.Type = SplitPolygon_Result::Types::NoSplit;
        return result;
      }
    }

    result.NewCell1DIndex = mesh.Cell1DAppend(1);
    const unsigned int& newCell1DIndex = result.NewCell1DIndex;
    mesh.Cell1DInsertExtremes(newCell1DIndex,
                              mesh.Cell2DVertex(cell2DIndex, toVertex),
                              mesh.Cell2DVertex(cell2DIndex, fromVertex));
    mesh.Cell1DSetMarker(newCell1DIndex, 0);
    mesh.Cell1DSetState(newCell1DIndex, true);
    mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex, 2);

    // Create sub-cells edges
    {
      unsigned int v = 0;
      for (const unsigned int& index : firstPolygonIndices)
      {
        subCells[0](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
        v++;
      }
      subCells[0](1, v) = newCell1DIndex;

      v = 0;
      for (const unsigned int& index : secondPolygonIndices)
      {
        subCells[1](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
        v++;
      }
      subCells[1](1, v) = newCell1DIndex;
    }

    result.NewCell2DsIndex = meshUtilities.SplitCell2D(cell2DIndex,
                                                       subCells,
                                                       mesh);

    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     0,
                                     result.NewCell2DsIndex[1]); // right
    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     1,
                                     result.NewCell2DsIndex[0]); // left

    result.Type = SplitPolygon_Result::Types::Split;

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::SplitPolygon_Result RefinementUtilities::SplitPolygon_NewVertexFrom(const unsigned int cell2DIndex,
                                                                                           const unsigned int cell2DNumVertices,
                                                                                           const unsigned int fromEdge,
                                                                                           const unsigned int toVertex,
                                                                                           const Eigen::Matrix3d& cell2DRotation,
                                                                                           const Eigen::Vector3d& cell2DTranslation,
                                                                                           const unsigned int fromNewCell0DIndex,
                                                                                           const std::vector<unsigned int>& fromSplitCell1DsIndex,
                                                                                           const bool& fromEdgeDirection,
                                                                                           IMeshDAO& mesh) const
  {
    SplitPolygon_Result result;

    // Create new polygons list
    std::list<unsigned int> firstPolygonIndices;
    std::list<unsigned int> secondPolygonIndices;

    for (unsigned int v = 0; ((fromEdge + 1 + v) % cell2DNumVertices) != toVertex; v++)
      firstPolygonIndices.push_back((fromEdge + 1 + v) % cell2DNumVertices);
    for (unsigned int v = 0; ((toVertex + v) % cell2DNumVertices) != fromEdge; v++)
      secondPolygonIndices.push_back((toVertex + v) % cell2DNumVertices);

    // Create sub-cells vertices
    std::vector<Eigen::MatrixXi> subCells(2);
    subCells[0].resize(2, firstPolygonIndices.size() + 2);
    subCells[1].resize(2, secondPolygonIndices.size() + 2);

    {
      unsigned int v = 0;
      subCells[0](0, v) = fromNewCell0DIndex;
      v++;
      for (const unsigned int& index : firstPolygonIndices)
      {
        subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
        v++;
      }
      subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, toVertex);

      // Check split areas
      if (!SplitPolygon_IsAreaPositive(subCells[0].row(0).transpose(),
                                       cell2DRotation,
                                       cell2DTranslation,
                                       mesh))
      {
        result.Type = SplitPolygon_Result::Types::NoSplit;
        return result;
      }

      v = 0;
      for (const unsigned int& index : secondPolygonIndices)
      {
        subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
        v++;
      }
      subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, fromEdge);

      v++;
      subCells[1](0, v) = fromNewCell0DIndex;
    }

    if (!SplitPolygon_IsAreaPositive(subCells[1].row(0).transpose(),
                                     cell2DRotation,
                                     cell2DTranslation,
                                     mesh))
    {
      result.Type = SplitPolygon_Result::Types::NoSplit;
      return result;
    }

    // Create new cell1D from new vertex to vertex
    {
      // check new edge length
      if (geometryUtilities.IsValueZero(geometryUtilities.SegmentLength(mesh.Cell2DVertexCoordinates(cell2DIndex, toVertex),
                                                                        mesh.Cell0DCoordinates(fromNewCell0DIndex)),
                                        geometryUtilities.Tolerance1D()))
      {
        result.Type = SplitPolygon_Result::Types::NoSplit;
        return result;
      }
    }

    result.NewCell1DIndex = mesh.Cell1DAppend(1);
    const unsigned int& newCell1DIndex = result.NewCell1DIndex;
    mesh.Cell1DInsertExtremes(newCell1DIndex,
                              mesh.Cell2DVertex(cell2DIndex, toVertex),
                              fromNewCell0DIndex);
    mesh.Cell1DSetMarker(newCell1DIndex, 0);
    mesh.Cell1DSetState(newCell1DIndex, true);
    mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex, 2);

    // Create sub-cells edges
    {
      unsigned int v = 0;
      subCells[0](1, v) = fromEdgeDirection ? fromSplitCell1DsIndex[1] : fromSplitCell1DsIndex[0];
      v++;
      for (const unsigned int& index : firstPolygonIndices)
      {
        subCells[0](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
        v++;
      }
      subCells[0](1, v) = newCell1DIndex;

      v = 0;
      for (const unsigned int& index : secondPolygonIndices)
      {
        subCells[1](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
        v++;
      }
      subCells[1](1, v) = fromEdgeDirection ? fromSplitCell1DsIndex[0] : fromSplitCell1DsIndex[1];
      v++;

      subCells[1](1, v) = newCell1DIndex;
    }

    result.NewCell2DsIndex = meshUtilities.SplitCell2D(cell2DIndex,
                                                       subCells,
                                                       mesh);

    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     0,
                                     result.NewCell2DsIndex[1]); // right
    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     1,
                                     result.NewCell2DsIndex[0]); // left

    result.Type = SplitPolygon_Result::Types::Split;

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::SplitPolygon_Result RefinementUtilities::SplitPolygon_NewVertexTo(const unsigned int cell2DIndex,
                                                                                         const unsigned int cell2DNumVertices,
                                                                                         const unsigned int fromVertex,
                                                                                         const unsigned int toEdge,
                                                                                         const Eigen::Matrix3d& cell2DRotation,
                                                                                         const Eigen::Vector3d& cell2DTranslation,
                                                                                         const unsigned int toNewCell0DIndex,
                                                                                         const std::vector<unsigned int>& toSplitCell1DsIndex,
                                                                                         const bool& toEdgeDirection,
                                                                                         IMeshDAO& mesh) const
  {
    SplitPolygon_Result result;

    // Create new polygons list
    std::list<unsigned int> firstPolygonIndices;
    std::list<unsigned int> secondPolygonIndices;

    for (unsigned int v = 0; ((fromVertex + v) % cell2DNumVertices) != toEdge; v++)
      firstPolygonIndices.push_back((fromVertex + v) % cell2DNumVertices);
    for (unsigned int v = 0; ((toEdge + 1 + v) % cell2DNumVertices) != fromVertex; v++)
      secondPolygonIndices.push_back((toEdge + 1 + v) % cell2DNumVertices);

    // Create sub-cells vertices
    std::vector<Eigen::MatrixXi> subCells(2);
    subCells[0].resize(2, firstPolygonIndices.size() + 2);
    subCells[1].resize(2, secondPolygonIndices.size() + 2);

    {
      unsigned int v = 0;
      for (const unsigned int& index : firstPolygonIndices)
      {
        subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
        v++;
      }
      subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, toEdge);
      v++;
      subCells[0](0, v) = toNewCell0DIndex;

      if (!SplitPolygon_IsAreaPositive(subCells[0].row(0).transpose(),
                                       cell2DRotation,
                                       cell2DTranslation,
                                       mesh))
      {
        result.Type = SplitPolygon_Result::Types::NoSplit;
        return result;
      }

      v = 0;
      subCells[1](0, v) = toNewCell0DIndex;
      v++;
      for (const unsigned int& index : secondPolygonIndices)
      {
        subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
        v++;
      }
      subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, fromVertex);
    }

    if (!SplitPolygon_IsAreaPositive(subCells[1].row(0).transpose(),
                                     cell2DRotation,
                                     cell2DTranslation,
                                     mesh))
    {
      result.Type = SplitPolygon_Result::Types::NoSplit;
      return result;
    }

    // Create new cell1D from vertex to new vertex
    {
      // check new edge length
      if (geometryUtilities.IsValueZero(geometryUtilities.SegmentLength(mesh.Cell0DCoordinates(toNewCell0DIndex),
                                                                        mesh.Cell2DVertexCoordinates(cell2DIndex, fromVertex)),
                                        geometryUtilities.Tolerance1D()))
      {
        result.Type = SplitPolygon_Result::Types::NoSplit;
        return result;
      }
    }

    result.NewCell1DIndex = mesh.Cell1DAppend(1);
    const unsigned int& newCell1DIndex = result.NewCell1DIndex;
    mesh.Cell1DInsertExtremes(newCell1DIndex,
                              toNewCell0DIndex,
                              mesh.Cell2DVertex(cell2DIndex, fromVertex));
    mesh.Cell1DSetMarker(newCell1DIndex, 0);
    mesh.Cell1DSetState(newCell1DIndex, true);
    mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex, 2);

    // Create sub-cells edges
    {
      unsigned int v = 0;
      for (const unsigned int& index : firstPolygonIndices)
      {
        subCells[0](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
        v++;
      }
      subCells[0](1, v) = toEdgeDirection ? toSplitCell1DsIndex[0] : toSplitCell1DsIndex[1];
      v++;
      subCells[0](1, v) = newCell1DIndex;

      v = 0;
      subCells[1](1, v) = toEdgeDirection ? toSplitCell1DsIndex[1] : toSplitCell1DsIndex[0];
      v++;
      for (const unsigned int& index : secondPolygonIndices)
      {
        subCells[1](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
        v++;
      }
      subCells[1](1, v) = newCell1DIndex;
    }

    result.NewCell2DsIndex = meshUtilities.SplitCell2D(cell2DIndex,
                                                       subCells,
                                                       mesh);

    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     0,
                                     result.NewCell2DsIndex[1]); // right
    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     1,
                                     result.NewCell2DsIndex[0]); // left

    result.Type = SplitPolygon_Result::Types::Split;

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::SplitPolygon_Result RefinementUtilities::SplitPolygon_NewVertices(const unsigned int cell2DIndex,
                                                                                         const unsigned int cell2DNumVertices,
                                                                                         const unsigned int fromEdge,
                                                                                         const unsigned int toEdge,
                                                                                         const Eigen::Matrix3d& cell2DRotation,
                                                                                         const Eigen::Vector3d& cell2DTranslation,
                                                                                         const unsigned int fromNewCell0DIndex,
                                                                                         const unsigned int toNewCell0DIndex,
                                                                                         const std::vector<unsigned int>& fromSplitCell1DsIndex,
                                                                                         const std::vector<unsigned int>& toSplitCell1DsIndex,
                                                                                         const bool& fromEdgeDirection,
                                                                                         const bool& toEdgeDirection,
                                                                                         IMeshDAO& mesh) const
  {
    SplitPolygon_Result result;

    // Create new polygons list
    std::list<unsigned int> firstPolygonIndices;
    std::list<unsigned int> secondPolygonIndices;

    for (unsigned int v = 0; ((fromEdge + 1 + v) % cell2DNumVertices) != toEdge; v++)
      firstPolygonIndices.push_back((fromEdge + 1 + v) % cell2DNumVertices);
    for (unsigned int v = 0; ((toEdge + 1 + v) % cell2DNumVertices) != fromEdge; v++)
      secondPolygonIndices.push_back((toEdge + 1 + v) % cell2DNumVertices);

    // Create sub-cells vertices
    std::vector<Eigen::MatrixXi> subCells(2);
    subCells[0].resize(2, firstPolygonIndices.size() + 3);
    subCells[1].resize(2, secondPolygonIndices.size() + 3);

    {
      unsigned int v = 0;
      subCells[0](0, v) = fromNewCell0DIndex;
      v++;
      for (const unsigned int& index : firstPolygonIndices)
      {
        subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
        v++;
      }
      subCells[0](0, v) = mesh.Cell2DVertex(cell2DIndex, toEdge);
      v++;
      subCells[0](0, v) = toNewCell0DIndex;

      // Check split areas
      if (!SplitPolygon_IsAreaPositive(subCells[0].row(0).transpose(),
                                       cell2DRotation,
                                       cell2DTranslation,
                                       mesh))
      {
        result.Type = SplitPolygon_Result::Types::NoSplit;
        return result;
      }

      v = 0;
      subCells[1](0, v) = toNewCell0DIndex;
      v++;
      for (const unsigned int& index : secondPolygonIndices)
      {
        subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, index);
        v++;
      }
      subCells[1](0, v) = mesh.Cell2DVertex(cell2DIndex, fromEdge);
      v++;
      subCells[1](0, v) = fromNewCell0DIndex;
    }

    if (!SplitPolygon_IsAreaPositive(subCells[1].row(0).transpose(),
                                     cell2DRotation,
                                     cell2DTranslation,
                                     mesh))
    {
      result.Type = SplitPolygon_Result::Types::NoSplit;
      return result;
    }

    // Create new cell1D from new vertex to new vertex
    {
      // check new edge length
      if (geometryUtilities.IsValueZero(geometryUtilities.SegmentLength(mesh.Cell0DCoordinates(toNewCell0DIndex),
                                                                        mesh.Cell0DCoordinates(fromNewCell0DIndex)),
                                        geometryUtilities.Tolerance1D()))
      {
        result.Type = SplitPolygon_Result::Types::NoSplit;
        return result;
      }
    }

    result.NewCell1DIndex = mesh.Cell1DAppend(1);
    const unsigned int& newCell1DIndex = result.NewCell1DIndex;
    mesh.Cell1DInsertExtremes(newCell1DIndex,
                              toNewCell0DIndex,
                              fromNewCell0DIndex);
    mesh.Cell1DSetMarker(newCell1DIndex, 0);
    mesh.Cell1DSetState(newCell1DIndex, true);
    mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex, 2);

    // Create sub-cells edges
    {
      unsigned int v = 0;
      subCells[0](1, v) = fromEdgeDirection ? fromSplitCell1DsIndex[1] : fromSplitCell1DsIndex[0];
      v++;
      for (const unsigned int& index : firstPolygonIndices)
      {
        subCells[0](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
        v++;
      }
      subCells[0](1, v) = toEdgeDirection ? toSplitCell1DsIndex[0] : toSplitCell1DsIndex[1];
      v++;
      subCells[0](1, v) = newCell1DIndex;

      v = 0;
      subCells[1](1, v) = toEdgeDirection ? toSplitCell1DsIndex[1] : toSplitCell1DsIndex[0];
      v++;
      for (const unsigned int& index : secondPolygonIndices)
      {
        subCells[1](1, v) = mesh.Cell2DEdge(cell2DIndex, index);
        v++;
      }
      subCells[1](1, v) = fromEdgeDirection ? fromSplitCell1DsIndex[0] : fromSplitCell1DsIndex[1];
      v++;
      subCells[1](1, v) = newCell1DIndex;
    }

    result.NewCell2DsIndex = meshUtilities.SplitCell2D(cell2DIndex,
                                                       subCells,
                                                       mesh);

    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     0,
                                     result.NewCell2DsIndex[1]); // right
    mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                     1,
                                     result.NewCell2DsIndex[0]); // left

    result.Type = SplitPolygon_Result::Types::Split;

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::TriangleMaxEdgeDirection RefinementUtilities::ComputeTriangleMaxEdgeDirection(const Eigen::VectorXd& edgesLength) const
  {
    TriangleMaxEdgeDirection result;

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
    std::iota(std::begin(cell2DsIndex),
              std::end(cell2DsIndex), 0);

    geometricData.Cell1Ds.Aligned.resize(mesh.Cell1DTotalNumber());
    std::iota(std::begin(geometricData.Cell1Ds.Aligned),
              std::end(geometricData.Cell1Ds.Aligned),
              1);
    geometricData.Cell1Ds.MaxAligned = mesh.Cell1DTotalNumber();

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
    geometricData.Cell1Ds.Aligned.resize(mesh.Cell1DTotalNumber(), 0);

    geometricData.Cell2Ds.Vertices.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.Area.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.Centroid.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.EdgesDirection.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.EdgesNormal.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.EdgesLength.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.Triangulations.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.Inertia.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.UnalignedVerticesIndex.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.UnalignedVertices.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.UnalignedEdgesLength.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.InRadius.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.CentroidEdgesDistance.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.CentroidVerticesDistance.resize(mesh.Cell2DTotalNumber());
    geometricData.Cell2Ds.Quality.resize(mesh.Cell2DTotalNumber());

    for (unsigned int c = 0; c < cell2DsIndex.size(); c++)
    {
      const unsigned int cell2DIndex = cell2DsIndex.at(c);

      const unsigned int numCell2DEdges = mesh.Cell2DNumberEdges(cell2DIndex);
      Eigen::MatrixXd& convexCell2DVertices = geometricData.Cell2Ds.Vertices.at(cell2DIndex);
      double& convexCell2DArea = geometricData.Cell2Ds.Area.at(cell2DIndex);
      std::vector<Eigen::Matrix3d>& convexCell2DTriangulationPoints = geometricData.Cell2Ds.Triangulations.at(cell2DIndex);
      Eigen::Vector3d& convexCell2DCentroid = geometricData.Cell2Ds.Centroid.at(cell2DIndex);
      Eigen::MatrixXd& convexCell2DUnalignedVertices = geometricData.Cell2Ds.UnalignedVertices.at(cell2DIndex);
      std::vector<unsigned int>& convexCell2DUnalignedVerticesFilter = geometricData.Cell2Ds.UnalignedVerticesIndex.at(cell2DIndex);

      convexCell2DVertices = mesh.Cell2DVerticesCoordinates(cell2DIndex);

      // compute original cell2D triangulation
      convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(convexCell2DVertices);
      convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(convexCell2DVertices,
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

      geometricData.Cell2Ds.EdgesDirection[cell2DIndex].resize(numCell2DEdges);
      for (unsigned int e = 0; e < numCell2DEdges; e++)
      {
        const unsigned int origin = mesh.Cell2DVertex(cell2DIndex, e);
        const unsigned int end = mesh.Cell2DVertex(cell2DIndex,
                                                   (e + 1) % numCell2DEdges);

        geometricData.Cell2Ds.EdgesDirection[cell2DIndex][e] = (mesh.Cell2DFindEdgeByExtremes(cell2DIndex,
                                                                                              origin,
                                                                                              end) == e);
      }

      geometricData.Cell2Ds.UnalignedEdgesLength[cell2DIndex] = geometryUtilities.PolygonEdgeLengths(convexCell2DUnalignedVertices);
      geometricData.Cell2Ds.EdgesLength[cell2DIndex] = geometryUtilities.PolygonEdgeLengths(convexCell2DVertices);

      geometricData.Cell2Ds.EdgesNormal[cell2DIndex] = geometryUtilities.PolygonEdgeNormals(convexCell2DVertices);
      geometricData.Cell2Ds.CentroidEdgesDistance[cell2DIndex] = geometryUtilities.PolygonCentroidEdgesDistance(convexCell2DVertices,
                                                                                                                convexCell2DCentroid,
                                                                                                                geometricData.Cell2Ds.EdgesNormal[cell2DIndex]);
      geometricData.Cell2Ds.CentroidVerticesDistance[cell2DIndex] = geometryUtilities.PolygonCentroidVerticesDistance(convexCell2DVertices,
                                                                                                                      convexCell2DCentroid);
      geometricData.Cell2Ds.InRadius[cell2DIndex] = geometryUtilities.PolygonInRadius(geometricData.Cell2Ds.CentroidEdgesDistance[cell2DIndex]);
      geometricData.Cell2Ds.Inertia[cell2DIndex] = geometryUtilities.PolygonInertia(convexCell2DCentroid,
                                                                                    convexCell2DTriangulationPoints);

      geometricData.Cell2Ds.Quality[cell2DIndex] = std::min(geometricData.Cell2Ds.EdgesLength[cell2DIndex].minCoeff(),
                                                            geometricData.Cell2Ds.InRadius[cell2DIndex]);

      for (unsigned int e = 0; e < numCell2DEdges; e++)
      {
        const unsigned int cell1DIndex = mesh.Cell2DEdge(cell2DIndex, e);

        if (mesh.Cell1DHasOriginalCell1D(cell1DIndex))
        {
          const unsigned int originalCell1DIndex = mesh.Cell1DOriginalCell1D(cell1DIndex);
          geometricData.Cell1Ds.Aligned[cell1DIndex] = geometricData.Cell1Ds.Aligned[originalCell1DIndex];
        }
        else if (geometricData.Cell1Ds.Aligned[cell1DIndex] == 0)
          geometricData.Cell1Ds.Aligned[cell1DIndex] = ++(geometricData.Cell1Ds.MaxAligned);
      }
    }
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolygon_CheckResult::Cell1DToSplit RefinementUtilities::RefinePolygonCell_IsCell1DToSplit(const unsigned int& cell1DIndex,
                                                                                                                       const unsigned int& cell2DIndex,
                                                                                                                       const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersection,
                                                                                                                       const std::vector<Eigen::VectorXd>& cell2DsEdgesLength,
                                                                                                                       const double& cell1DsQualityWeight,
                                                                                                                       const double& cell1DsAlignedWeight,
                                                                                                                       const std::vector<double>& cell2DsQuality,
                                                                                                                       const std::vector<unsigned int>& cell1DsAligned,
                                                                                                                       const IMeshDAO& mesh) const
  {
    RefinePolygon_CheckResult::Cell1DToSplit result;

    result.Cell2DEdgeIndex = edgeIntersection.Index;

    // check if the edge intersection is inside edge
    result.IsIntersectionInside = (edgeIntersection.Type ==
                                   GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);

    if (!result.IsIntersectionInside)
    {
      result.Type = RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside;
      return result;
    }

    const double cell1DLength = cell2DsEdgesLength[cell2DIndex][edgeIntersection.Index];

    // check if the edge length is enough
    result.IsEdgeLengthEnough = !geometryUtilities.IsValueZero(0.5 * cell1DLength,
                                                               geometryUtilities.Tolerance1D());

    if (!result.IsEdgeLengthEnough)
    {
      result.Type = RefinePolygon_CheckResult::Cell1DToSplit::Types::EdgeLengthNotEnough;
      return result;
    }

    // check if the edge quality is enough
    result.IsLocalQualityEnough = true;
    result.IsQualityEnough = true;
    result.IsNeighQualityEnough.resize(mesh.Cell1DNumberNeighbourCell2D(cell1DIndex),
                                       true);

    for (unsigned int c2Dn = 0; c2Dn < mesh.Cell1DNumberNeighbourCell2D(cell1DIndex); c2Dn++)
    {
      if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, c2Dn))
        continue;

      const unsigned int cell2DNeighIndex = mesh.Cell1DNeighbourCell2D(cell1DIndex, c2Dn);

      result.IsNeighQualityEnough[c2Dn] = geometryUtilities.IsValueGreaterOrEqual(0.5 * cell1DLength,
                                                                                  cell1DsQualityWeight * cell2DsQuality[cell2DNeighIndex],
                                                                                  geometryUtilities.Tolerance1D());

      if (cell2DNeighIndex == cell2DIndex &&
          !result.IsNeighQualityEnough[c2Dn])
      {
        if (!result.IsQualityEnough)
          result.Type = RefinePolygon_CheckResult::Cell1DToSplit::Types::BothQualityNotEnough;
        else
          result.Type = RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyLocalQualityNotEnough;

        result.IsLocalQualityEnough = false;
      }

      if (!result.IsNeighQualityEnough[c2Dn])
      {
        if (!result.IsLocalQualityEnough)
          result.Type = RefinePolygon_CheckResult::Cell1DToSplit::Types::BothQualityNotEnough;
        else
          result.Type = RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyNeighQualityNotEnough;

        result.IsQualityEnough = false;
      }
    }

    if (!result.IsQualityEnough)
      return result;

    // check if edge is part of aligned edges in cell2D or neighbours
    const unsigned int aligned = cell1DsAligned.at(cell1DIndex);

    result.IsLocalAlignedRespect = true;
    result.IsAlignedRespect = true;
    result.IsNeighAlignedRespect.resize(mesh.Cell1DNumberNeighbourCell2D(cell1DIndex),
                                        true);

    for (unsigned int c2Dn = 0; c2Dn < mesh.Cell1DNumberNeighbourCell2D(cell1DIndex); c2Dn++)
    {
      if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, c2Dn))
        continue;

      const unsigned int cell2DNeighIndex = mesh.Cell1DNeighbourCell2D(cell1DIndex, c2Dn);

      unsigned int numAligned = 0;
      double alignedEdgesLength = 0.0;
      for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(cell2DNeighIndex); ed++)
      {
        const unsigned int cell1DEdgeIndex = mesh.Cell2DEdge(cell2DNeighIndex, ed);

        if (cell1DsAligned.at(cell1DEdgeIndex) == aligned)
        {
          numAligned++;
          alignedEdgesLength += cell2DsEdgesLength[cell2DNeighIndex][ed];
        }
      }

      // check if the length of the aligned splitted edge (|l|) is higher than the mean length (|L|) after the cut
      // |l| / 2 > c_alg * |L| / (N_aligned + 1)
      result.IsNeighAlignedRespect[c2Dn] = geometryUtilities.IsValueGreaterOrEqual(0.5 * cell1DLength,
                                                                                   cell1DsAlignedWeight * alignedEdgesLength / double(numAligned + 1),
                                                                                   geometryUtilities.Tolerance1D());


      if (cell2DNeighIndex == cell2DIndex &&
          !result.IsNeighAlignedRespect[c2Dn])
      {
        if (!result.IsAlignedRespect)
          result.Type = RefinePolygon_CheckResult::Cell1DToSplit::Types::BothAlignedNotRespect;
        else
          result.Type = RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyLocalAlignedNotRespect;

        result.IsLocalAlignedRespect = false;
      }

      if (!result.IsNeighAlignedRespect[c2Dn])
      {
        if (!result.IsLocalAlignedRespect)
          result.Type = RefinePolygon_CheckResult::Cell1DToSplit::Types::BothAlignedNotRespect;
        else
          result.Type = RefinePolygon_CheckResult::Cell1DToSplit::Types::OnlyNeighAlignedNotRespect;

        result.IsAlignedRespect = false;
      }
    }

    if (!result.IsAlignedRespect)
      return result;

    result.Type = RefinePolygon_CheckResult::Cell1DToSplit::Types::ToSplit;
    result.IsToSplit = true;
    return result;
  }
  // ***************************************************************************
  RefinementUtilities::PolygonDirection RefinementUtilities::ComputePolygonMaxDiameterDirection(const Eigen::MatrixXd unalignedVertices,
                                                                                                const Eigen::Vector3d& centroid) const
  {
    PolygonDirection result;

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
  RefinementUtilities::PolygonDirection RefinementUtilities::ComputePolygonMaxInertiaDirection(const Eigen::MatrixXd& unalignedVertices,
                                                                                               const Eigen::VectorXd& unalignedEdgesLength,
                                                                                               const Eigen::Vector3d& centroid,
                                                                                               const Eigen::Matrix3d& inertia) const
  {
    PolygonDirection result;

    if (unalignedVertices.cols() == 3)
    {
      Eigen::VectorXd::Index maxEdgeLocalIndex;
      unalignedEdgesLength.maxCoeff(&maxEdgeLocalIndex);
      const unsigned int oppositeVertexIndex = (maxEdgeLocalIndex + 2) % 3;

      result.LineTangent = geometryUtilities.SegmentTangent(unalignedVertices.col(oppositeVertexIndex),
                                                            0.5 * (unalignedVertices.col(maxEdgeLocalIndex) +
                                                                   unalignedVertices.col((maxEdgeLocalIndex + 1) % 3))).normalized();
      result.LineOrigin = unalignedVertices.col(oppositeVertexIndex);
    }
    else
    {
      // get maximux direction of inertia tensor
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver;

      Eigen::Matrix2d eig_inertia = inertia.block(0, 0, 2, 2);
      if (eig_inertia.isDiagonal())
      {
        eig_inertia(0, 1) = 0.0;
        eig_inertia(1, 0) = 0.0;
      }

      eigensolver.computeDirect(eig_inertia);
      if (eigensolver.info() != Eigen::Success)
        throw std::runtime_error("Inertia not correct");

      result.LineTangent.setZero();
      result.LineTangent.segment(0, 2) = eigensolver.eigenvectors().col(1);
      result.LineOrigin = centroid;
    }

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolygon_Result RefinementUtilities::RefineTriangleCell_ByEdge(const unsigned int& cell2DIndex,
                                                                                           const unsigned int& edgeIndex,
                                                                                           const unsigned int& oppositeVertexIndex,
                                                                                           const std::vector<bool>& cell2DEdgesDirection,
                                                                                           const double& cell2DArea,
                                                                                           const Eigen::Matrix3d& cell2DRotation,
                                                                                           const Eigen::Vector3d& cell2DTranslation,
                                                                                           const Eigen::VectorXd& cell2DEdgesLength,
                                                                                           IMeshDAO& mesh) const
  {
    RefinePolygon_Result result;

    if (mesh.Cell2DHasUpdatedCell2Ds(cell2DIndex))
    {
      result.ResultType = RefinePolygon_Result::ResultTypes::Cell2DAlreadySplitted;
      return result;
    }

    // Check if the cell refinement creates null cell2Ds or cell1Ds
    if (geometryUtilities.IsValueZero(0.5 * cell2DArea,
                                      geometryUtilities.Tolerance2D()) ||
        geometryUtilities.IsValueZero(0.5 * cell2DEdgesLength[edgeIndex],
                                      geometryUtilities.Tolerance1D()))
    {
      result.ResultType = RefinePolygon_Result::ResultTypes::Cell2DSplitUnderTolerance;
      return result;
    }

    // Add new mesh vertex, middle point of edge
    const unsigned int cell1DIndex = mesh.Cell2DEdge(cell2DIndex, edgeIndex);
    const bool cell1DDirection = cell2DEdgesDirection[edgeIndex];

    const SplitCell1D_Result splitCell1D = SplitCell1D_MiddlePoint(cell1DIndex,
                                                                   mesh);

    const SplitPolygon_Result splitResult = SplitPolygon_NewVertexTo(cell2DIndex,
                                                                     3,
                                                                     oppositeVertexIndex,
                                                                     edgeIndex,
                                                                     cell2DRotation,
                                                                     cell2DTranslation,
                                                                     splitCell1D.NewCell0DIndex,
                                                                     splitCell1D.NewCell1DsIndex,
                                                                     cell1DDirection,
                                                                     mesh);
    switch (splitResult.Type)
    {
      case SplitPolygon_Result::Types::Unknown:
        throw std::runtime_error("Unknown SplitPolygon_Result Type");
      case SplitPolygon_Result::Types::NoSplit:
        mesh.Cell1DRemove(splitCell1D.NewCell1DsIndex[1]);
        mesh.Cell1DRemove(splitCell1D.NewCell1DsIndex[0]);
        mesh.Cell0DRemove(splitCell1D.NewCell0DIndex);
        mesh.Cell1DSetState(cell1DIndex, true);
        result.ResultType = RefinePolygon_Result::ResultTypes::Cell2DSplitUnderTolerance;
        return result;
      default:
        break;
    }

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
                                                                const std::vector<Eigen::Matrix3d>& cell2DsRotation,
                                                                const std::vector<Eigen::Vector3d>& cell2DsTranslation,
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

      const Eigen::Matrix3d& cell2DRotation = cell2DsRotation.at(n);
      const Eigen::Vector3d& cell2DTranslation = cell2DsTranslation.at(n);

      SplitPolygon_NewVertexTo(neighCell2DIndex,
                               3,
                               neighOppositeVertexIndex,
                               neighEdgeIndex,
                               cell2DRotation,
                               cell2DTranslation,
                               newCell0DIndex,
                               splitCell1DsIndex,
                               !cell2DEdgeDirection,
                               mesh);
    }
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolygon_CheckResult RefinementUtilities::RefinePolygonCell_CheckRefinement(const unsigned int& cell2DIndex,
                                                                                                        const Eigen::MatrixXd& cell2DVertices,
                                                                                                        const Eigen::Vector3d& lineTangent,
                                                                                                        const Eigen::Vector3d& lineOrigin,
                                                                                                        const std::vector<double>& cell2DsQuality,
                                                                                                        const std::vector<unsigned int>& cell1DsAligned,
                                                                                                        const double& cell1DsQualityWeight,
                                                                                                        const double& cell1DsAlignedWeight,
                                                                                                        const double& cell2DArea,
                                                                                                        const std::vector<Eigen::VectorXd>& cell2DsEdgesLength,
                                                                                                        const std::vector<bool>& cell2DEdgesDirection,
                                                                                                        IMeshDAO& mesh) const
  {
    RefinePolygon_CheckResult result;

    if (mesh.Cell2DHasUpdatedCell2Ds(cell2DIndex))
    {
      result.ResultType = RefinePolygon_CheckResult::ResultTypes::Cell2DAlreadySplitted;
      return result;
    }

    // Check if the cell refinement creates null cell2Ds
    if (geometryUtilities.IsValueZero(0.5 * cell2DArea,
                                      geometryUtilities.Tolerance2D()))
    {
      result.ResultType = RefinePolygon_CheckResult::ResultTypes::Cell2DSplitUnderTolerance;
      return result;
    }

    GeometryUtilities::LinePolygonPositionResult linePolygonPosition = geometryUtilities.LinePolygonPosition(lineTangent,
                                                                                                             lineOrigin,
                                                                                                             cell2DVertices);
    Output::Assert(linePolygonPosition.Type != GeometryUtilities::LinePolygonPositionResult::Types::Unknown);

    if (linePolygonPosition.Type != GeometryUtilities::LinePolygonPositionResult::Types::Intersecting)
    {
      result.ResultType = RefinePolygon_CheckResult::ResultTypes::SplitDirectionNotInsideCell2D;
      return result;
    }

    if (linePolygonPosition.EdgeIntersections.size() < 2)
      throw std::runtime_error("Not enough intersections found");

    const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersectionOne = linePolygonPosition.EdgeIntersections[0];
    Output::Assert(edgeIntersectionOne.Type != GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::Unknown);
    if (edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::Parallel)
    {
      result.ResultType = RefinePolygon_CheckResult::ResultTypes::SplitDirectionNotInsideCell2D;
      return result;
    }

    const unsigned int cell1DIndexOne = mesh.Cell2DEdge(cell2DIndex,
                                                        edgeIntersectionOne.Index);

    GeometryUtilities::LinePolygonPositionResult::EdgeIntersection* findEdgeIntersectionTwo = nullptr;

    for (unsigned int i = 1; i < linePolygonPosition.EdgeIntersections.size(); i++)
    {
      const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersectionTwo = linePolygonPosition.EdgeIntersections[i];

      switch (edgeIntersectionTwo.Type)
      {
        case GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::Parallel:
        {
          result.ResultType = RefinePolygon_CheckResult::ResultTypes::SplitDirectionNotInsideCell2D;
          return result;
        }
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

    const RefinePolygon_CheckResult::Cell1DToSplit createNewVertexOne = RefinePolygonCell_IsCell1DToSplit(cell1DIndexOne,
                                                                                                          cell2DIndex,
                                                                                                          edgeIntersectionOne,
                                                                                                          cell2DsEdgesLength,
                                                                                                          cell1DsQualityWeight,
                                                                                                          cell1DsAlignedWeight,
                                                                                                          cell2DsQuality,
                                                                                                          cell1DsAligned,
                                                                                                          mesh);
    const RefinePolygon_CheckResult::Cell1DToSplit createNewVertexTwo = RefinePolygonCell_IsCell1DToSplit(cell1DIndexTwo,
                                                                                                          cell2DIndex,
                                                                                                          edgeIntersectionTwo,
                                                                                                          cell2DsEdgesLength,
                                                                                                          cell1DsQualityWeight,
                                                                                                          cell1DsAlignedWeight,
                                                                                                          cell2DsQuality,
                                                                                                          cell1DsAligned,
                                                                                                          mesh);

    result.Cell1DsIndex = { cell1DIndexOne, cell1DIndexTwo };
    result.Cell1DsIntersection = { edgeIntersectionOne, edgeIntersectionTwo };
    result.Cell1DsToSplit = { createNewVertexOne, createNewVertexTwo };

    result.ResultType = RefinePolygon_CheckResult::ResultTypes::Cell2DToBeSplitted;

    return result;
  }
  // ***************************************************************************
  RefinementUtilities::RefinePolygon_Result RefinementUtilities::RefinePolygonCell_ByDirection(const unsigned int& cell2DIndex,
                                                                                               const Gedim::GeometryUtilities::PolygonTypes& cell2DPolygonType,
                                                                                               const Gedim::GeometryUtilities::PolygonTypes& cell2DUnalignedPolygonType,
                                                                                               const Eigen::MatrixXd& cell2DVertices,
                                                                                               const RefinePolygon_CheckResult& cell2DCheckToRefine,
                                                                                               const Eigen::Matrix3d& cell2DRotation,
                                                                                               const Eigen::Vector3d& cell2DTranslation,
                                                                                               const std::vector<bool>& cell2DEdgesDirection,
                                                                                               const bool& extendToNeighbours,
                                                                                               IMeshDAO& mesh) const
  {
    RefinePolygon_Result result;

    switch (cell2DCheckToRefine.ResultType)
    {
      case RefinePolygon_CheckResult::ResultTypes::Cell2DAlreadySplitted:
      {
        result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
        result.ResultType = RefinePolygon_Result::ResultTypes::Cell2DAlreadySplitted;
        return result;
      }
      case RefinePolygon_CheckResult::ResultTypes::Cell2DSplitUnderTolerance:
      {
        result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
        result.ResultType = RefinePolygon_Result::ResultTypes::Cell2DSplitUnderTolerance;
        return result;
      }
      case RefinePolygon_CheckResult::ResultTypes::SplitDirectionNotInsideCell2D:
      {
        result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
        result.ResultType = RefinePolygon_Result::ResultTypes::SplitDirectionNotInsideCell2D;
        return result;
      }
      case RefinePolygon_CheckResult::ResultTypes::Cell2DToBeSplitted:
        break;
      default:
        throw std::runtime_error("Unmanaged RefinePolygon_CheckResult type");
    }

    const unsigned int cell2DNumVertices = cell2DVertices.cols();
    Gedim::Output::Assert(cell2DCheckToRefine.Cell1DsIndex.size() == 2);
    const unsigned int cell1DIndexOne = cell2DCheckToRefine.Cell1DsIndex[0];
    const unsigned int cell1DIndexTwo = cell2DCheckToRefine.Cell1DsIndex[1];

    Gedim::Output::Assert(cell2DCheckToRefine.Cell1DsIntersection.size() == 2);
    const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersectionOne = cell2DCheckToRefine.Cell1DsIntersection[0];
    const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersectionTwo = cell2DCheckToRefine.Cell1DsIntersection[1];

    Gedim::Output::Assert(cell2DCheckToRefine.Cell1DsToSplit.size() == 2);
    const RefinePolygon_CheckResult::Cell1DToSplit& createNewVertexOne = cell2DCheckToRefine.Cell1DsToSplit[0];
    const RefinePolygon_CheckResult::Cell1DToSplit& createNewVertexTwo = cell2DCheckToRefine.Cell1DsToSplit[1];

    if (extendToNeighbours &&
        cell2DUnalignedPolygonType != Gedim::GeometryUtilities::PolygonTypes::Triangle &&
        !SplitPolygon_CheckIsNotToExtend(createNewVertexOne,
                                         createNewVertexTwo))
    {
      result.ResultType = RefinePolygon_Result::ResultTypes::SplitQualityCheckCell2DFailed;
      result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
      return result;
    }

    const CheckSplitType_Result splitTypeResult = SplitPolygon_CheckSplitType(cell2DPolygonType,
                                                                              cell2DUnalignedPolygonType,
                                                                              cell2DVertices,
                                                                              cell2DCheckToRefine);

    switch (splitTypeResult.Type)
    {
      case CheckSplitType_Result::SplitTypes::NoNewVertices:
      {
        // no new vertices
        const unsigned int fromVertex = splitTypeResult.NoNewVerticesIndex[0];
        const unsigned int toVertex = splitTypeResult.NoNewVerticesIndex[1];

        if (AreVerticesAligned(cell2DVertices, fromVertex, toVertex))
        {
          if (extendToNeighbours &&
              !SplitPolygon_CheckIsToSplit_Relaxed(createNewVertexOne,
                                                   createNewVertexTwo))
          {
            result.ResultType = RefinePolygon_Result::ResultTypes::SplitQualityCheckCell2DFailed;
            result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
            return result;
          }

          result.ResultType = RefinePolygon_Result::ResultTypes::SplitDirectionNotInsideCell2D;
          result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
          throw std::runtime_error("Case to manage");
          return result;
        }

        const SplitPolygon_Result splitResult = SplitPolygon_NoNewVertices(cell2DIndex,
                                                                           cell2DNumVertices,
                                                                           fromVertex,
                                                                           toVertex,
                                                                           cell2DRotation,
                                                                           cell2DTranslation,
                                                                           mesh);

        switch (splitResult.Type)
        {
          case SplitPolygon_Result::Types::Unknown:
            throw std::runtime_error("Unknown SplitPolygon_Result Type");
          case SplitPolygon_Result::Types::NoSplit:
          {
            result.ResultType = RefinePolygon_Result::ResultTypes::Cell2DSplitUnderTolerance;
            result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
            return result;
          }
          default:
            result.SplitType = CheckSplitType_Result::SplitTypes::NoNewVertices;
            break;
        }

        result.NewCell1DsIndex.resize(1);
        result.NewCell1DsIndex[0].Type = RefinePolygon_Result::RefinedCell1D::Types::New;
        result.NewCell1DsIndex[0].NewCell1DsIndex = { splitResult.NewCell1DIndex };

        result.NewCell2DsIndex = splitResult.NewCell2DsIndex;
      }
        break;
      case CheckSplitType_Result::SplitTypes::NewVertexFrom:
      {
        // check if edge can be splitted
        if (createNewVertexOne.Type ==
            RefinePolygon_CheckResult::Cell1DToSplit::Types::EdgeLengthNotEnough)
        {
          result.ResultType = RefinePolygon_Result::ResultTypes::Cell2DSplitUnderTolerance;
          result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
          return result;
        }

        // new vertex one
        Gedim::Output::Assert(edgeIntersectionOne.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);

        unsigned int toVertex = geometryUtilities.IsValueGreaterOrEqual(0.5, edgeIntersectionTwo.CurvilinearCoordinate,
                                                                        geometryUtilities.Tolerance1D()) ?
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
        {
          // try to flip vertices
          if (createNewVertexTwo.Type !=
              RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside)
          {
            // no fixed vertices, flip toVertex
            toVertex = (toVertex == edgeIntersectionTwo.Index) ?
                         (edgeIntersectionTwo.Index + 1) % cell2DNumVertices :
                         edgeIntersectionTwo.Index;
          }

          if (geometryUtilities.PointIsAligned(cell2DVertices.col(edgeIntersectionOne.Index),
                                               cell2DVertices.col(toVertex),
                                               middleCoordinate) ||
              geometryUtilities.PointIsAligned(cell2DVertices.col((edgeIntersectionOne.Index + 1) % cell2DNumVertices),
                                               cell2DVertices.col(toVertex),
                                               middleCoordinate))
          {
            result.ResultType = RefinePolygon_Result::ResultTypes::SplitDirectionNotInsideCell2D;
            result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
            throw std::runtime_error("Case to manage");
            return result;
          }
        }

        const SplitCell1D_Result splitCell1DOne = SplitCell1D_MiddlePoint(cell1DIndexOne,
                                                                          mesh);

        const SplitPolygon_Result splitResult = SplitPolygon_NewVertexFrom(cell2DIndex,
                                                                           cell2DNumVertices,
                                                                           edgeIntersectionOne.Index,
                                                                           toVertex,
                                                                           cell2DRotation,
                                                                           cell2DTranslation,
                                                                           splitCell1DOne.NewCell0DIndex,
                                                                           splitCell1DOne.NewCell1DsIndex,
                                                                           cell2DEdgesDirection.at(edgeIntersectionOne.Index),
                                                                           mesh);

        switch (splitResult.Type)
        {
          case SplitPolygon_Result::Types::Unknown:
            throw std::runtime_error("Unknown SplitPolygon_Result Type");
          case SplitPolygon_Result::Types::NoSplit:
          {
            mesh.Cell1DRemove(splitCell1DOne.NewCell1DsIndex[1]);
            mesh.Cell1DRemove(splitCell1DOne.NewCell1DsIndex[0]);
            mesh.Cell0DRemove(splitCell1DOne.NewCell0DIndex);
            mesh.Cell1DSetState(cell1DIndexOne, true);
            result.ResultType = RefinePolygon_Result::ResultTypes::Cell2DSplitUnderTolerance;
            result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
            return result;
          }
          default:
            result.SplitType = CheckSplitType_Result::SplitTypes::NewVertexFrom;
            break;
        }

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
        break;
      case CheckSplitType_Result::SplitTypes::NewVertexTo:
      {
        // check if edge can be divided
        if (createNewVertexTwo.Type ==
            RefinePolygon_CheckResult::Cell1DToSplit::Types::EdgeLengthNotEnough)
        {
          result.ResultType = RefinePolygon_Result::ResultTypes::Cell2DSplitUnderTolerance;
          result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
          return result;
        }

        // new vertex two
        Gedim::Output::Assert(edgeIntersectionTwo.Type == GeometryUtilities::LinePolygonPositionResult::EdgeIntersection::Types::InsideEdge);

        unsigned int fromVertex = geometryUtilities.IsValueGreaterOrEqual(0.5, edgeIntersectionOne.CurvilinearCoordinate,
                                                                          geometryUtilities.Tolerance1D()) ?
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
        {
          // try to flip vertices
          if (createNewVertexOne.Type !=
              RefinePolygon_CheckResult::Cell1DToSplit::Types::NotInside)
          {
            // no fixed vertices, flip toVertex
            fromVertex = (fromVertex == edgeIntersectionOne.Index) ?
                           (edgeIntersectionOne.Index + 1) % cell2DNumVertices :
                           edgeIntersectionOne.Index;
          }

          if (geometryUtilities.PointIsAligned(cell2DVertices.col(fromVertex),
                                               cell2DVertices.col(edgeIntersectionTwo.Index),
                                               middleCoordinate) ||
              geometryUtilities.PointIsAligned(cell2DVertices.col(fromVertex),
                                               cell2DVertices.col((edgeIntersectionTwo.Index + 1) % cell2DNumVertices),
                                               middleCoordinate))
          {
            result.ResultType = RefinePolygon_Result::ResultTypes::SplitDirectionNotInsideCell2D;
            result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
            throw std::runtime_error("Case to manage");
            return result;
          }
        }

        const SplitCell1D_Result splitCell1DTwo = SplitCell1D_MiddlePoint(cell1DIndexTwo,
                                                                          mesh);

        const SplitPolygon_Result splitResult = SplitPolygon_NewVertexTo(cell2DIndex,
                                                                         cell2DNumVertices,
                                                                         fromVertex,
                                                                         edgeIntersectionTwo.Index,
                                                                         cell2DRotation,
                                                                         cell2DTranslation,
                                                                         splitCell1DTwo.NewCell0DIndex,
                                                                         splitCell1DTwo.NewCell1DsIndex,
                                                                         cell2DEdgesDirection.at(edgeIntersectionTwo.Index),
                                                                         mesh);

        switch (splitResult.Type)
        {
          case SplitPolygon_Result::Types::Unknown:
            throw std::runtime_error("Unknown SplitPolygon_Result Type");
          case SplitPolygon_Result::Types::NoSplit:
          {
            mesh.Cell1DRemove(splitCell1DTwo.NewCell1DsIndex[1]);
            mesh.Cell1DRemove(splitCell1DTwo.NewCell1DsIndex[0]);
            mesh.Cell0DRemove(splitCell1DTwo.NewCell0DIndex);
            mesh.Cell1DSetState(cell1DIndexTwo, true);
            result.ResultType = RefinePolygon_Result::ResultTypes::Cell2DSplitUnderTolerance;
            result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
            return result;
          }
          default:
            result.SplitType = CheckSplitType_Result::SplitTypes::NewVertexTo;
            break;
        }

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
        break;
      case CheckSplitType_Result::SplitTypes::NewVertices:
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
                                                                         cell2DRotation,
                                                                         cell2DTranslation,
                                                                         splitCell1DOne.NewCell0DIndex,
                                                                         splitCell1DTwo.NewCell0DIndex,
                                                                         splitCell1DOne.NewCell1DsIndex,
                                                                         splitCell1DTwo.NewCell1DsIndex,
                                                                         cell2DEdgesDirection.at(edgeIntersectionOne.Index),
                                                                         cell2DEdgesDirection.at(edgeIntersectionTwo.Index),
                                                                         mesh);

        switch (splitResult.Type)
        {
          case SplitPolygon_Result::Types::Unknown:
            throw std::runtime_error("Unknown SplitPolygon_Result Type");
          case SplitPolygon_Result::Types::NoSplit:
          {
            mesh.Cell1DRemove(splitCell1DTwo.NewCell1DsIndex[1]);
            mesh.Cell1DRemove(splitCell1DTwo.NewCell1DsIndex[0]);
            mesh.Cell1DRemove(splitCell1DOne.NewCell1DsIndex[1]);
            mesh.Cell1DRemove(splitCell1DOne.NewCell1DsIndex[0]);
            mesh.Cell0DRemove(splitCell1DTwo.NewCell0DIndex);
            mesh.Cell0DRemove(splitCell1DOne.NewCell0DIndex);
            mesh.Cell1DSetState(cell1DIndexTwo, true);
            mesh.Cell1DSetState(cell1DIndexOne, true);
            result.ResultType = RefinePolygon_Result::ResultTypes::Cell2DSplitUnderTolerance;
            result.SplitType = CheckSplitType_Result::SplitTypes::NoSplit;
            return result;
          }
          default:
            result.SplitType = CheckSplitType_Result::SplitTypes::NewVertices;
            break;
        }

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
        break;
      default:
        throw std::runtime_error("Split type not managed case!");
    }

    result.ResultType = RefinePolygon_Result::ResultTypes::Successfull;
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
