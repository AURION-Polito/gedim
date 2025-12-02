// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#include "RefinementUtilities.hpp"
#include "VTKUtilities.hpp"

namespace Gedim
{
// ***************************************************************************
std::vector<unsigned int> RefinementUtilities::refine_mesh_2D_triangles(const Gedim::GeometryUtilities &geometry_utilities,
                                                                        const std::vector<unsigned int> &cell2DsToRefineIndex,
                                                                        Gedim::RefinementUtilities::Cell2Ds_GeometricData &meshGeometricData,
                                                                        Gedim::IMeshDAO &mesh) const
{
    struct Cell2DsToUpdateGeometricData final
    {
        unsigned int OriginalCell2DIndex = 0;
        unsigned int NewCell2DIndex = 0;
    };

    std::list<unsigned int> toRefineCells(cell2DsToRefineIndex.begin(), cell2DsToRefineIndex.end());
    std::list<unsigned int> refinedCells;

    while (!toRefineCells.empty())
    {
        const unsigned int cell2DIndex = toRefineCells.front();

        toRefineCells.pop_front();

        std::list<Cell2DsToUpdateGeometricData> cell2DsToUpdateGeometricData;

        std::list<unsigned int> updatedCell2Ds;
        mesh.Cell2DUpdatedCell2Ds(cell2DIndex, updatedCell2Ds);

        if (updatedCell2Ds.size() > 1)
            continue;

        const unsigned int cell2DToRefineIndex = updatedCell2Ds.size() == 0 ? cell2DIndex : updatedCell2Ds.front();

        if (meshGeometricData.Cell2Ds.Status[cell2DToRefineIndex] ==
            Gedim::RefinementUtilities::Cell2Ds_GeometricData::Cell2D_GeometricData::StatusTypes::Inactive)
            continue;

        const auto cell2DUnalignedPolygonType = (meshGeometricData.Cell2Ds.UnalignedVertices[cell2DToRefineIndex].cols() == 3)
                                                    ? Gedim::GeometryUtilities::PolygonTypes::Triangle
                                                    : Gedim::GeometryUtilities::PolygonTypes::Generic_Convex;

        if (mesh.Cell2DHasUpdatedCell2Ds(cell2DToRefineIndex))
        {
            Gedim::Output::Assert(!mesh.Cell2DIsActive(cell2DToRefineIndex));
            continue;
        }

        if (geometry_utilities.IsValueZero(0.5 * meshGeometricData.Cell2Ds.Area.at(cell2DToRefineIndex),
                                           geometry_utilities.Tolerance2D()))
        {
            meshGeometricData.Cell2Ds.Status[cell2DToRefineIndex] =
                Gedim::RefinementUtilities::Cell2Ds_GeometricData::Cell2D_GeometricData::StatusTypes::Inactive;
            continue;
        }

        const auto direction = ComputeTriangleMaxEdgeDirection(meshGeometricData.Cell2Ds.EdgesLength.at(cell2DToRefineIndex));

        const auto refineResult = RefineTriangleCell_ByEdge(cell2DToRefineIndex,
                                                            direction.MaxEdgeIndex,
                                                            direction.OppositeVertexIndex,
                                                            meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex),
                                                            meshGeometricData.Cell2Ds.Area.at(cell2DToRefineIndex),
                                                            Eigen::Matrix3d::Identity(),
                                                            Eigen::Vector3d::Zero(),
                                                            meshGeometricData.Cell2Ds.EdgesLength.at(cell2DToRefineIndex),
                                                            mesh);

        switch (refineResult.ResultType)
        {
        case Gedim::RefinementUtilities::RefinePolygon_Result::ResultTypes::SplitDirectionNotInsideCell2D:
            continue;
        case Gedim::RefinementUtilities::RefinePolygon_Result::ResultTypes::Cell2DAlreadySplitted: {
            Gedim::Output::Assert(!mesh.Cell2DIsActive(cell2DToRefineIndex));
        }
            continue;
        case Gedim::RefinementUtilities::RefinePolygon_Result::ResultTypes::Cell2DSplitUnderTolerance: {
            meshGeometricData.Cell2Ds.Status[cell2DToRefineIndex] =
                Gedim::RefinementUtilities::Cell2Ds_GeometricData::Cell2D_GeometricData::StatusTypes::Inactive;
        }
            continue;
        case Gedim::RefinementUtilities::RefinePolygon_Result::ResultTypes::Successfull: {
            Gedim::Output::Assert(cell2DUnalignedPolygonType == Gedim::GeometryUtilities::PolygonTypes::Triangle);

            const unsigned int cell1DIndex = mesh.Cell2DEdge(cell2DIndex, direction.MaxEdgeIndex);
            meshGeometricData.Cell1Ds.Status[cell1DIndex] =
                Gedim::RefinementUtilities::Cell2Ds_GeometricData::Cell1D_GeometricData::StatusTypes::QualityToCheck;

            refinedCells.push_back(cell2DToRefineIndex);
        }
        break;
        default:
            throw std::runtime_error("Refine result not managed");
        }

        for (unsigned int rnc = 0; rnc < refineResult.NewCell2DsIndex.size(); rnc++)
            cell2DsToUpdateGeometricData.push_back({cell2DToRefineIndex, refineResult.NewCell2DsIndex[rnc]});

        for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
        {
            if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
                continue;

            const auto newNeighboursCell2DsIndex = RefineTriangleCell_UpdateNeighbours(
                cell2DToRefineIndex,
                refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                {},
                {},
                mesh);

            for (unsigned int rnc = 0; rnc < newNeighboursCell2DsIndex.UpdatedCell2Ds.size(); rnc++)
                cell2DsToUpdateGeometricData.push_back({newNeighboursCell2DsIndex.UpdatedCell2Ds[rnc].OriginalCell2DIndex,
                                                        newNeighboursCell2DsIndex.UpdatedCell2Ds[rnc].NewCell2DIndex});
        }

        // update cell2Ds information
        std::vector<unsigned int> cell2DsToUpdate(cell2DsToUpdateGeometricData.size());
        unsigned int cPosition = 0;
        for (const Cell2DsToUpdateGeometricData &cell2DIndexToUpdate : cell2DsToUpdateGeometricData)
        {
            cell2DsToUpdate[cPosition] = cell2DIndexToUpdate.NewCell2DIndex;
            cPosition++;
        }

        RefinePolygonCell_UpdateGeometricData(mesh, cell2DsToUpdate, meshGeometricData);

        {
            // extend elements with cell1D marked
            std::unordered_set<unsigned int> cell2DsChecked;
            for (const unsigned int cell2DIndexToCheck : cell2DsToUpdate)
            {
                for (unsigned int ce = 0; ce < mesh.Cell2DNumberEdges(cell2DIndexToCheck); ce++)
                {
                    const unsigned int cell1DIndex = mesh.Cell2DEdge(cell2DIndexToCheck, ce);
                    const unsigned int aligned = meshGeometricData.Cell1Ds.Aligned.at(cell1DIndex);

                    if (meshGeometricData.Cell1Ds.Status[cell1DIndex] ==
                        Gedim::RefinementUtilities::Cell2Ds_GeometricData::Cell1D_GeometricData::StatusTypes::QualityToCheck)
                    {
                        unsigned int cell1DNumActiveNeighs = 0, cell1DQualityCheckNeighs = 0;
                        for (unsigned int c1Dn = 0; c1Dn < mesh.Cell1DNumberNeighbourCell2D(cell1DIndex); c1Dn++)
                        {
                            if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, c1Dn))
                                continue;

                            const unsigned int cell2DNeighIndex = mesh.Cell1DNeighbourCell2D(cell1DIndex, c1Dn);
                            const bool isNeighTriangle = (meshGeometricData.Cell2Ds.Vertices[cell2DNeighIndex].cols() == 3);

                            unsigned int numAligned = 0;
                            if (!isNeighTriangle)
                            {
                                for (unsigned int c1Dne = 0; c1Dne < mesh.Cell2DNumberEdges(cell2DNeighIndex); c1Dne++)
                                {
                                    const unsigned int cell1DEdgeIndex = mesh.Cell2DEdge(cell2DNeighIndex, c1Dne);

                                    if (meshGeometricData.Cell1Ds.Aligned.at(cell1DEdgeIndex) == aligned)
                                        numAligned++;
                                }

                                Gedim::Output::Assert(numAligned > 0);
                            }

                            cell1DNumActiveNeighs++;

                            if (!isNeighTriangle)
                            {
                                if (cell2DsChecked.find(cell2DNeighIndex) == cell2DsChecked.end())
                                {
                                    cell2DsChecked.insert(cell2DNeighIndex);
                                    toRefineCells.push_back(cell2DNeighIndex);
                                }
                            }
                            else
                                cell1DQualityCheckNeighs++;
                        }

                        if (cell1DNumActiveNeighs == cell1DQualityCheckNeighs)
                        {
                            meshGeometricData.Cell1Ds.Status[cell1DIndex] =
                                Gedim::RefinementUtilities::Cell2Ds_GeometricData::Cell1D_GeometricData::StatusTypes::QualityChecked;

                            unsigned int originalCell1DIndex = cell1DIndex;
                            while (mesh.Cell1DHasOriginalCell1D(originalCell1DIndex))
                            {
                                originalCell1DIndex = mesh.Cell1DOriginalCell1D(cell1DIndex);

                                if (meshGeometricData.Cell1Ds.Status[originalCell1DIndex] ==
                                    Gedim::RefinementUtilities::Cell2Ds_GeometricData::Cell1D_GeometricData::StatusTypes::QualityChecked)
                                    break;

                                std::list<unsigned int> updatedCell1Ds;
                                mesh.Cell1DUpdatedCell1Ds(originalCell1DIndex, updatedCell1Ds);

                                unsigned int childStatus = 0;
                                for (const unsigned int updateCell1D : updatedCell1Ds)
                                {
                                    if (meshGeometricData.Cell1Ds.Status[updateCell1D] ==
                                        Gedim::RefinementUtilities::Cell2Ds_GeometricData::Cell1D_GeometricData::StatusTypes::QualityToCheck)
                                        break;

                                    childStatus++;
                                }

                                if (childStatus == updatedCell1Ds.size())
                                {
                                    meshGeometricData.Cell1Ds.Status[originalCell1DIndex] =
                                        Gedim::RefinementUtilities::Cell2Ds_GeometricData::Cell1D_GeometricData::StatusTypes::QualityChecked;
                                }
                                else
                                    break;
                            }
                        }
                    }
                }
            }
        }
    }

    mesh.Compress();

    return std::vector<unsigned int>(refinedCells.begin(), refinedCells.end());
}
// ***************************************************************************
} // namespace Gedim
