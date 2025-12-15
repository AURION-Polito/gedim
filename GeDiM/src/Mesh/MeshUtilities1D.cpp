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

#include "MeshUtilities.hpp"

#include "ImportExportUtilities.hpp"
#include <fstream>
#include <set>

namespace Gedim
{
// ***************************************************************************
void MeshUtilities::FillMesh1D(const Gedim::GeometryUtilities &,
                               const Eigen::Vector3d &segmentOrigin,
                               const Eigen::Vector3d &segmentTangent,
                               const std::vector<double> &coordinates,
                               Gedim::IMeshDAO &mesh) const
{
    if (coordinates.size() == 0)
        return;

    mesh.InitializeDimension(1);

    const unsigned int &numCell0Ds = coordinates.size();
    mesh.Cell0DsInitialize(numCell0Ds);
    for (unsigned int c = 0; c < numCell0Ds; c++)
    {
        mesh.Cell0DSetState(c, true);
        mesh.Cell0DInsertCoordinates(c, segmentOrigin + coordinates[c] * segmentTangent);
    }
    mesh.Cell0DSetMarker(0, 1);
    mesh.Cell0DSetMarker(numCell0Ds - 1, 2);

    const unsigned int numCell1Ds = numCell0Ds - 1;
    mesh.Cell1DsInitialize(numCell1Ds);
    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
        mesh.Cell1DInsertExtremes(e, e, e + 1);
        mesh.Cell1DSetState(e, true);
        mesh.Cell1DSetMarker(e, 0);
    }
}
// ***************************************************************************
MeshUtilities::FilterMeshData MeshUtilities::FilterMesh1D(const std::vector<unsigned int> &cell1DsFilter,
                                                          const Gedim::IMeshDAO &mesh) const
{
    std::list<unsigned int> cell1Ds;
    std::set<unsigned int> cell0Ds;

    for (const unsigned int cell1DIndex : cell1DsFilter)
    {
        if (!mesh.Cell1DIsActive(cell1DIndex))
            continue;

        cell1Ds.push_back(cell1DIndex);

        for (unsigned int v = 0; v < 2; v++)
        {
            const unsigned int cell0DIndex = mesh.Cell1DVertex(cell1DIndex, v);

            if (cell0Ds.find(cell0DIndex) == cell0Ds.end())
                cell0Ds.insert(cell0DIndex);
        }
    }

    MeshUtilities::FilterMeshData result;

    result.Cell0Ds = std::vector<unsigned int>(cell0Ds.begin(), cell0Ds.end());
    result.Cell1Ds = std::vector<unsigned int>(cell1Ds.begin(), cell1Ds.end());

    return result;
}
// ***************************************************************************
MeshUtilities::ExtractMeshData MeshUtilities::ExtractMesh1D(const std::vector<unsigned int> &cell0DsFilter,
                                                            const std::vector<unsigned int> &cell1DsFilter,
                                                            const Gedim::IMeshDAO &originalMesh,
                                                            Gedim::IMeshDAO &mesh) const
{
    ExtractMeshData result;
    result.NewCell0DToOldCell0D.resize(cell0DsFilter.size(), std::numeric_limits<unsigned int>::max());
    result.NewCell1DToOldCell1D.resize(cell1DsFilter.size(), std::numeric_limits<unsigned int>::max());
    result.OldCell0DToNewCell0D.resize(originalMesh.Cell0DTotalNumber(), std::numeric_limits<unsigned int>::max());
    result.OldCell1DToNewCell1D.resize(originalMesh.Cell1DTotalNumber(), std::numeric_limits<unsigned int>::max());

    Eigen::MatrixXd newCell0Ds(3, cell0DsFilter.size());
    for (unsigned int v = 0; v < cell0DsFilter.size(); v++)
    {
        const unsigned int oldCell0DIndex = cell0DsFilter[v];
        result.NewCell0DToOldCell0D[v] = oldCell0DIndex;
        result.OldCell0DToNewCell0D[oldCell0DIndex] = v;

        newCell0Ds.col(v) = originalMesh.Cell0DCoordinates(oldCell0DIndex);
    }

    Eigen::MatrixXi newCell1Ds(2, cell1DsFilter.size());
    for (unsigned int e = 0; e < cell1DsFilter.size(); e++)
    {
        const unsigned int oldCell1DIndex = cell1DsFilter[e];
        result.NewCell1DToOldCell1D[e] = oldCell1DIndex;
        result.OldCell1DToNewCell1D[oldCell1DIndex] = e;

        const Eigen::VectorXi cell1DExtremes = originalMesh.Cell1DExtremes(oldCell1DIndex);

        newCell1Ds(0, e) = result.OldCell0DToNewCell0D.at(cell1DExtremes[0]);
        newCell1Ds(1, e) = result.OldCell0DToNewCell0D.at(cell1DExtremes[1]);
    }

    mesh.InitializeDimension(1);

    const unsigned int &numCell0Ds = newCell0Ds.cols();
    mesh.Cell0DsInitialize(numCell0Ds);
    mesh.Cell0DsInsertCoordinates(newCell0Ds);
    for (unsigned int c = 0; c < numCell0Ds; c++)
        mesh.Cell0DSetState(c, true);

    const unsigned int numCell1Ds = newCell1Ds.cols();
    mesh.Cell1DsInitialize(numCell1Ds);
    mesh.Cell1DsInsertExtremes(newCell1Ds);
    for (unsigned int e = 0; e < numCell1Ds; e++)
        mesh.Cell1DSetState(e, true);

    return result;
}
// ***************************************************************************
void MeshUtilities::Mesh1DFromSegment(const Gedim::GeometryUtilities &geometryUtilities,
                                      const Eigen::MatrixXd &segmentVertices,
                                      const std::vector<unsigned int> vertexMarkers,
                                      Gedim::IMeshDAO &mesh) const
{
    FillMesh1D(geometryUtilities, segmentVertices.col(0), segmentVertices.col(1), {0.0, 1.0}, mesh);

    mesh.Cell0DSetMarker(0, vertexMarkers[0]);
    mesh.Cell0DSetMarker(mesh.Cell0DTotalNumber() - 1, vertexMarkers[1]);
}
// ***************************************************************************
MeshUtilities::MeshGeometricData1D MeshUtilities::FillMesh1DGeometricData(const Gedim::GeometryUtilities &geometryUtilities,
                                                                          const Gedim::IMeshDAO &convexMesh) const
{
    MeshGeometricData1D result;

    result.Cell1DsBoundingBox.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsVertices.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsTangents.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsLengths.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsSquaredLengths.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsCentroids.resize(convexMesh.Cell1DTotalNumber());

    for (unsigned int c = 0; c < convexMesh.Cell1DTotalNumber(); c++)
    {
        if (!convexMesh.Cell1DIsActive(c))
            continue;

        // Extract original cell1D geometric information

        result.Cell1DsVertices[c] = convexMesh.Cell1DCoordinates(c);
        result.Cell1DsBoundingBox[c] = geometryUtilities.PointsBoundingBox(result.Cell1DsVertices[c]);
        result.Cell1DsTangents[c] =
            geometryUtilities.SegmentTangent(result.Cell1DsVertices[c].col(0), result.Cell1DsVertices[c].col(1));
        result.Cell1DsLengths[c] =
            geometryUtilities.SegmentLength(result.Cell1DsVertices[c].col(0), result.Cell1DsVertices[c].col(1));
        result.Cell1DsSquaredLengths[c] = result.Cell1DsLengths[c] * result.Cell1DsLengths[c];
        result.Cell1DsCentroids[c] = geometryUtilities.SimplexBarycenter(result.Cell1DsVertices[c]);
    }

    return result;
}
// ***************************************************************************
MeshUtilities::MeshGeometricData1D MeshUtilities::ImportMeshGeometricData1DFromTxt(const std::string &file_path)
{
    using Gedim_ImportExport_Utilities::operator>>;

    MeshGeometricData1D mesh_geometric_data;

    std::ifstream inFile;
    inFile.open(file_path.c_str());

    if (!inFile.is_open())
        throw std::runtime_error("Import file '" + file_path + "' not opened");

    inFile >> mesh_geometric_data.Cell1DsBoundingBox;
    inFile >> mesh_geometric_data.Cell1DsVertices;
    inFile >> mesh_geometric_data.Cell1DsTangents;
    inFile >> mesh_geometric_data.Cell1DsLengths;
    inFile >> mesh_geometric_data.Cell1DsSquaredLengths;
    inFile >> mesh_geometric_data.Cell1DsCentroids;

    return mesh_geometric_data;
}
// ***************************************************************************
void MeshUtilities::ExportMeshGeometricData1DToTxt(const MeshGeometricData1D &mesh_geometric_data, const std::string &file_path)
{
    using Gedim_ImportExport_Utilities::operator<<;

    std::ofstream file(file_path);

    if (!file.is_open())
        throw std::runtime_error("Export file '" + file_path + "' not opened");

    file.precision(16);

    file << std::scientific << mesh_geometric_data.Cell1DsBoundingBox << std::endl;
    file << std::scientific << mesh_geometric_data.Cell1DsVertices << std::endl;
    file << std::scientific << mesh_geometric_data.Cell1DsTangents << std::endl;
    file << std::scientific << mesh_geometric_data.Cell1DsLengths << std::endl;
    file << std::scientific << mesh_geometric_data.Cell1DsSquaredLengths << std::endl;
    file << std::scientific << mesh_geometric_data.Cell1DsCentroids << std::endl;

    file.close();
}
// ***************************************************************************
std::vector<unsigned int> MeshUtilities::SplitCell1D(const unsigned int &cell1DIndex,
                                                     const Eigen::MatrixXi subCell1Ds,
                                                     Gedim::IMeshDAO &mesh) const
{
    const unsigned int numSubCells = subCell1Ds.cols();
    unsigned int newCell1DsStartingIndex = mesh.Cell1DAppend(numSubCells);

    mesh.Cell1DSetState(cell1DIndex, false);

    std::vector<unsigned int> newCell1DsIndex(numSubCells);

    for (unsigned int c = 0; c < numSubCells; c++)
    {
        newCell1DsIndex[c] = newCell1DsStartingIndex + c;

        const unsigned int &newCell1DIndex = newCell1DsIndex[c];

        mesh.Cell1DInsertExtremes(newCell1DIndex, subCell1Ds(0, c), subCell1Ds(1, c));

        mesh.Cell1DSetMarker(newCell1DIndex, mesh.Cell1DMarker(cell1DIndex));
        mesh.Cell1DSetState(newCell1DIndex, true);

        mesh.Cell1DInsertUpdatedCell1D(cell1DIndex, newCell1DIndex);

        for (unsigned int v = 0; v < 2; v++)
        {
            const unsigned int cell0DIndex = mesh.Cell1DVertex(newCell1DIndex, v);

            for (unsigned int n = 0; n < mesh.Cell0DNumberNeighbourCell1D(cell0DIndex); n++)
            {
                if (!mesh.Cell0DHasNeighbourCell1D(cell0DIndex, n))
                    continue;

                if (mesh.Cell0DNeighbourCell1D(cell0DIndex, n) == cell1DIndex)
                    mesh.Cell0DInsertNeighbourCell1D(cell0DIndex, n, newCell1DIndex);
            }
        }

        const unsigned int numCell1DNumberNeighbourCell2D = mesh.Cell1DNumberNeighbourCell2D(cell1DIndex);

        if (numCell1DNumberNeighbourCell2D > 0)
        {

            mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex, numCell1DNumberNeighbourCell2D);

            for (unsigned int n = 0; n < numCell1DNumberNeighbourCell2D; n++)
            {
                if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, n))
                    continue;

                mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex, n, mesh.Cell1DNeighbourCell2D(cell1DIndex, n));
            }
        }

        const unsigned int numCell1DNumberNeighbourCell3D = mesh.Cell1DNumberNeighbourCell3D(cell1DIndex);

        if (numCell1DNumberNeighbourCell3D > 0)
        {
            mesh.Cell1DInitializeNeighbourCell3Ds(newCell1DIndex, numCell1DNumberNeighbourCell3D);

            for (unsigned int n = 0; n < numCell1DNumberNeighbourCell3D; n++)
            {
                if (!mesh.Cell1DHasNeighbourCell3D(cell1DIndex, n))
                    continue;

                mesh.Cell1DInsertNeighbourCell3D(newCell1DIndex, n, mesh.Cell1DNeighbourCell3D(cell1DIndex, n));
            }
        }
    }

    return newCell1DsIndex;
}
// ***************************************************************************
MeshUtilities::AgglomerateCell1DInformation MeshUtilities::AgglomerateCell1Ds(const std::unordered_set<unsigned int> &cell1DsIndex,
                                                                              const Gedim::IMeshDAO &mesh) const
{
    AgglomerateCell1DInformation result;

    if (cell1DsIndex.size() < 2)
        return result;

    std::unordered_set<unsigned int> agglomerated_cell0Ds;
    std::unordered_set<unsigned int> removed_cell0Ds;

    for (const unsigned int c1D_index : cell1DsIndex)
    {

        for (unsigned int v = 0; v < 2; v++)
        {
            const unsigned int c0D_index = mesh.Cell1DVertex(c1D_index, v);

            if (removed_cell0Ds.find(c0D_index) != removed_cell0Ds.end())
                continue;

            bool to_remove = true;

            unsigned int num_neigh_cell1Ds = 0;
            for (unsigned int n = 0; n < mesh.Cell0DNumberNeighbourCell1D(c0D_index); n++)
            {
                if (!mesh.Cell0DHasNeighbourCell1D(c0D_index, n))
                {
                    to_remove = false;
                    break;
                }

                const unsigned int c1D_neigh_index = mesh.Cell0DNeighbourCell1D(c0D_index, n);

                if (cell1DsIndex.find(c1D_neigh_index) == cell1DsIndex.end())
                {
                    to_remove = false;
                    break;
                }

                num_neigh_cell1Ds++;
            }

            if (num_neigh_cell1Ds < 2)
                to_remove = false;

            if (to_remove)
            {
                if (removed_cell0Ds.find(c0D_index) == removed_cell0Ds.end())
                    removed_cell0Ds.insert(c0D_index);

                continue;
            }

            if (agglomerated_cell0Ds.find(c0D_index) != agglomerated_cell0Ds.end())
                continue;

            agglomerated_cell0Ds.insert(c0D_index);
        }
    }

    result.AgglomerateCell1DVertices = std::vector<unsigned int>(agglomerated_cell0Ds.begin(), agglomerated_cell0Ds.end());
    result.SubCell1DsRemovedVertices = std::vector<unsigned int>(removed_cell0Ds.begin(), removed_cell0Ds.end());

    return result;
}
// ***************************************************************************
unsigned int MeshUtilities::AgglomerateCell1Ds(const std::unordered_set<unsigned int> &subCell1DsIndex,
                                               const std::vector<unsigned int> &agglomerateCell1DVertices,
                                               const std::vector<unsigned int> &subCell1DsRemovedCell0Ds,
                                               Gedim::IMeshDAO &mesh,
                                               std::vector<std::vector<unsigned int>> &meshCell1DsOriginalCell1Ds,
                                               const bool mantain_neigh2D_order) const

{
    if (subCell1DsIndex.size() == 0)
        return mesh.Cell1DTotalNumber();

    if (subCell1DsIndex.size() == 1)
        return *subCell1DsIndex.begin();

    const unsigned int agglomeratedCell1DIndex = mesh.Cell1DAppend(1);

    unsigned int max_marker = 0;
    std::list<unsigned int> agglomeratedCell1DConvexCells;
    std::unordered_map<unsigned int, unsigned int> neigh_cell2Ds;
    std::unordered_set<unsigned int> neigh_cell3Ds;

    for (const auto subCell1DIndex : subCell1DsIndex)
    {
        mesh.Cell1DSetState(subCell1DIndex, false);

        if (mesh.Cell1DMarker(subCell1DIndex) > max_marker)
            max_marker = mesh.Cell1DMarker(subCell1DIndex);

        const auto &convexCells = meshCell1DsOriginalCell1Ds.at(subCell1DIndex);
        for (const auto &convexCell : convexCells)
            agglomeratedCell1DConvexCells.push_back(convexCell);

        for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell2D(subCell1DIndex); n++)
        {
            if (!mesh.Cell1DHasNeighbourCell2D(subCell1DIndex, n))
                continue;

            const unsigned int neighCell2DIndex = mesh.Cell1DNeighbourCell2D(subCell1DIndex, n);

            if (neigh_cell2Ds.find(neighCell2DIndex) == neigh_cell2Ds.end())
                neigh_cell2Ds.insert(std::make_pair(neighCell2DIndex, n));
        }

        for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell3D(subCell1DIndex); n++)
        {
            if (!mesh.Cell1DHasNeighbourCell3D(subCell1DIndex, n))
                continue;

            const unsigned int neighCell3DIndex = mesh.Cell1DNeighbourCell3D(subCell1DIndex, n);

            if (neigh_cell3Ds.find(neighCell3DIndex) == neigh_cell3Ds.end())
                neigh_cell3Ds.insert(neighCell3DIndex);
        }
    }
    meshCell1DsOriginalCell1Ds.resize(meshCell1DsOriginalCell1Ds.size() + 1);
    meshCell1DsOriginalCell1Ds.at(agglomeratedCell1DIndex) =
        std::vector<unsigned int>(agglomeratedCell1DConvexCells.begin(), agglomeratedCell1DConvexCells.end());

    Gedim::Output::Assert(agglomerateCell1DVertices.size() == 2);

    mesh.Cell1DInsertExtremes(agglomeratedCell1DIndex, agglomerateCell1DVertices.at(0), agglomerateCell1DVertices.at(1));

    mesh.Cell1DSetMarker(agglomeratedCell1DIndex, max_marker);
    mesh.Cell1DSetState(agglomeratedCell1DIndex, true);

    if (!neigh_cell2Ds.empty())
    {
        if (mantain_neigh2D_order)
        {
            mesh.Cell1DInitializeNeighbourCell2Ds(agglomeratedCell1DIndex, 2);
            for (const auto &neigh : neigh_cell2Ds)
            {
                mesh.Cell1DInsertNeighbourCell2D(agglomeratedCell1DIndex, neigh.second, neigh.first);
            }
        }
        else
        {
            std::vector<unsigned int> neighs_2D;
            neighs_2D.reserve(neigh_cell2Ds.size());
            for (const auto &neigh : neigh_cell2Ds)
                neighs_2D.push_back(neigh.first);

            mesh.Cell1DInitializeNeighbourCell2Ds(agglomeratedCell1DIndex, neighs_2D);
        }
    }

    if (!neigh_cell3Ds.empty())
    {
        mesh.Cell1DInitializeNeighbourCell3Ds(agglomeratedCell1DIndex,
                                              std::vector<unsigned int>(neigh_cell3Ds.begin(), neigh_cell3Ds.end()));
    }

    for (unsigned int v = 0; v < 2; v++)
    {
        const unsigned int cell0DIndex = mesh.Cell1DVertex(agglomeratedCell1DIndex, v);

        std::list<unsigned int> neighCell1Ds;

        for (unsigned int n = 0; n < mesh.Cell0DNumberNeighbourCell1D(cell0DIndex); n++)
        {
            if (!mesh.Cell0DHasNeighbourCell1D(cell0DIndex, n))
                continue;

            const unsigned int neighCell1DIndex = mesh.Cell0DNeighbourCell1D(cell0DIndex, n);

            if (subCell1DsIndex.find(neighCell1DIndex) == subCell1DsIndex.end())
                neighCell1Ds.push_back(neighCell1DIndex);
        }

        if (neighCell1Ds.empty())
            continue;

        neighCell1Ds.push_back(agglomeratedCell1DIndex);
        mesh.Cell0DInitializeNeighbourCell1Ds(cell0DIndex, std::vector<unsigned int>(neighCell1Ds.begin(), neighCell1Ds.end()));
    }

    for (const auto cell0DIndex : subCell1DsRemovedCell0Ds)
        mesh.Cell0DSetState(cell0DIndex, false);

    return agglomeratedCell1DIndex;
}
// ***************************************************************************
bool MeshUtilities::CollapseCell1D(const unsigned int cell1D_index,
                                   Gedim::IMeshDAO &mesh) const
{
  if (cell1D_index > mesh.Cell1DTotalNumber())
    return false;
  if (!mesh.Cell1DIsActive(cell1D_index))
    return false;

  const auto cell1D_origin_index = mesh.Cell1DOrigin(cell1D_index);
  const auto cell1D_end_index = mesh.Cell1DEnd(cell1D_index);

  const auto cell1D_cell2Ds = mesh.Cell1DNeighbourCell2Ds(cell1D_index);

  if (cell1D_cell2Ds.empty())
    return false;

  const auto cell1D_cell3Ds = mesh.Cell1DNeighbourCell3Ds(cell1D_index);

  if (mesh.Dimension() == 3 &&
      cell1D_cell3Ds.empty())
    return false;

  const auto cell1D_origin_cell1Ds_neigh_vec = mesh.Cell0DNeighbourCell1Ds(cell1D_origin_index);
  std::set<unsigned int> cell1D_origin_cell1Ds_neigh(cell1D_origin_cell1Ds_neigh_vec.begin(),
                                                     cell1D_origin_cell1Ds_neigh_vec.end());
  const auto cell1D_end_cell1Ds_neigh = mesh.Cell0DNeighbourCell1Ds(cell1D_end_index);

  if (cell1D_origin_cell1Ds_neigh.empty() ||
      cell1D_end_cell1Ds_neigh.empty())
    return false;

  const auto cell1D_end_index_cell2Ds = mesh.Cell0DNeighbourCell2Ds(cell1D_end_index);

  const auto cell1D_origin_cell2Ds_neigh_vec = mesh.Cell0DNeighbourCell2Ds(cell1D_origin_index);
  std::set<unsigned int> cell1D_origin_cell2Ds_neigh(cell1D_origin_cell2Ds_neigh_vec.begin(),
                                                     cell1D_origin_cell2Ds_neigh_vec.end());

  if (cell1D_end_index_cell2Ds.empty() ||
      cell1D_origin_cell2Ds_neigh.empty())
    return false;

  const auto cell1D_end_index_cell3Ds = mesh.Cell0DNeighbourCell3Ds(cell1D_end_index);

  const auto cell1D_origin_cell3Ds_neigh_vec = mesh.Cell0DNeighbourCell3Ds(cell1D_origin_index);
  std::set<unsigned int> cell1D_origin_cell3Ds_neigh(cell1D_origin_cell3Ds_neigh_vec.begin(),
                                                     cell1D_origin_cell3Ds_neigh_vec.end());

  if (mesh.Dimension() == 3 &&
      (cell1D_end_index_cell3Ds.empty() ||
      cell1D_origin_cell3Ds_neigh.empty()))
    return false;

  for (const auto cell2D_index : cell1D_cell2Ds)
  {
    if (cell2D_index > mesh.Cell2DTotalNumber())
      continue;

    if (mesh.Cell2DNumberVertices(cell2D_index) == 3)
      return false;
  }

  for (const auto cell3D_index : cell1D_cell3Ds)
  {
    if (cell3D_index > mesh.Cell3DTotalNumber())
      continue;


    if (mesh.Cell3DNumberVertices(cell3D_index) == 4)
      return false;
  }

  // remove cell1D
  mesh.Cell1DSetState(cell1D_index, false);
  mesh.Cell0DSetState(cell1D_end_index, false);

  // remove cell1D_end_index from cell1Ds
  std::unordered_map<unsigned int, unsigned int> replace_cell1Ds_index;
  std::unordered_map<unsigned int, unsigned int> new_cell1Ds_index;
  for (const auto cell1D_neigh_index : cell1D_end_cell1Ds_neigh)
  {
    if (cell1D_neigh_index > mesh.Cell1DTotalNumber())
      continue;

    if (cell1D_neigh_index == cell1D_index)
      continue;

    if (!mesh.Cell1DIsActive(cell1D_neigh_index))
      continue;

    const auto cell1D_n_extremes = mesh.Cell1DExtremes(cell1D_neigh_index);

    Eigen::MatrixXi new_cell1D_extremes(2, 1);
    for (unsigned int e = 0; e < 2; ++e)
    {
      new_cell1D_extremes(e, 0) =
          static_cast<unsigned int>(cell1D_n_extremes[e]) == cell1D_end_index ?
            cell1D_origin_index :
            cell1D_n_extremes[e];
    }

    const auto new_cell1D_index = SplitCell1D(cell1D_neigh_index,
                                              new_cell1D_extremes,
                                              mesh);

    new_cell1Ds_index.insert(std::make_pair(cell1D_neigh_index,
                                            new_cell1D_index.at(0)));
    cell1D_origin_cell1Ds_neigh.insert(new_cell1D_index.at(0));
  }

  mesh.Cell0DInitializeNeighbourCell1Ds(cell1D_origin_index,
                                        cell1D_origin_cell1Ds_neigh.size());
  unsigned int c0D_1D_n = 0;
  for (const auto c1D_n : cell1D_origin_cell1Ds_neigh)
  {
    mesh.Cell0DInsertNeighbourCell1D(cell1D_origin_index,
                                     c0D_1D_n++,
                                     c1D_n);
  }

  // update cell0D and cell1D from cell2Ds
  std::unordered_map<unsigned int, unsigned int> new_cell2Ds_index;
  for (const auto cell2D_index : cell1D_end_index_cell2Ds)
  {
    if (cell2D_index > mesh.Cell2DTotalNumber())
      continue;

    if (!mesh.Cell2DIsActive(cell2D_index))
      continue;

    const bool has_edge_to_remove = std::find(cell1D_cell2Ds.begin(),
                                              cell1D_cell2Ds.end(),
                                              cell2D_index) != cell1D_cell2Ds.end();

    const auto cell2D_num_vertices = mesh.Cell2DNumberVertices(cell2D_index);
    const unsigned int num_new_vertices = has_edge_to_remove ?
                                            cell2D_num_vertices - 1 :
                                            cell2D_num_vertices;

    Eigen::MatrixXi new_cell2D_extremes(2,
                                        num_new_vertices);
    unsigned int n_v = 0;
    for (unsigned int v = 0; v < cell2D_num_vertices; ++v)
    {
      const auto cell2D_vertex_index = mesh.Cell2DVertex(cell2D_index,
                                                     v);
      const auto cell2D_edge_index = mesh.Cell2DEdge(cell2D_index,
                                                     v);

      if (cell2D_edge_index == cell1D_index)
        continue;

      const auto edge_found = new_cell1Ds_index.find(cell2D_edge_index);

      if (cell2D_vertex_index == cell1D_end_index)      
        new_cell2D_extremes(0, n_v) = cell1D_origin_index;
      else
        new_cell2D_extremes(0, n_v) = cell2D_vertex_index;

      if (edge_found == new_cell1Ds_index.end())
        new_cell2D_extremes(1, n_v) = cell2D_edge_index;
      else
        new_cell2D_extremes(1, n_v) = edge_found->second;

      n_v++;
    }

    const auto new_cell2D_index = SplitCell2D(cell2D_index,
                                              { new_cell2D_extremes },
                                              mesh);
    new_cell2Ds_index.insert(std::make_pair(cell2D_index,
                                            new_cell2D_index.at(0)));
    cell1D_origin_cell2Ds_neigh.insert(new_cell2D_index.at(0));
  }

  mesh.Cell0DInitializeNeighbourCell2Ds(cell1D_origin_index,
                                        cell1D_origin_cell2Ds_neigh.size());
  unsigned int c0D_2D_n = 0;
  for (const auto c2D_n : cell1D_origin_cell2Ds_neigh)
  {
    mesh.Cell0DInsertNeighbourCell2D(cell1D_origin_index,
                                     c0D_2D_n++,
                                     c2D_n);
  }

  if (mesh.Dimension() < 3)
    return true;

  // update cell0D and cell1D from cell3Ds
  for (const auto cell3D_index : cell1D_end_index_cell3Ds)
  {
    if (cell3D_index > mesh.Cell3DTotalNumber())
      continue;

    if (!mesh.Cell3DIsActive(cell3D_index))
      continue;

    const bool has_edge_to_remove = std::find(cell1D_cell3Ds.begin(),
                                              cell1D_cell3Ds.end(),
                                              cell3D_index) != cell1D_cell3Ds.end();

    std::vector<unsigned int> new_vertices;
    std::vector<unsigned int> new_edges;
    std::vector<unsigned int> new_faces;

    {
      const auto cell3D_num_vertices = mesh.Cell3DNumberVertices(cell3D_index);
      const unsigned int num_new_vertices = has_edge_to_remove ?
                                              cell3D_num_vertices - 1 :
                                              cell3D_num_vertices;
      new_vertices.resize(num_new_vertices);

      unsigned int n_v = 0;
      for (unsigned int v = 0; v < cell3D_num_vertices; ++v)
      {
        const auto cell3D_vertex_index = mesh.Cell3DVertex(cell3D_index,
                                                       v);

        if (cell3D_vertex_index == cell1D_end_index)
        {
          if (!has_edge_to_remove)
            new_vertices[n_v] = cell1D_origin_index;
          else
            continue;
        }
        else
          new_vertices[n_v] = cell3D_vertex_index;

        n_v++;
      }
    }

    {
      const auto cell3D_num_edges = mesh.Cell3DNumberEdges(cell3D_index);
      const unsigned int num_new_edges = has_edge_to_remove ?
                                              cell3D_num_edges - 1 :
                                              cell3D_num_edges;
      new_edges.resize(num_new_edges);

      unsigned int n_e = 0;
      for (unsigned int e = 0; e < cell3D_num_edges; ++e)
      {
        const auto cell3D_edge_index = mesh.Cell3DEdge(cell3D_index,
                                                       e);

        if (cell3D_edge_index == cell1D_index)
          continue;

        const auto edge_found = new_cell1Ds_index.find(cell3D_edge_index);


        if (edge_found == new_cell1Ds_index.end())
          new_edges[n_e] = cell3D_edge_index;
        else
          new_edges[n_e] = edge_found->second;

        n_e++;
      }
    }

    {
      const auto cell3D_num_faces = mesh.Cell3DNumberFaces(cell3D_index);
      new_faces.resize(cell3D_num_faces);

      for (unsigned int f = 0; f < cell3D_num_faces; ++f)
      {
        const auto cell3D_face_index = mesh.Cell3DFace(cell3D_index,
                                                       f);

        const auto face_found = new_cell2Ds_index.find(cell3D_face_index);


        if (face_found == new_cell2Ds_index.end())
          new_faces[f] = cell3D_face_index;
        else
          new_faces[f] = face_found->second;
      }
    }

    const auto new_cell3D_index = SplitCell3D(cell3D_index,
                                              { new_vertices },
                                              { new_edges },
                                              { new_faces },
                                              mesh);
    cell1D_origin_cell3Ds_neigh.insert(new_cell3D_index.at(0));
  }

  mesh.Cell0DInitializeNeighbourCell3Ds(cell1D_origin_index,
                                        cell1D_origin_cell3Ds_neigh.size());
  unsigned int c0D_3D_n = 0;
  for (const auto c3D_n : cell1D_origin_cell3Ds_neigh)
  {
    mesh.Cell0DInsertNeighbourCell3D(cell1D_origin_index,
                                     c0D_3D_n++,
                                     c3D_n);
  }

  return true;
}
// ***************************************************************************
} // namespace Gedim
