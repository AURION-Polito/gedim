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

#include "CommonUtilities.hpp"
#include "GeometryUtilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  /*GeometryUtilities::MergePolyhedronsResult GeometryUtilities::MergePolyhedrons(const std::array<GeometryUtilities::Polyhedron, 2>& polyhedrons) const
{
  MergePolyhedronsResult result;

  const unsigned int max_value = std::numeric_limits<unsigned int>::max();

  // vertices
  result.OriginalToMergedVertices[0].resize(polyhedrons[0].Vertices.cols());
  result.OriginalToMergedVertices[1].resize(polyhedrons[1].Vertices.cols());
  std::iota(std::begin(result.OriginalToMergedVertices[0]),
            std::end(result.OriginalToMergedVertices[0]),
            0);
  std::iota(std::begin(result.OriginalToMergedVertices[1]),
            std::end(result.OriginalToMergedVertices[1]),
            result.OriginalToMergedVertices[0].size());

  result.MergedToOriginalVertices.resize(polyhedrons[0].Vertices.cols() +
      polyhedrons[1].Vertices.cols());
  std::iota(std::begin(result.MergedToOriginalVertices),
            std::begin(result.MergedToOriginalVertices) + polyhedrons[0].Vertices.cols(),
            0);
  std::iota(std::begin(result.MergedToOriginalVertices) + polyhedrons[0].Vertices.cols(),
            std::end(result.MergedToOriginalVertices),
            0);

  Polyhedron& merged_polyhedron = result.MergedPolyhedron;

  merged_polyhedron.Vertices.resize(3,
                                    polyhedrons[0].Vertices.cols() +
                                    polyhedrons[1].Vertices.cols());
  merged_polyhedron.Vertices.block(0,
                                   0,
                                   3,
                                   polyhedrons[0].Vertices.cols())<< polyhedrons[0].Vertices;
  merged_polyhedron.Vertices.block(0,
                                   polyhedrons[0].Vertices.cols(),
                                   3,
                                   polyhedrons[1].Vertices.cols())<< polyhedrons[1].Vertices;

  // edges
  result.OriginalToMergedEdges[0].resize(polyhedrons[0].Edges.cols());
  result.OriginalToMergedEdges[1].resize(polyhedrons[1].Edges.cols());
  std::iota(std::begin(result.OriginalToMergedEdges[0]),
            std::end(result.OriginalToMergedEdges[0]),
            0);
  std::iota(std::begin(result.OriginalToMergedEdges[1]),
            std::end(result.OriginalToMergedEdges[1]),
            result.OriginalToMergedEdges[0].size());
  result.MergedToOriginalEdges.resize(polyhedrons[0].Edges.cols() +
      polyhedrons[1].Edges.cols());
  std::iota(std::begin(result.MergedToOriginalEdges),
            std::begin(result.MergedToOriginalEdges) + polyhedrons[0].Edges.cols(),
            0);
  std::iota(std::begin(result.MergedToOriginalEdges) + polyhedrons[0].Edges.cols(),
            std::end(result.MergedToOriginalEdges),
            0);

  merged_polyhedron.Edges.resize(2,
                                 polyhedrons[0].Edges.cols() +
                                 polyhedrons[1].Edges.cols());
  for (unsigned int e = 0; e < polyhedrons[0].Edges.cols(); ++e)
  {
    merged_polyhedron.Edges(0, e) = result.OriginalToMergedVertices[0][polyhedrons[0].Edges(0, e)];
    merged_polyhedron.Edges(1, e) = result.OriginalToMergedVertices[0][polyhedrons[0].Edges(1, e)];
  }
  for (unsigned int e = 0; e < polyhedrons[1].Edges.cols(); ++e)
  {
    merged_polyhedron.Edges(0, polyhedrons[0].Edges.cols() + e) = result.OriginalToMergedVertices[1][polyhedrons[1].Edges(0, e)];
    merged_polyhedron.Edges(1, polyhedrons[0].Edges.cols() + e) = result.OriginalToMergedVertices[1][polyhedrons[1].Edges(1, e)];
  }

  // faces
  result.OriginalToMergedFaces[0].resize(polyhedrons[0].Faces.size());
  result.OriginalToMergedFaces[1].resize(polyhedrons[1].Faces.size());
  std::iota(std::begin(result.OriginalToMergedFaces[0]),
            std::end(result.OriginalToMergedFaces[0]),
            0);
  std::iota(std::begin(result.OriginalToMergedFaces[1]),
            std::end(result.OriginalToMergedFaces[1]),
            result.OriginalToMergedFaces[0].size());
  result.MergedToOriginalFaces.resize(polyhedrons[0].Faces.size() +
      polyhedrons[1].Faces.size());
  std::iota(std::begin(result.MergedToOriginalFaces),
            std::begin(result.MergedToOriginalFaces) + polyhedrons[0].Faces.size(),
            0);
  std::iota(std::begin(result.MergedToOriginalFaces) + polyhedrons[0].Faces.size(),
            std::end(result.MergedToOriginalFaces),
            0);

  merged_polyhedron.Faces.reserve(polyhedrons[0].Faces.size() +
                                  polyhedrons[1].Faces.size());
  for (unsigned int f = 0; f < polyhedrons[0].Faces.size(); ++f)
  {
    const auto& original_face = polyhedrons[0].Faces[f];
    merged_polyhedron.Faces.push_back(original_face);
    for (unsigned int f_v = 0; f_v < original_face.cols(); ++f_v)
    {
      merged_polyhedron.Faces[f](0, f_v) = result.OriginalToMergedVertices[0][original_face(0, f_v)];
      merged_polyhedron.Faces[f](1, f_v) = result.OriginalToMergedEdges[0][original_face(1, f_v)];
    }
  }
  for (unsigned int f = 0; f < polyhedrons[1].Faces.size(); ++f)
  {
    const auto& original_face = polyhedrons[1].Faces[f];
    merged_polyhedron.Faces.push_back(original_face);
    for (unsigned int f_v = 0; f_v < original_face.cols(); ++f_v)
    {
      merged_polyhedron.Faces[polyhedrons[0].Faces.size() + f](0, f_v) = result.OriginalToMergedVertices[1][original_face(0, f_v)];
      merged_polyhedron.Faces[polyhedrons[0].Faces.size() + f](1, f_v) = result.OriginalToMergedEdges[1][original_face(1, f_v)];
    }
  }

  return result;
}*/
  // ***************************************************************************
  GeometryUtilities::MergePolyhedronsResult GeometryUtilities::MergePolyhedrons(const std::array<Polyhedron, 2>& polyhedrons,
                                                                                const MergePolyhedronsInput& merge_information) const
  {
    MergePolyhedronsResult result;

    const auto& poly_one = polyhedrons.at(0);
    const auto& poly_two = polyhedrons.at(1);

    // vertices
    auto& original_to_merged_vertices_one = result.OriginalToMergedVertices.at(0);
    auto& original_to_merged_vertices_two = result.OriginalToMergedVertices.at(1);
    original_to_merged_vertices_one.resize(poly_one.Vertices.cols(),
                                           MergePolyhedronsResult::none);
    original_to_merged_vertices_two.resize(poly_two.Vertices.cols(),
                                           MergePolyhedronsResult::none);

    std::list<std::array<unsigned int, 2>> merged_to_original_vertices;

    // vertices
    const bool has_vertices_type = (merge_information.Vertices_Type.at(0).size() ==  poly_one.Vertices.cols());

    unsigned int common_vertices = 0;
    for (unsigned int v = 0; v < poly_one.Vertices.cols(); ++v)
    {
      const MergePolyhedronsInput::MergeTypes vertex_type = !has_vertices_type ?
                                                              MergePolyhedronsInput::MergeTypes::None :
                                                              merge_information.Vertices_Type.at(0).at(v).first;

      switch (vertex_type)
      {
        case MergePolyhedronsInput::MergeTypes::None:
        {
          merged_to_original_vertices.push_back({ v, MergePolyhedronsResult::none });
          original_to_merged_vertices_one[v] = merged_to_original_vertices.size();
        }
          break;
        case MergePolyhedronsInput::MergeTypes::Common:
        {
          const unsigned int common_vertex_index = merge_information.Vertices_Type.at(0).at(v).second;
          const auto& common_vertices = merge_information.Common_vertices.at(common_vertex_index);
          original_to_merged_vertices_one[common_vertices[0]] = merged_to_original_vertices.size();
          original_to_merged_vertices_two[common_vertices[1]] = merged_to_original_vertices.size();
          merged_to_original_vertices.push_back(common_vertices);
        }
          break;
        case MergePolyhedronsInput::MergeTypes::Remove:
        {
        }
          break;
        default:
          throw std::runtime_error("unknown vertex type");
      }
    }


    return result;
  }
  // ***************************************************************************
} // namespace Gedim
