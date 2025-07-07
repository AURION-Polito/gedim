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
GeometryUtilities::MergePolyhedronsResult GeometryUtilities::MergePolyhedrons(const std::array<GeometryUtilities::Polyhedron, 2>& polyhedrons) const
{
  MergePolyhedronsResult result;

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
}
// ***************************************************************************
GeometryUtilities::MergePolyhedronsResult GeometryUtilities::MergePolyhedronsByFace(const std::array<GeometryUtilities::Polyhedron, 2>& polyhedrons,
                                                                        const std::array<unsigned int, 2>& polyhedrons_face_index) const
{
  MergePolyhedronsResult result;

  return result;
}
// ***************************************************************************
} // namespace Gedim
