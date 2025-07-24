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
  GeometryUtilities::MergePolyhedronsInput GeometryUtilities::MergePolyhedronByFace(const std::array<Polyhedron, 2>& polyhedrons,
                                                                                    const std::array<unsigned int, 2> polyhedrons_common_face_index,
                                                                                    const bool remove_common_face) const
  {
    MergePolyhedronsInput result;

    const auto& poly_one = polyhedrons.at(0);
    const auto& poly_two = polyhedrons.at(1);

    // vertices
    result.Vertices_Type[0].resize(poly_one.Vertices.cols(),
                                   {
                                     MergePolyhedronsInput::MergeTypes::None,
                                     MergePolyhedronsInput::none
                                   });
    result.Vertices_Type[1].resize(poly_two.Vertices.cols(),
                                   {
                                     MergePolyhedronsInput::MergeTypes::None,
                                     MergePolyhedronsInput::none
                                   });

    // edges
    result.Edges_Type[0].resize(poly_one.Edges.cols(),
                                {
                                  MergePolyhedronsInput::MergeTypes::None,
                                  MergePolyhedronsInput::none
                                });
    result.Edges_Type[1].resize(poly_two.Edges.cols(),
                                {
                                  MergePolyhedronsInput::MergeTypes::None,
                                  MergePolyhedronsInput::none
                                });

    // faces
    result.Faces_Type[0].resize(poly_one.Faces.size(),
                                {
                                  MergePolyhedronsInput::MergeTypes::None,
                                  MergePolyhedronsInput::none
                                });
    result.Faces_Type[1].resize(poly_two.Faces.size(),
                                {
                                  MergePolyhedronsInput::MergeTypes::None,
                                  MergePolyhedronsInput::none
                                });

    // common face
    const auto& common_face_one = poly_one.Faces.at(polyhedrons_common_face_index.at(0));
    const auto& common_face_two = poly_one.Faces.at(polyhedrons_common_face_index.at(1));

    Gedim::Output::Assert(common_face_one.cols() == common_face_two.cols());
    const unsigned int size_common_face = common_face_one.cols();

    result.Common_vertices.resize(size_common_face);
    result.Common_edges.resize(size_common_face);

    for (unsigned int f_v = 0; f_v < size_common_face; ++f_v)
    {
      const unsigned int vertex_index_one = common_face_one(0, f_v);
      const unsigned int vertex_index_two = common_face_two(0, f_v);

      result.Vertices_Type[0][vertex_index_one] = { MergePolyhedronsInput::MergeTypes::Common, f_v };
      result.Vertices_Type[1][vertex_index_two] = { MergePolyhedronsInput::MergeTypes::Common, f_v };

      result.Common_vertices[f_v] = { vertex_index_one, vertex_index_two };

      const unsigned int edge_index_one = common_face_one(1, f_v);
      const unsigned int edge_index_two = common_face_two(1, f_v);

      result.Edges_Type[0][edge_index_one] = { MergePolyhedronsInput::MergeTypes::Common, f_v };
      result.Edges_Type[1][edge_index_two] = { MergePolyhedronsInput::MergeTypes::Common, f_v };

      result.Common_edges[f_v] = { edge_index_one, edge_index_two };
    }

    if (remove_common_face)
    {
      result.Faces_Type[0][polyhedrons_common_face_index.at(0)] = { MergePolyhedronsInput::MergeTypes::Remove, MergePolyhedronsInput::none };
      result.Faces_Type[1][polyhedrons_common_face_index.at(1)] = { MergePolyhedronsInput::MergeTypes::Remove, MergePolyhedronsInput::none };
    }
    else
    {
      result.Faces_Type[0][polyhedrons_common_face_index.at(0)] = { MergePolyhedronsInput::MergeTypes::Common, 0 };
      result.Faces_Type[1][polyhedrons_common_face_index.at(1)] = { MergePolyhedronsInput::MergeTypes::Common, 0 };
      result.Common_faces.resize(1);
      result.Common_faces[0] = polyhedrons_common_face_index;
    }

    return result;
  }
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

    const bool poly_one_has_vertices_type = (static_cast<unsigned int>(merge_information.Vertices_Type.at(0).size()) ==
                                             static_cast<unsigned int>(poly_one.Vertices.cols()));

    for (unsigned int v = 0; v < poly_one.Vertices.cols(); ++v)
    {
      const MergePolyhedronsInput::MergeTypes vertex_type = !poly_one_has_vertices_type ?
                                                              MergePolyhedronsInput::MergeTypes::None :
                                                              merge_information.Vertices_Type.at(0).at(v).first;

      switch (vertex_type)
      {
        case MergePolyhedronsInput::MergeTypes::None:
        {
          original_to_merged_vertices_one[v] = merged_to_original_vertices.size();
          merged_to_original_vertices.push_back({ v, MergePolyhedronsResult::none });
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
          break;
        default:
          throw std::runtime_error("unknown vertex type");
      }
    }

    const bool poly_two_has_vertices_type = (static_cast<unsigned int>(merge_information.Vertices_Type.at(1).size()) ==
                                             static_cast<unsigned int>(poly_two.Vertices.cols()));

    for (unsigned int v = 0; v < poly_two.Vertices.cols(); ++v)
    {
      const MergePolyhedronsInput::MergeTypes vertex_type = !poly_two_has_vertices_type ?
                                                              MergePolyhedronsInput::MergeTypes::None :
                                                              merge_information.Vertices_Type.at(1).at(v).first;

      switch (vertex_type)
      {
        case MergePolyhedronsInput::MergeTypes::None:
        {
          original_to_merged_vertices_two[v] = merged_to_original_vertices.size();
          merged_to_original_vertices.push_back({ MergePolyhedronsResult::none, v });
        }
          break;
        case MergePolyhedronsInput::MergeTypes::Common:
          break;
        case MergePolyhedronsInput::MergeTypes::Remove:
          break;
        default:
          throw std::runtime_error("unknown vertex type");
      }
    }

    result.MergedToOriginalVertices = std::vector<std::array<unsigned int, 2>>(merged_to_original_vertices.begin(),
                                                                               merged_to_original_vertices.end());
    result.MergedPolyhedron.Vertices.resize(3, merged_to_original_vertices.size());
    for (unsigned int m_v = 0; m_v < result.MergedToOriginalVertices.size(); ++m_v)
    {
      const auto& merged_vertex = result.MergedToOriginalVertices.at(m_v);

      if (merged_vertex.at(0) != MergePolyhedronsResult::none)
        result.MergedPolyhedron.Vertices.col(m_v)<< poly_one.Vertices.col(merged_vertex.at(0));
      else if (merged_vertex.at(1) != MergePolyhedronsResult::none)
        result.MergedPolyhedron.Vertices.col(m_v)<< poly_two.Vertices.col(merged_vertex.at(1));
      else
        throw std::runtime_error("Error on merged vertex creation");
    }

    // edges
    auto& original_to_merged_edges_one = result.OriginalToMergedEdges.at(0);
    auto& original_to_merged_edges_two = result.OriginalToMergedEdges.at(1);
    original_to_merged_edges_one.resize(poly_one.Edges.cols(),
                                        MergePolyhedronsResult::none);
    original_to_merged_edges_two.resize(poly_two.Edges.cols(),
                                        MergePolyhedronsResult::none);

    std::list<std::array<unsigned int, 2>> merged_to_original_edges;

    const bool poly_one_has_edges_type = (static_cast<unsigned int>(merge_information.Edges_Type.at(0).size()) ==
                                          static_cast<unsigned int>(poly_one.Edges.cols()));

    for (unsigned int e = 0; e < poly_one.Edges.cols(); ++e)
    {
      const MergePolyhedronsInput::MergeTypes edge_type = !poly_one_has_edges_type ?
                                                            MergePolyhedronsInput::MergeTypes::None :
                                                            merge_information.Edges_Type.at(0).at(e).first;

      switch (edge_type)
      {
        case MergePolyhedronsInput::MergeTypes::None:
        {
          original_to_merged_edges_one[e] = merged_to_original_edges.size();
          merged_to_original_edges.push_back({ e, MergePolyhedronsResult::none });
        }
          break;
        case MergePolyhedronsInput::MergeTypes::Common:
        {
          const unsigned int common_edge_index = merge_information.Edges_Type.at(0).at(e).second;
          const auto& common_edges = merge_information.Common_edges.at(common_edge_index);

          original_to_merged_edges_one[common_edges[0]] = merged_to_original_edges.size();
          original_to_merged_edges_two[common_edges[1]] = merged_to_original_edges.size();
          merged_to_original_edges.push_back(common_edges);
        }
          break;
        case MergePolyhedronsInput::MergeTypes::Remove:
          break;
        default:
          throw std::runtime_error("unknown edge type");
      }
    }

    const bool poly_two_has_edges_type = (static_cast<unsigned int>(merge_information.Edges_Type.at(1).size()) ==
                                          static_cast<unsigned int>(poly_two.Edges.cols()));

    for (unsigned int e = 0; e < poly_two.Edges.cols(); ++e)
    {
      const MergePolyhedronsInput::MergeTypes edge_type = !poly_two_has_edges_type ?
                                                            MergePolyhedronsInput::MergeTypes::None :
                                                            merge_information.Edges_Type.at(1).at(e).first;

      switch (edge_type)
      {
        case MergePolyhedronsInput::MergeTypes::None:
        {
          original_to_merged_edges_two[e] = merged_to_original_edges.size();
          merged_to_original_edges.push_back({ MergePolyhedronsResult::none, e });
        }
          break;
        case MergePolyhedronsInput::MergeTypes::Common:
          break;
        case MergePolyhedronsInput::MergeTypes::Remove:
          break;
        default:
          throw std::runtime_error("unknown edge type");
      }
    }

    result.MergedToOriginalEdges = std::vector<std::array<unsigned int, 2>>(merged_to_original_edges.begin(),
                                                                            merged_to_original_edges.end());
    result.MergedPolyhedron.Edges.resize(2, merged_to_original_edges.size());
    for (unsigned int m_e = 0; m_e < result.MergedToOriginalEdges.size(); ++m_e)
    {
      const auto& merged_edge = result.MergedToOriginalEdges.at(m_e);

      if (merged_edge.at(0) != MergePolyhedronsResult::none)
      {
        const Eigen::VectorXi original_edge = poly_one.Edges.col(merged_edge.at(0));
        result.MergedPolyhedron.Edges(0, m_e) = original_to_merged_vertices_one.at(original_edge[0]);
        result.MergedPolyhedron.Edges(1, m_e) = original_to_merged_vertices_one.at(original_edge[1]);
      }
      else if (merged_edge.at(1) != MergePolyhedronsResult::none)
      {
        const Eigen::VectorXi original_edge = poly_two.Edges.col(merged_edge.at(1));
        result.MergedPolyhedron.Edges(0, m_e) = original_to_merged_vertices_two.at(original_edge[0]);
        result.MergedPolyhedron.Edges(1, m_e) = original_to_merged_vertices_two.at(original_edge[1]);
      }
      else
        throw std::runtime_error("Error on merged edge creation");
    }

    // faces
    auto& original_to_merged_faces_one = result.OriginalToMergedFaces.at(0);
    auto& original_to_merged_faces_two = result.OriginalToMergedFaces.at(1);
    original_to_merged_faces_one.resize(poly_one.Faces.size(),
                                        MergePolyhedronsResult::none);
    original_to_merged_faces_two.resize(poly_two.Faces.size(),
                                        MergePolyhedronsResult::none);

    std::list<std::array<unsigned int, 2>> merged_to_original_faces;

    const bool poly_one_has_faces_type = (static_cast<unsigned int>(merge_information.Faces_Type.at(0).size()) ==
                                          static_cast<unsigned int>(poly_one.Faces.size()));

    for (unsigned int f = 0; f < poly_one.Faces.size(); ++f)
    {
      const MergePolyhedronsInput::MergeTypes face_type = !poly_one_has_faces_type ?
                                                            MergePolyhedronsInput::MergeTypes::None :
                                                            merge_information.Faces_Type.at(0).at(f).first;

      switch (face_type)
      {
        case MergePolyhedronsInput::MergeTypes::None:
        {
          original_to_merged_faces_one[f] = merged_to_original_faces.size();
          merged_to_original_faces.push_back({ f, MergePolyhedronsResult::none });
        }
          break;
        case MergePolyhedronsInput::MergeTypes::Common:
        {
          const unsigned int common_face_index = merge_information.Faces_Type.at(0).at(f).second;
          const auto& common_faces = merge_information.Common_faces.at(common_face_index);

          original_to_merged_faces_one[common_faces[0]] = merged_to_original_faces.size();
          original_to_merged_faces_two[common_faces[1]] = merged_to_original_faces.size();
          merged_to_original_faces.push_back(common_faces);
        }
        case MergePolyhedronsInput::MergeTypes::Remove:
          break;
        default:
          throw std::runtime_error("unknown face type");
      }
    }

    const bool poly_two_has_faces_type = (static_cast<unsigned int>(merge_information.Faces_Type.at(1).size()) ==
                                          static_cast<unsigned int>(poly_two.Faces.size()));

    for (unsigned int f = 0; f < poly_two.Faces.size(); ++f)
    {
      const MergePolyhedronsInput::MergeTypes face_type = !poly_two_has_faces_type ?
                                                            MergePolyhedronsInput::MergeTypes::None :
                                                            merge_information.Faces_Type.at(1).at(f).first;

      switch (face_type)
      {
        case MergePolyhedronsInput::MergeTypes::None:
        {
          original_to_merged_faces_two[f] = merged_to_original_faces.size();
          merged_to_original_faces.push_back({ MergePolyhedronsResult::none, f });
        }
          break;
        case MergePolyhedronsInput::MergeTypes::Common:
          break;
        case MergePolyhedronsInput::MergeTypes::Remove:
          break;
        default:
          throw std::runtime_error("unknown face type");
      }
    }

    result.MergedToOriginalFaces = std::vector<std::array<unsigned int, 2>>(merged_to_original_faces.begin(),
                                                                            merged_to_original_faces.end());
    result.MergedPolyhedron.Faces.resize(merged_to_original_faces.size());
    for (unsigned int m_f = 0; m_f < result.MergedToOriginalFaces.size(); ++m_f)
    {
      const auto& merged_face = result.MergedToOriginalFaces.at(m_f);
      auto& face = result.MergedPolyhedron.Faces.at(m_f);

      if (merged_face.at(0) != MergePolyhedronsResult::none)
      {
        const auto original_face = poly_one.Faces.at(merged_face.at(0));
        face.resize(2, original_face.cols());
        for (unsigned int f_v = 0; f_v < original_face.cols(); ++f_v)
        {
          face(0, f_v) = original_to_merged_vertices_one.at(original_face(0, f_v));
          face(1, f_v) = original_to_merged_edges_one.at(original_face(1, f_v));
        }
      }
      else if (merged_face.at(1) != MergePolyhedronsResult::none)
      {
        const auto original_face = poly_two.Faces.at(merged_face.at(1));
        face.resize(2, original_face.cols());
        for (unsigned int f_v = 0; f_v < original_face.cols(); ++f_v)
        {
          face(0, f_v) = original_to_merged_vertices_two.at(original_face(0, f_v));
          face(1, f_v) = original_to_merged_edges_two.at(original_face(1, f_v));
        }
      }
      else
        throw std::runtime_error("Error on merged face creation");
    }

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
