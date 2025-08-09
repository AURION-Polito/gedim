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

#ifndef __TEST_GEOMETRY_MERGED_POLYHEDRON_H
#define __TEST_GEOMETRY_MERGED_POLYHEDRON_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"
#include "PlatonicSolid.hpp"
#include "VTKUtilities.hpp"

namespace GedimUnitTesting
{

  TEST(TestGeometryUtilities, Test_MergePolyhedrons_no_intersections)
  {
    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities meshUtilities;

    const Gedim::PlatonicSolid platonicSolid = Gedim::PlatonicSolid(geometryUtilities,
                                                                    meshUtilities);

    const auto polyhedron_one = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(-2.0, -2.0, -2.0),
                                                                       4.0);
    const auto polyhedron_two = platonicSolid.dodecahedron();
    const std::array<Gedim::GeometryUtilities::Polyhedron, 2> polyhedrons = {
      polyhedron_one,
      polyhedron_two
    };

    const auto merged_polyhedron = geometryUtilities.MergePolyhedrons(polyhedrons);
    std::string exportFolder = "./Export/TestGeometryUtilities/Test_MergePolyhedrons_no_intersections";
    Gedim::Output::CreateFolder(exportFolder);
    Gedim::Output::CreateFolder(exportFolder + "/polyhedron_one");
    geometryUtilities.ExportPolyhedronToVTU(polyhedron_one,
                                            exportFolder + "/polyhedron_one");
    Gedim::Output::CreateFolder(exportFolder + "/polyhedron_two");
    geometryUtilities.ExportPolyhedronToVTU(polyhedron_two,
                                            exportFolder + "/polyhedron_two");
    Gedim::Output::CreateFolder(exportFolder + "/merged_polyhedron");
    geometryUtilities.ExportPolyhedronToVTU(merged_polyhedron.MergedPolyhedron,
                                            exportFolder + "/merged_polyhedron");


    ASSERT_EQ(merged_polyhedron.MergedToOriginalVertices.size(),
              polyhedrons[0].Vertices.cols() +
        polyhedrons[1].Vertices.cols());
    ASSERT_EQ(merged_polyhedron.MergedToOriginalEdges.size(),
              polyhedrons[0].Edges.cols() +
        polyhedrons[1].Edges.cols());
    ASSERT_EQ(merged_polyhedron.MergedToOriginalFaces.size(),
              polyhedrons[0].Faces.size() +
        polyhedrons[1].Faces.size());
    ASSERT_EQ(merged_polyhedron.MergedPolyhedron.Vertices.cols(),
              polyhedrons[0].Vertices.cols() +
        polyhedrons[1].Vertices.cols());
    ASSERT_EQ(merged_polyhedron.MergedPolyhedron.Edges.cols(),
              polyhedrons[0].Edges.cols() +
        polyhedrons[1].Edges.cols());
    ASSERT_EQ(merged_polyhedron.MergedPolyhedron.Faces.size(),
              polyhedrons[0].Faces.size() +
        polyhedrons[1].Faces.size());
    for (unsigned int p = 0; p < polyhedrons.size(); ++p)
    {
      const auto& polyhedron = polyhedrons[p];
      ASSERT_EQ(merged_polyhedron.OriginalToMergedVertices[p].size(),
                polyhedron.Vertices.cols());
      ASSERT_EQ(merged_polyhedron.OriginalToMergedEdges[p].size(),
                polyhedron.Edges.cols());
      ASSERT_EQ(merged_polyhedron.OriginalToMergedFaces[p].size(),
                polyhedron.Faces.size());

      const unsigned int shift_vertices = (p == 0) ? 0 :
                                                     polyhedrons[p - 1].Vertices.cols();
      const unsigned int shift_edges = (p == 0) ? 0 :
                                                  polyhedrons[p - 1].Edges.cols();
      const unsigned int shift_faces = (p == 0) ? 0 :
                                                  polyhedrons[p - 1].Faces.size();

      for (unsigned int v = 0; v < polyhedron.Vertices.cols(); ++v)
      {
        ASSERT_EQ(shift_vertices + v,
                  merged_polyhedron.OriginalToMergedVertices[p][v]);
        ASSERT_EQ(v,
                  merged_polyhedron.MergedToOriginalVertices[shift_vertices + v][p]);
      }
      for (unsigned int e = 0; e < polyhedron.Edges.cols(); ++e)
      {
        ASSERT_EQ(shift_edges + e,
                  merged_polyhedron.OriginalToMergedEdges[p][e]);
        ASSERT_EQ(e,
                  merged_polyhedron.MergedToOriginalEdges[shift_edges + e][p]);
      }
      for (unsigned int f = 0; f < polyhedron.Faces.size(); ++f)
      {
        ASSERT_EQ(shift_faces + f,
                  merged_polyhedron.OriginalToMergedFaces[p][f]);
        ASSERT_EQ(f,
                  merged_polyhedron.MergedToOriginalFaces[shift_faces + f][p]);
      }
    }
  }

  TEST(TestGeometryUtilities, Test_MergePolyhedrons_face_intersections)
  {
    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const auto polyhedron_one = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                       1.0);
    const auto polyhedron_two = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                       1.0);
    const std::array<Gedim::GeometryUtilities::Polyhedron, 2> polyhedrons = {
      polyhedron_one,
      polyhedron_two
    };

    const auto merged_polyhedron_input = geometryUtilities.MergePolyhedronByFace(polyhedrons,
                                                                                 { 1, 0 },
                                                                                 true);

    const std::vector<std::array<unsigned int, 2>> expected_common_vertices = { { 4u, 0u }, { 5u, 1u }, { 6u, 2u }, { 7u, 3u } };
    const std::vector<std::array<unsigned int, 2>> expected_common_edges = { { 4u, 0u }, { 5u, 1u }, { 6u, 2u }, { 7u, 3u } };

    for (unsigned int p = 0; p < polyhedrons.size(); ++p)
    {
      const auto& polyhedron = polyhedrons[p];

      ASSERT_EQ(merged_polyhedron_input.Vertices_Type[p].size(),
                polyhedron.Vertices.cols());
      ASSERT_EQ(merged_polyhedron_input.Edges_Type[p].size(),
                polyhedron.Edges.cols());
      ASSERT_EQ(merged_polyhedron_input.Faces_Type[p].size(),
                polyhedron.Faces.size());

      for (unsigned int c_v = 0; c_v < expected_common_vertices.size(); ++c_v)
      {
        const auto& expected_common_vertex = expected_common_vertices[c_v];
        ASSERT_EQ(merged_polyhedron_input.Vertices_Type[p][expected_common_vertex.at(p)],
            std::make_pair(Gedim::GeometryUtilities::MergePolyhedronsInput::MergeTypes::Common, c_v));
      }

      for (unsigned int c_e = 0; c_e < expected_common_edges.size(); ++c_e)
      {
        const auto& expected_common_edge = expected_common_edges[c_e];
        ASSERT_EQ(merged_polyhedron_input.Edges_Type[p][expected_common_edge.at(p)],
            std::make_pair(Gedim::GeometryUtilities::MergePolyhedronsInput::MergeTypes::Common, c_e));
      }
    }

    ASSERT_EQ(merged_polyhedron_input.Common_vertices,
              expected_common_vertices);
    ASSERT_EQ(merged_polyhedron_input.Common_edges,
              expected_common_edges);
    ASSERT_EQ(merged_polyhedron_input.Faces_Type[0][1],
        std::make_pair(Gedim::GeometryUtilities::MergePolyhedronsInput::MergeTypes::Remove, Gedim::GeometryUtilities::MergePolyhedronsInput::none));
    ASSERT_EQ(merged_polyhedron_input.Faces_Type[1][0],
        std::make_pair(Gedim::GeometryUtilities::MergePolyhedronsInput::MergeTypes::Remove, Gedim::GeometryUtilities::MergePolyhedronsInput::none));


    const auto merged_polyhedron = geometryUtilities.MergePolyhedrons(polyhedrons,
                                                                      merged_polyhedron_input);
    std::string exportFolder = "./Export/TestGeometryUtilities/Test_MergePolyhedrons_face_intersections";
    Gedim::Output::CreateFolder(exportFolder);
    Gedim::Output::CreateFolder(exportFolder + "/polyhedron_one");
    geometryUtilities.ExportPolyhedronToVTU(polyhedron_one,
                                            exportFolder + "/polyhedron_one");
    Gedim::Output::CreateFolder(exportFolder + "/polyhedron_two");
    geometryUtilities.ExportPolyhedronToVTU(polyhedron_two,
                                            exportFolder + "/polyhedron_two");
    Gedim::Output::CreateFolder(exportFolder + "/merged_polyhedron");
    geometryUtilities.ExportPolyhedronToVTU(merged_polyhedron.MergedPolyhedron,
                                            exportFolder + "/merged_polyhedron");


    ASSERT_EQ(merged_polyhedron.MergedPolyhedron.Vertices.cols(),
              polyhedrons[0].Vertices.cols() +
        polyhedrons[1].Vertices.cols() -
        merged_polyhedron_input.Common_vertices.size());
    ASSERT_EQ(merged_polyhedron.MergedPolyhedron.Edges.cols(),
              polyhedrons[0].Edges.cols() +
        polyhedrons[1].Edges.cols() -
        merged_polyhedron_input.Common_edges.size());
    ASSERT_EQ(merged_polyhedron.MergedPolyhedron.Faces.size(),
              polyhedrons[0].Faces.size() +
        polyhedrons[1].Faces.size() - 2);

    std::array<std::vector<unsigned int>, 2> expected_original_to_merged_vertices;
    expected_original_to_merged_vertices[0] = { 0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u };
    expected_original_to_merged_vertices[1] = { 4u, 5u, 6u, 7u, 8u, 9u, 10u, 11u };

    std::array<std::vector<unsigned int>, 2> expected_original_to_merged_edges;
    expected_original_to_merged_edges[0] = { 0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u, 11u };
    expected_original_to_merged_edges[1] = { 4u, 5u, 6u, 7u, 12u, 13u, 14u, 15u, 16u, 17u, 18u, 19u };

    std::array<std::vector<unsigned int>, 2> expected_original_to_merged_faces;
    expected_original_to_merged_faces[0] = { 0u, Gedim::GeometryUtilities::MergePolyhedronsResult::none, 1u, 2u, 3u, 4u };
    expected_original_to_merged_faces[1] = { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 5u, 6u, 7u, 8u, 9u };

    ASSERT_EQ(merged_polyhedron.OriginalToMergedVertices,
              expected_original_to_merged_vertices);
    ASSERT_EQ(merged_polyhedron.OriginalToMergedEdges,
              expected_original_to_merged_edges);
    ASSERT_EQ(merged_polyhedron.OriginalToMergedFaces,
              expected_original_to_merged_faces);

    const std::vector<std::array<unsigned int, 2>> expected_merged_to_original_vertices =
    {
      { 0u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 1u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 2u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 3u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 4u, 0u },
      { 5u, 1u },
      { 6u, 2u },
      { 7u, 3u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 4u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 5u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 6u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 7u }
    };
    const std::vector<std::array<unsigned int, 2>> expected_merged_to_original_edges =
    {
      { 0u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 1u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 2u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 3u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 4u, 0u },
      { 5u, 1u },
      { 6u, 2u },
      { 7u, 3u },
      { 8u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 9u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 10u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 11u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 4u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 5u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 6u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 7u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 8u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 9u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 10u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 11u }
    };
    const std::vector<std::array<unsigned int, 2>> expected_merged_to_original_faces =
    {
      { 0u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 2u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 3u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 4u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { 5u, Gedim::GeometryUtilities::MergePolyhedronsResult::none },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 1u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 2u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 3u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 4u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 5u }
    };

    ASSERT_EQ(merged_polyhedron.MergedToOriginalVertices,
              expected_merged_to_original_vertices);
    ASSERT_EQ(merged_polyhedron.MergedToOriginalEdges,
              expected_merged_to_original_edges);
    ASSERT_EQ(merged_polyhedron.MergedToOriginalFaces,
              expected_merged_to_original_faces);

  }

  TEST(TestGeometryUtilities, Test_MergePolyhedrons_edge_removed)
  {
    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const auto polyhedron_one = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                       1.0);
    const auto polyhedron_two = geometryUtilities.CreateTetrahedronWithVertices(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                Eigen::Vector3d(0.0, 1.0, 0.0),
                                                                                Eigen::Vector3d(0.0, 0.0, 1.0));
    const std::array<Gedim::GeometryUtilities::Polyhedron, 2> polyhedrons = {
      polyhedron_one,
      polyhedron_two
    };

    using merge_input_type = Gedim::GeometryUtilities::MergePolyhedronsInput;

    merge_input_type merged_polyhedrons_input;

    merged_polyhedrons_input.Vertices_Type[0] =
    {
      { merge_input_type::MergeTypes::Common, 0 },
      { merge_input_type::MergeTypes::Common, 1 },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Common, 2 },
      { merge_input_type::MergeTypes::Common, 3 },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none }
    };
    merged_polyhedrons_input.Vertices_Type[1] =
    {
      { merge_input_type::MergeTypes::Common, 0 },
      { merge_input_type::MergeTypes::Common, 1 },
      { merge_input_type::MergeTypes::Common, 2 },
      { merge_input_type::MergeTypes::Common, 3 }
    };
    merged_polyhedrons_input.Common_vertices =
    { { 0, 0 }, { 1, 1 }, { 3, 2 }, { 4, 3 } };

    merged_polyhedrons_input.Edges_Type[0] =
    {
      { merge_input_type::MergeTypes::Common, 0 },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Common, 1 },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Common, 2 },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none }
    };
    merged_polyhedrons_input.Edges_Type[1] =
    {
      { merge_input_type::MergeTypes::Common, 0 },
      { merge_input_type::MergeTypes::Common, 1 },
      { merge_input_type::MergeTypes::None, merge_input_type::none },
      { merge_input_type::MergeTypes::Common, 2 },
      { merge_input_type::MergeTypes::None, merge_input_type::none },
      { merge_input_type::MergeTypes::None, merge_input_type::none }
    };
    merged_polyhedrons_input.Common_edges =
    { { 0, 0 }, { 3, 1 }, { 8, 3 } };

    merged_polyhedrons_input.Faces_Type[0] =
    {
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none },
      { merge_input_type::MergeTypes::Remove, merge_input_type::none }
    };
    merged_polyhedrons_input.Faces_Type[1] =
    {
      { merge_input_type::MergeTypes::None, merge_input_type::none },
      { merge_input_type::MergeTypes::None, merge_input_type::none },
      { merge_input_type::MergeTypes::None, merge_input_type::none },
      { merge_input_type::MergeTypes::None, merge_input_type::none }
    };
    merged_polyhedrons_input.Common_faces = { };

    const auto merged_polyhedron = geometryUtilities.MergePolyhedrons(polyhedrons,
                                                                      merged_polyhedrons_input);
    std::string exportFolder = "./Export/TestGeometryUtilities/Test_MergePolyhedrons_edge_removed";
    Gedim::Output::CreateFolder(exportFolder);
    Gedim::Output::CreateFolder(exportFolder + "/polyhedron_one");
    geometryUtilities.ExportPolyhedronToVTU(polyhedron_one,
                                            exportFolder + "/polyhedron_one");
    Gedim::Output::CreateFolder(exportFolder + "/polyhedron_two");
    geometryUtilities.ExportPolyhedronToVTU(polyhedron_two,
                                            exportFolder + "/polyhedron_two");
    Gedim::Output::CreateFolder(exportFolder + "/merged_polyhedron");
    geometryUtilities.ExportPolyhedronToVTU(merged_polyhedron.MergedPolyhedron,
                                            exportFolder + "/merged_polyhedron");

    ASSERT_EQ(merged_polyhedron.MergedPolyhedron.Vertices.cols(),
              4);
    ASSERT_EQ(merged_polyhedron.MergedPolyhedron.Edges.cols(),
              6);
    ASSERT_EQ(merged_polyhedron.MergedPolyhedron.Faces.size(),
              4);

    std::array<std::vector<unsigned int>, 2> expected_original_to_merged_vertices;
    expected_original_to_merged_vertices[0] =
    {
      0u,
      1u,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      2u,
      3u,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
    };
    expected_original_to_merged_vertices[1] =
    {
      0u,
      1u,
      2u,
      3u
    };

    std::array<std::vector<unsigned int>, 2> expected_original_to_merged_edges;
    expected_original_to_merged_edges[0] =
    {
      0u,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      1u,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      2u,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none
    };
    expected_original_to_merged_edges[1] =
    {
      0u,
      1u,
      3u,
      2u,
      4u,
      5u
    };

    std::array<std::vector<unsigned int>, 2> expected_original_to_merged_faces;
    expected_original_to_merged_faces[0] =
    {
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none,
      Gedim::GeometryUtilities::MergePolyhedronsResult::none
    };
    expected_original_to_merged_faces[1] =
    {
      0u,
      1u,
      2u,
      3u
    };

    ASSERT_EQ(merged_polyhedron.OriginalToMergedVertices,
              expected_original_to_merged_vertices);
    ASSERT_EQ(merged_polyhedron.OriginalToMergedEdges,
              expected_original_to_merged_edges);
    ASSERT_EQ(merged_polyhedron.OriginalToMergedFaces,
              expected_original_to_merged_faces);

    const std::vector<std::array<unsigned int, 2>> expected_merged_to_original_vertices =
    {
      { 0u, 0u },
      { 1u, 1u },
      { 3u, 2u },
      { 4u, 3u }
    };
    const std::vector<std::array<unsigned int, 2>> expected_merged_to_original_edges =
    {
      { 0u, 0u },
      { 3u, 1u },
      { 8u, 3u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 2u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 4u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 5u }
    };
    const std::vector<std::array<unsigned int, 2>> expected_merged_to_original_faces =
    {
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 0u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 1u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 2u },
      { Gedim::GeometryUtilities::MergePolyhedronsResult::none, 3u }
    };

    ASSERT_EQ(merged_polyhedron.MergedToOriginalVertices,
              expected_merged_to_original_vertices);
    ASSERT_EQ(merged_polyhedron.MergedToOriginalEdges,
              expected_merged_to_original_edges);
    ASSERT_EQ(merged_polyhedron.MergedToOriginalFaces,
              expected_merged_to_original_faces);

  }
} // namespace GedimUnitTesting

#endif // __TEST_GEOMETRY_MERGED_POLYHEDRON_H
