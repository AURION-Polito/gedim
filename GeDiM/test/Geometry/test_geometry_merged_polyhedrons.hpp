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

  TEST(TestGeometryUtilities, Test_MergePolyhedrons)
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
    std::string exportFolder = "./Export/TestGeometryUtilities/Test_MergePolyhedrons";
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

} // namespace GedimUnitTesting

#endif // __TEST_GEOMETRY_MERGED_POLYHEDRON_H
