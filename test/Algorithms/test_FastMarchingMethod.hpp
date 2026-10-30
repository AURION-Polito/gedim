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

#ifndef __TEST_FastMarchingMethod_H
#define __TEST_FastMarchingMethod_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "FastMarchingMethod.hpp"
#include "MeshUtilities.hpp"
#include "TriangleInterface.hpp"
#include "VTKUtilities.hpp"
#include "VtkMeshExporter.hpp"

namespace UnitTesting
{
  TEST(TestFastMarchingMethod, TestFastMarchingMethod_Simple)
  {
    std::string exportFolder = "./Export/TestFastMarchingMethod/TestFastMarchingMethod_Simple";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
      const auto square = geometry_utilities.CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);

      Gedim::External::MeshUtilities::CreateTriangularMesh(square,
                                                           0.01,
                                                           mesh);
      Gedim::MeshUtilities::ComputeCell0DCell1DNeighbours(mesh);
    }

    {
      const std::string mesh_folder = exportFolder + "/mesh";
      Gedim::Output::CreateFolder(mesh_folder);
      Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_folder, "mesh");
    }

    const double speed = 1.0;
    const unsigned int source = mesh.Cell0DTotalNumber() > 3 ? 4 : 0;
    const auto result = Gedim::Algorithms::FastMarchingMethod::fast_marching_method_2D(geometry_utilities,
                                                                                       mesh,
                                                                                       speed,
                                                                                       source);

    {
      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
        {
          double cell0D_status = 0.0;
          switch (result.vertices_status.at(c))
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes::Far:
              cell0D_status = -1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes::Trial:
              cell0D_status = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes::Verified:
              cell0D_status = +1.0;
              break;
            default:
              throw std::runtime_error("status not managed");
          }

          exporter.AddPoint(
                mesh.Cell0DCoordinates(c),
                {{"status", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_status},
                 {"time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &result.vertices_time.at(c)}});
        }

        exporter.Export(exportFolder + "/cell0Ds_time.vtu");
      }

    }
  }

} // namespace UnitTesting

#endif // __TEST_FastMarchingMethod_H
