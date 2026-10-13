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

#ifndef __TEST_meshUtilities3D_MarchingCubes_H
#define __TEST_meshUtilities3D_MarchingCubes_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "VtkMeshExporter.hpp"

namespace GedimUnitTesting
{

  TEST(TestMeshUtilities, TestMeshUtilities_3D_MarchingCubes)
  {
    // TPMS EQUATION

    // N.B. all variables (x,y,z) are intended as already scaled (x->x*2*Pi/L and so on).

    /*
  Gyroid:
  double val = std::sin(x) * std::cos(y) +
               std::sin(y) * std::cos(z) +
               std::sin(z) * std::cos(x)-c;

  Diamond:
  double val = std::sin(x) * std::sin(y) * std::sin(z) +
               std::sin(x) * std::cos(y) * std::cos(z) +
               std::cos(x) * std::sin(y) * std::cos(z) +
               std::cos(x) * std::cos(y) * std::sin(z) - c;


  Split-P:
  double val = 1.1*(std::sin(2*x) * std::sin(z) * std::cos(y) +
                    std::sin(2*y) * std::sin(x) * std::cos(z) +
                    std::sin(2*z)*std::sin(y)*std::cos(x)) -
               0.2*(std::cos(2*x) * std::cos(2*y) + std::cos(2*y) * std::cos(2*z) +
                    std::cos(2*z) * std::cos(2*x)) -
               0.4*(std::cos(2*x) * std::cos(2*y) * std::cos(2*z)) - c;


  Other (instable) TPMS:

  Lidinoid:
  double val = std::sin(2*x) * std::cos(y) * std::sin(z) +
               std::sin(2*y) * std::cos(z) * std::sin(x) +
               std::sin(2*z) * std::cos(x) * std::sin(y) -
               std::cos(2*x) * std::cos(2*y) -
               std::cos(2*y) * std::cos(2*z) -
               std::cos(2*z) * std::cos(2*x) + 0.3;

  Neovius:
  double val = 3.0*(std::cos(x) +std::cos(y) +std::cos(z)) +
               4 * std::cos(x) *std::cos(y) *std::cos(z);


  */

    std::string exportFolder = "./Export/TestMeshUtilities/TestMeshUtilities_3D_MarchingCubes";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
      const Gedim::GeometryUtilities::Polyhedron cube =
          geometry_utilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
      const Eigen::Vector3d parallelepipedsLengthTangent = cube.Vertices.col(1) - cube.Vertices.col(0);
      const Eigen::Vector3d parallelepipedsHeightTangent = cube.Vertices.col(3) - cube.Vertices.col(0);
      const Eigen::Vector3d parallelepipedsWidthTangent = cube.Vertices.col(4) - cube.Vertices.col(0);

      const std::vector<double> lengthMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4, 0.0, 1.0, 1);
      const std::vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(5, 0.0, 1.0, 1);
      const std::vector<double> widthMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(6, 0.0, 1.0, 1);

      Gedim::MeshUtilities::CreateParallelepipedMesh(cube.Vertices.col(0),
                                                     parallelepipedsLengthTangent,
                                                     parallelepipedsHeightTangent,
                                                     parallelepipedsWidthTangent,
                                                     lengthMeshCurvilinearCoordinates,
                                                     heightMeshCurvilinearCoordinates,
                                                     widthMeshCurvilinearCoordinates,
                                                     mesh);
    }

    Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, exportFolder, "original_mesh");

    auto surface_implicit_equation = [](const Eigen::MatrixXd &points) {
      float L = 1.0f; // splitp = 5.0f, nevious = 3.0f.
      float scale = (2.0 * std::numbers::pi) / L;
      float c = 0.0f;

      Eigen::VectorXd surface_values(points.cols());
      for (unsigned int p = 0; p < surface_values.size(); p++)
      {
        const Eigen::Vector3d point = points.col(p);
        double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
        const double suface_value = std::sin(x) * std::sin(y) * std::sin(z) + std::sin(x) * std::cos(y) * std::cos(z) +
                                    std::cos(x) * std::sin(y) * std::cos(z) + std::cos(x) * std::cos(y) * std::sin(z) - c;

        surface_values[p] = suface_value;
      }

      return surface_values;
    };

    const auto marching_cubes_result = Gedim::MeshUtilities::MarchingCubes(geometry_utilities, surface_implicit_equation, mesh);

    {
      //        {
      //            Gedim::VTKUtilities exporter;

      //            auto random_points = geometry_utilities.generate_uniform_random_points_in_sphere(100000,
      //                                                                                             0.5 * std::sqrt(3.0));
      //            random_points.colwise() += Eigen::Vector3d(0.5, 0.5, 0.5);
      //            const auto surface_values = surface_implicit_equation(random_points);

      //            for (unsigned int c = 0; c < random_points.cols(); ++c)
      //            {
      //                const double surface_value = surface_values[c];

      //                if (!geometry_utilities.IsValueLowerOrEqual(std::abs(surface_value),
      //                                                            1.0e-2,
      //                                                            geometry_utilities.Tolerance1D()))
      //                    continue;

      //                exporter.AddPoint(random_points.col(c),
      //                                  {{"surface_value", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &surface_value}});
      //            }

      //            exporter.Export(exportFolder + "/surface.vtu");
      //        }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
        {
          const double cell0D_value = marching_cubes_result.cell0Ds_surface_value[c];
          const auto cell0D_mark = marching_cubes_result.cell0Ds_surface_position.at(c);

          double cell0D_sign = 0.0;
          switch (cell0D_mark)
          {
            case Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::OnSurface:
              cell0D_sign = 0.0;
              break;
            case Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::Positive:
              cell0D_sign = +1.0;
              break;
            case Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::Negative:
              cell0D_sign = -1.0;
              break;
            default:
              throw std::runtime_error("Position not managed");
          }

          exporter.AddPoint(mesh.Cell0DCoordinates(c),
                            {
                              {
                                "sign",
                                Gedim::VTPProperty::Formats::Cells,
                                static_cast<unsigned int>(1),
                                &cell0D_sign
                              },
                              {
                                "value",
                                Gedim::VTPProperty::Formats::Cells,
                                static_cast<unsigned int>(1),
                                &cell0D_value
                              }
                            }
                            );
        }

        exporter.Export(exportFolder + "/cell0Ds_surface_position.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int p = 0; p < marching_cubes_result.surface_intersections.surface_intersection_points.size(); ++p)
        {
          const auto& surface_point = marching_cubes_result.surface_intersections.surface_intersection_points.at(p);
          const double surface_id = surface_point.mesh_id;

          double surface_type = 0.0;
          switch (surface_point.Type)
          {
            case Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex:
              surface_type = 1.0;
              break;
            case Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Edge:
              surface_type = 2.0;
              break;
            default:
              throw std::runtime_error("Surface Type not managed");
          }

          exporter.AddPoint(surface_point.coordinate,
                            {
                              {
                                "type",
                                Gedim::VTPProperty::Formats::Cells,
                                static_cast<unsigned int>(1),
                                &surface_type
                              },
                              {
                                "id",
                                Gedim::VTPProperty::Formats::Cells,
                                static_cast<unsigned int>(1),
                                &surface_id
                              }
                            });
        }

        exporter.Export(exportFolder + "/marching_cubes_points.vtu");
      }

      {
        Gedim::VTKUtilities exporter;
;
        for (unsigned int p = 0; p < marching_cubes_result.surface_triangles.size(); ++p)
        {
          const auto& triangle_ids = marching_cubes_result.surface_triangles.at(p);
          Eigen::Matrix3d triangle;
          triangle.col(0) = marching_cubes_result.surface_intersections.surface_intersection_points.at(triangle_ids.at(0)).coordinate;
          triangle.col(1) = marching_cubes_result.surface_intersections.surface_intersection_points.at(triangle_ids.at(1)).coordinate;
          triangle.col(2) = marching_cubes_result.surface_intersections.surface_intersection_points.at(triangle_ids.at(2)).coordinate;

          exporter.AddPolygon(triangle);
        }

        exporter.Export(exportFolder + "/marching_cubes_triangles.vtu");
      }
    }
  }

} // namespace GedimUnitTesting

#endif // __TEST_meshUtilities3D_MarchingCubes_H
