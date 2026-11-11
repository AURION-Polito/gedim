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
#include "TetgenInterface.hpp"
#include "TriangleInterface.hpp"
#include "VTKUtilities.hpp"
#include "VtkMeshExporter.hpp"

namespace UnitTesting
{
  TEST(TestFastMarchingMethod, TestFastMarchingMethod_2D_Simple)
  {
    std::string exportFolder = "./Export/TestFastMarchingMethod/TestFastMarchingMethod_2D_Simple";
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
    const auto fmm_mesh_graph = Gedim::Algorithms::FastMarchingMethod::create_mesh_graph(mesh);
    const auto fmm_config = Gedim::Algorithms::FastMarchingMethod::create_configuration(geometry_utilities,
                                                                                        fmm_mesh_graph,
                                                                                        { source },
                                                                                        0.0,
                                                                                        speed,
                                                                                        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes::Linear);
    const auto result = Gedim::Algorithms::FastMarchingMethod::fast_marching_method(geometry_utilities,
                                                                                    fmm_mesh_graph,
                                                                                    fmm_config);

    {
      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
        {
          const auto& vertex_config = fmm_config.cell0Ds_initial_condition.at(c);

          double cell0D_init_type = 0.0;
          switch (vertex_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Source:
              cell0D_init_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Initialized:
              cell0D_init_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unreachable:
              cell0D_init_type = std::numeric_limits<double>::infinity();
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unknown:
              cell0D_init_type = -1.0;
              break;
            default:
              throw std::runtime_error("status not managed");
          }

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
                {
                  {"init_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_init_type},
                  {"init_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.initial_time},
                  {"border_distance", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.border_distance},
                  {"status", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_status},
                  {"final_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &result.vertices_time.at(c)}
                });
        }

        exporter.Export(exportFolder + "/cell0Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
        {
          const auto& edge_config = fmm_config.cell1Ds_restriction.at(e);

          double cell1D_restriction_type = 0.0;
          switch (edge_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Free:
              cell1D_restriction_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Constrained:
              cell1D_restriction_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Unreachable:
              cell1D_restriction_type = std::numeric_limits<double>::infinity();
              break;
            default:
              throw std::runtime_error("restriction not managed");
          }

          exporter.AddSegment(
                mesh.Cell1DOriginCoordinates(e),
                mesh.Cell1DEndCoordinates(e),
                {
                  {"restriction_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell1D_restriction_type},
                });
        }

        exporter.Export(exportFolder + "/cell1Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                             mesh.Cell2DsVertices(),
                             {{"time", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(result.vertices_time.size()), result.vertices_time.data()}});

        exporter.Export(exportFolder + "/cell2Ds_time.vtu");
      }

    }
  }

  TEST(TestFastMarchingMethod, TestFastMarchingMethod_3D_Simple)
  {
    std::string exportFolder = "./Export/TestFastMarchingMethod/TestFastMarchingMethod_3D_Simple";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
      const auto cube = geometry_utilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);

      Gedim::External::MeshUtilities::CreateTetrahedralMesh(cube.Vertices,
                                                            cube.Edges,
                                                            cube.Faces,
                                                            0.001,
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
    const auto fmm_mesh_graph = Gedim::Algorithms::FastMarchingMethod::create_mesh_graph(mesh);
    const auto fmm_config = Gedim::Algorithms::FastMarchingMethod::create_configuration(geometry_utilities,
                                                                                        fmm_mesh_graph,
                                                                                        { source },
                                                                                        0.0,
                                                                                        speed,
                                                                                        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes::Linear);
    const auto result = Gedim::Algorithms::FastMarchingMethod::fast_marching_method(geometry_utilities,
                                                                                    fmm_mesh_graph,
                                                                                    fmm_config);

    {
      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
        {
          const auto& vertex_config = fmm_config.cell0Ds_initial_condition.at(c);

          double cell0D_init_type = 0.0;
          switch (vertex_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Source:
              cell0D_init_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Initialized:
              cell0D_init_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unreachable:
              cell0D_init_type = std::numeric_limits<double>::infinity();
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unknown:
              cell0D_init_type = -1.0;
              break;
            default:
              throw std::runtime_error("status not managed");
          }

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
                {
                  {"init_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_init_type},
                  {"init_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.initial_time},
                  {"border_distance", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.border_distance},
                  {"status", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_status},
                  {"final_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &result.vertices_time.at(c)}
                });
        }

        exporter.Export(exportFolder + "/cell0Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
        {
          const auto& edge_config = fmm_config.cell1Ds_restriction.at(e);

          double cell1D_restriction_type = 0.0;
          switch (edge_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Free:
              cell1D_restriction_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Constrained:
              cell1D_restriction_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Unreachable:
              cell1D_restriction_type = std::numeric_limits<double>::infinity();
              break;
            default:
              throw std::runtime_error("restriction not managed");
          }

          exporter.AddSegment(
                mesh.Cell1DOriginCoordinates(e),
                mesh.Cell1DEndCoordinates(e),
                {
                  {"restriction_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell1D_restriction_type},
                });
        }

        exporter.Export(exportFolder + "/cell1Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                             mesh.Cell2DsVertices(),
                             {{"time", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(result.vertices_time.size()), result.vertices_time.data()}});

        exporter.Export(exportFolder + "/cell2Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPolyhedrons(mesh.Cell0DsCoordinates(),
                             mesh.Cell3DsFacesVertices(),
                             {{"time", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(result.vertices_time.size()), result.vertices_time.data()}});

        exporter.Export(exportFolder + "/cell3Ds_time.vtu");
      }

    }
  }

  TEST(TestFastMarchingMethod, TestFastMarchingMethod_Simple_Circle)
  {
    std::string exportFolder = "./Export/TestFastMarchingMethod/TestFastMarchingMethod_Simple_Circle";
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

    const double R = 0.01;
    const double speed = 1.0;
    const unsigned int source = mesh.Cell0DTotalNumber() > 3 ? 4 : 0;
    const auto fmm_mesh_graph = Gedim::Algorithms::FastMarchingMethod::create_mesh_graph(mesh);
    const auto fmm_config = Gedim::Algorithms::FastMarchingMethod::create_configuration(geometry_utilities,
                                                                                        fmm_mesh_graph,
                                                                                        { source },
                                                                                        R,
                                                                                        speed,
                                                                                        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes::Linear);
    const auto result = Gedim::Algorithms::FastMarchingMethod::fast_marching_method(geometry_utilities,
                                                                                    fmm_mesh_graph,
                                                                                    fmm_config);


    {
      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
        {
          const auto& vertex_config = fmm_config.cell0Ds_initial_condition.at(c);

          double cell0D_init_type = 0.0;
          switch (vertex_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Source:
              cell0D_init_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Initialized:
              cell0D_init_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unreachable:
              cell0D_init_type = std::numeric_limits<double>::infinity();
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unknown:
              cell0D_init_type = -1.0;
              break;
            default:
              throw std::runtime_error("status not managed");
          }

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
                {
                  {"init_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_init_type},
                  {"init_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.initial_time},
                  {"border_distance", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.border_distance},
                  {"status", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_status},
                  {"final_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &result.vertices_time.at(c)}
                });
        }

        exporter.Export(exportFolder + "/cell0Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
        {
          const auto& edge_config = fmm_config.cell1Ds_restriction.at(e);

          double cell1D_restriction_type = 0.0;
          switch (edge_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Free:
              cell1D_restriction_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Constrained:
              cell1D_restriction_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Unreachable:
              cell1D_restriction_type = std::numeric_limits<double>::infinity();
              break;
            default:
              throw std::runtime_error("restriction not managed");
          }

          exporter.AddSegment(
                mesh.Cell1DOriginCoordinates(e),
                mesh.Cell1DEndCoordinates(e),
                {
                  {"restriction_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell1D_restriction_type},
                });
        }

        exporter.Export(exportFolder + "/cell1Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                             mesh.Cell2DsVertices(),
                             {{"time", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(result.vertices_time.size()), result.vertices_time.data()}});

        exporter.Export(exportFolder + "/cell2Ds_time.vtu");
      }

    }
  }

  TEST(TestFastMarchingMethod, TestFastMarchingMethod_Curves)
  {
    std::string exportFolder = "./Export/TestFastMarchingMethod/TestFastMarchingMethod_Curves";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
      Gedim::MeshMatrices original_mesh_data;
      Gedim::MeshMatricesDAO original_mesh(original_mesh_data);

      {
        const double domain_area = 1.0;
        const double max_relative_area = 0.001;
        const double max_cell_edge = std::sqrt(domain_area * max_relative_area);
        const unsigned int num_cells_base = std::ceil(1.0 / max_cell_edge) + 1;
        const unsigned int num_cells_height = std::ceil(1.0 / max_cell_edge) + 1;

        const auto square = geometry_utilities.CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
        const Eigen::Vector3d square_LengthTangent = square.col(1) - square.col(0);
        const Eigen::Vector3d square_HeightTangent = square.col(3) - square.col(0);

        const auto lengthMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(num_cells_base, 0.0, 1.0, 1);
        const auto heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(num_cells_height, 0.0, 1.0, 1);

        Gedim::MeshUtilities::CreateRectangleMesh(square.col(0),
                                                  square_LengthTangent,
                                                  square_HeightTangent,
                                                  lengthMeshCurvilinearCoordinates,
                                                  heightMeshCurvilinearCoordinates,
                                                  original_mesh);
      }

      auto plane_x = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const Eigen::Vector3d p0 = Eigen::Vector3d(0.4, 0.4, 0.0);
        const Eigen::Vector3d n = Eigen::Vector3d(0.0, 1.0, 0.0).normalized();
        return (point - p0).dot(n);
      };

      auto generic_plane = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const Eigen::Vector3d p0 = Eigen::Vector3d(0.6, 0.4, 0.0);
        const Eigen::Vector3d n = Eigen::Vector3d(0.5, 1.0, 0.0).normalized();
        return (point - p0).dot(n);
      };

      auto circle = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = 1.0;
        const double radius = 0.2;
        const double x = scale * (point.x() - c - 0.5), y = scale * (point.y() - c - 0.5);
        return x * x + y * y - radius * radius;
      };

      auto gyroid = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = (2.0 * std::numbers::pi) / L;
        const double x = scale * point.x(), y = scale * point.y();
        return std::sin(x) * std::cos(y) + std::sin(y) - c * scale;
      };

      auto circle_shifted = [&circle](const Eigen::Vector3d &point, const double &L, const double &c) {
        return circle(point, L, c) * circle(point, L, c + 0.1);
      };

      auto gyroid_shifted = [&gyroid](const Eigen::Vector3d &point, const double &L, const double &c) {
        return gyroid(point, L, c + 0.1) * gyroid(point, L, c + 0.0);
      };

      auto surface_implicit_equation = [&plane_x, &generic_plane, &circle, &gyroid, &gyroid_shifted, &circle_shifted](const Eigen::MatrixXd &points) {
        Eigen::VectorXd surface_values(points.cols());
        for (unsigned int p = 0; p < surface_values.size(); p++)
          surface_values[p] = circle_shifted(points.col(p), 1.0, 0.0);

        return surface_values;
      };

      const auto marching_squares_result =
          Gedim::MeshUtilities::MarchingSquares(geometry_utilities,
                                                surface_implicit_equation,
                                                original_mesh);

      const auto surface_intersection_graph =
          Gedim::MeshUtilities::ComputeMeshSurfaceIntersectionGraph(marching_squares_result.surface_intersections,
                                                                    marching_squares_result.surface_segments,
                                                                    {},
                                                                    original_mesh);

      Gedim::MeshMatrices conformed_mesh_data;
      Gedim::MeshMatricesDAO conformed_mesh(conformed_mesh_data);

      const auto conform_result = Gedim::MeshUtilities::ConformMeshToSurface(original_mesh,
                                                                             surface_intersection_graph,
                                                                             conformed_mesh);

      std::vector<unsigned int> cell2D_to_mantain;
      cell2D_to_mantain.reserve(conformed_mesh.Cell2DTotalNumber());
      for (unsigned int c = 0; c < conformed_mesh.Cell2DTotalNumber(); ++c)
      {
        switch (conform_result.new_cell2Ds_surface_position.at(c))
        {
          case Gedim::MeshUtilities::SurfacePositions::OnSurface:
          case Gedim::MeshUtilities::SurfacePositions::Positive:
            cell2D_to_mantain.push_back(c);
            break;
          case Gedim::MeshUtilities::SurfacePositions::Negative:
            break;
          default:
            throw std::runtime_error("Position not managed");
        }
      }
      cell2D_to_mantain.shrink_to_fit();

      const auto mesh_filter = Gedim::MeshUtilities::FilterMesh2D(cell2D_to_mantain,
                                                                  conformed_mesh);

      Gedim::MeshMatrices filtered_mesh_data;
      Gedim::MeshMatricesDAO filtered_mesh(filtered_mesh_data);

      Gedim::MeshUtilities::ExtractMesh2D(mesh_filter.Cell0Ds,
                                          mesh_filter.Cell1Ds,
                                          mesh_filter.Cell2Ds,
                                          conformed_mesh,
                                          filtered_mesh);

      const auto cell2Ds_triangles = Gedim::MeshUtilities::ComputeMeshCell2DsTriangles(geometry_utilities,
                                                                                       Gedim::MeshUtilities::MeshCell2DsTrianglesData::TrianglesTypes::ByBarycenter,
                                                                                       filtered_mesh);

      const auto extract_tetra = Gedim::MeshUtilities::ExtractTriangularMesh(cell2Ds_triangles,
                                                                             filtered_mesh);

      Gedim::MeshUtilities::TetrahedralMeshDataToMeshDAO(extract_tetra,
                                                         filtered_mesh,
                                                         mesh);
      Gedim::MeshUtilities::ComputeCell0DCell1DNeighbours(mesh);
    }

    {
      const std::string mesh_folder = exportFolder + "/mesh";
      Gedim::Output::CreateFolder(mesh_folder);
      Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_folder, "mesh");
    }

    const double speed = 1.0;
    const unsigned int source = mesh.Cell0DTotalNumber() > 170 ? 170 : 0;
    const auto fmm_mesh_graph = Gedim::Algorithms::FastMarchingMethod::create_mesh_graph(mesh);
    const auto fmm_config = Gedim::Algorithms::FastMarchingMethod::create_configuration(geometry_utilities,
                                                                                        fmm_mesh_graph,
                                                                                        { source },
                                                                                        0.0,
                                                                                        speed,
                                                                                        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes::Linear);
    const auto result = Gedim::Algorithms::FastMarchingMethod::fast_marching_method(geometry_utilities,
                                                                                    fmm_mesh_graph,
                                                                                    fmm_config);
    {
      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
        {
          const auto& vertex_config = fmm_config.cell0Ds_initial_condition.at(c);

          double cell0D_init_type = 0.0;
          switch (vertex_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Source:
              cell0D_init_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Initialized:
              cell0D_init_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unreachable:
              cell0D_init_type = std::numeric_limits<double>::infinity();
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unknown:
              cell0D_init_type = -1.0;
              break;
            default:
              throw std::runtime_error("status not managed");
          }

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
                {
                  {"init_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_init_type},
                  {"init_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.initial_time},
                  {"border_distance", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.border_distance},
                  {"status", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_status},
                  {"final_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &result.vertices_time.at(c)}
                });
        }

        exporter.Export(exportFolder + "/cell0Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
        {
          const auto& edge_config = fmm_config.cell1Ds_restriction.at(e);

          double cell1D_restriction_type = 0.0;
          switch (edge_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Free:
              cell1D_restriction_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Constrained:
              cell1D_restriction_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Unreachable:
              cell1D_restriction_type = std::numeric_limits<double>::infinity();
              break;
            default:
              throw std::runtime_error("restriction not managed");
          }

          exporter.AddSegment(
                mesh.Cell1DOriginCoordinates(e),
                mesh.Cell1DEndCoordinates(e),
                {
                  {"restriction_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell1D_restriction_type},
                });
        }

        exporter.Export(exportFolder + "/cell1Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                             mesh.Cell2DsVertices(),
                             {{"time", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(result.vertices_time.size()), result.vertices_time.data()}});

        exporter.Export(exportFolder + "/cell2Ds_time.vtu");
      }

    }
  }

  TEST(TestFastMarchingMethod, TestFastMarchingMethod_Curves_Circle)
  {
    std::string exportFolder = "./Export/TestFastMarchingMethod/TestFastMarchingMethod_Curves_Circle";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
      Gedim::MeshMatrices original_mesh_data;
      Gedim::MeshMatricesDAO original_mesh(original_mesh_data);

      {
        const double domain_area = 1.0;
        const double max_relative_area = 0.001;
        const double max_cell_edge = std::sqrt(domain_area * max_relative_area);
        const unsigned int num_cells_base = std::ceil(1.0 / max_cell_edge) + 1;
        const unsigned int num_cells_height = std::ceil(1.0 / max_cell_edge) + 1;

        const auto square = geometry_utilities.CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
        const Eigen::Vector3d square_LengthTangent = square.col(1) - square.col(0);
        const Eigen::Vector3d square_HeightTangent = square.col(3) - square.col(0);

        const auto lengthMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(num_cells_base, 0.0, 1.0, 1);
        const auto heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(num_cells_height, 0.0, 1.0, 1);

        Gedim::MeshUtilities::CreateRectangleMesh(square.col(0),
                                                  square_LengthTangent,
                                                  square_HeightTangent,
                                                  lengthMeshCurvilinearCoordinates,
                                                  heightMeshCurvilinearCoordinates,
                                                  original_mesh);
      }

      auto plane_x = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const Eigen::Vector3d p0 = Eigen::Vector3d(0.4, 0.4, 0.0);
        const Eigen::Vector3d n = Eigen::Vector3d(0.0, 1.0, 0.0).normalized();
        return (point - p0).dot(n);
      };

      auto generic_plane = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const Eigen::Vector3d p0 = Eigen::Vector3d(0.6, 0.4, 0.0);
        const Eigen::Vector3d n = Eigen::Vector3d(0.5, 1.0, 0.0).normalized();
        return (point - p0).dot(n);
      };

      auto circle = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = 1.0;
        const double radius = 0.2;
        const double x = scale * (point.x() - c - 0.5), y = scale * (point.y() - c - 0.5);
        return x * x + y * y - radius * radius;
      };

      auto gyroid = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = (2.0 * std::numbers::pi) / L;
        const double x = scale * point.x(), y = scale * point.y();
        return std::sin(x) * std::cos(y) + std::sin(y) - c * scale;
      };

      auto circle_shifted = [&circle](const Eigen::Vector3d &point, const double &L, const double &c) {
        return circle(point, L, c) * circle(point, L, c + 0.1);
      };

      auto gyroid_shifted = [&gyroid](const Eigen::Vector3d &point, const double &L, const double &c) {
        return gyroid(point, L, c + 0.1) * gyroid(point, L, c + 0.0);
      };

      auto surface_implicit_equation = [&plane_x, &generic_plane, &circle, &gyroid, &gyroid_shifted, &circle_shifted](const Eigen::MatrixXd &points) {
        Eigen::VectorXd surface_values(points.cols());
        for (unsigned int p = 0; p < surface_values.size(); p++)
          surface_values[p] = circle_shifted(points.col(p), 1.0, 0.0);

        return surface_values;
      };

      unsigned int original_mesh_max_marker = 0;
      for (unsigned int v = 0; v < original_mesh.Cell0DTotalNumber(); ++v)
      {
        if (!original_mesh.Cell0DIsActive(v))
          continue;

        if (original_mesh.Cell0DMarker(v) > original_mesh_max_marker)
          original_mesh_max_marker = original_mesh.Cell0DMarker(v);
      }
      for (unsigned int e = 0; e < original_mesh.Cell1DTotalNumber(); ++e)
      {
        if (!original_mesh.Cell1DIsActive(e))
          continue;

        if (original_mesh.Cell1DMarker(e) > original_mesh_max_marker)
          original_mesh_max_marker = original_mesh.Cell1DMarker(e);
      }

      const auto marching_squares_result =
          Gedim::MeshUtilities::MarchingSquares(geometry_utilities,
                                                surface_implicit_equation,
                                                original_mesh);

      const auto surface_intersection_graph =
          Gedim::MeshUtilities::ComputeMeshSurfaceIntersectionGraph(marching_squares_result.surface_intersections,
                                                                    marching_squares_result.surface_segments,
                                                                    {},
                                                                    original_mesh);

      Gedim::MeshMatrices conformed_mesh_data;
      Gedim::MeshMatricesDAO conformed_mesh(conformed_mesh_data);

      const auto conform_result = Gedim::MeshUtilities::ConformMeshToSurface(original_mesh,
                                                                             surface_intersection_graph,
                                                                             conformed_mesh);

      for (unsigned int v = 0; v < conformed_mesh.Cell0DTotalNumber(); ++v)
      {
        if (!conformed_mesh.Cell0DIsActive(v))
          continue;

        if (conformed_mesh.Cell0DMarker(v) != 0)
          continue;

        switch (conform_result.new_cell0Ds_surface_position.at(v))
        {
          case Gedim::MeshUtilities::SurfacePositions::OnSurface:
            conformed_mesh.Cell0DSetMarker(v, original_mesh_max_marker + 1);
            break;
          case Gedim::MeshUtilities::SurfacePositions::Positive:
            break;
          case Gedim::MeshUtilities::SurfacePositions::Negative:
            break;
          default:
            throw std::runtime_error("Position not managed");
        }
      }

      for (unsigned int e = 0; e < conformed_mesh.Cell1DTotalNumber(); ++e)
      {
        if (!conformed_mesh.Cell1DIsActive(e))
          continue;

        if (conformed_mesh.Cell1DMarker(e) != 0)
          continue;

        switch (conform_result.new_cell1Ds_surface_position.at(e))
        {
          case Gedim::MeshUtilities::SurfacePositions::OnSurface:
            conformed_mesh.Cell1DSetMarker(e, original_mesh_max_marker + 1);
            break;
          case Gedim::MeshUtilities::SurfacePositions::Positive:
            break;
          case Gedim::MeshUtilities::SurfacePositions::Negative:
            break;
          default:
            throw std::runtime_error("Position not managed");
        }
      }

      std::vector<unsigned int> cell2D_to_mantain;
      cell2D_to_mantain.reserve(conformed_mesh.Cell2DTotalNumber());
      for (unsigned int c = 0; c < conformed_mesh.Cell2DTotalNumber(); ++c)
      {
        switch (conform_result.new_cell2Ds_surface_position.at(c))
        {
          case Gedim::MeshUtilities::SurfacePositions::OnSurface:
          case Gedim::MeshUtilities::SurfacePositions::Positive:
            cell2D_to_mantain.push_back(c);
            break;
          case Gedim::MeshUtilities::SurfacePositions::Negative:
            break;
          default:
            throw std::runtime_error("Position not managed");
        }
      }
      cell2D_to_mantain.shrink_to_fit();

      const auto mesh_filter = Gedim::MeshUtilities::FilterMesh2D(cell2D_to_mantain,
                                                                  conformed_mesh);

      Gedim::MeshMatrices filtered_mesh_data;
      Gedim::MeshMatricesDAO filtered_mesh(filtered_mesh_data);

      Gedim::MeshUtilities::ExtractMesh2D(mesh_filter.Cell0Ds,
                                          mesh_filter.Cell1Ds,
                                          mesh_filter.Cell2Ds,
                                          conformed_mesh,
                                          filtered_mesh);

      const auto cell2Ds_triangles = Gedim::MeshUtilities::ComputeMeshCell2DsTriangles(geometry_utilities,
                                                                                       Gedim::MeshUtilities::MeshCell2DsTrianglesData::TrianglesTypes::ByBarycenter,
                                                                                       filtered_mesh);

      const auto extract_tetra = Gedim::MeshUtilities::ExtractTriangularMesh(cell2Ds_triangles,
                                                                             filtered_mesh);

      Gedim::MeshUtilities::TetrahedralMeshDataToMeshDAO(extract_tetra,
                                                         filtered_mesh,
                                                         mesh);
      Gedim::MeshUtilities::ComputeCell0DCell1DNeighbours(mesh);
    }

    {
      const std::string mesh_folder = exportFolder + "/mesh";
      Gedim::Output::CreateFolder(mesh_folder);
      Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_folder, "mesh");
    }

    const double R = 0.05;
    const double speed = 1.0;
    const unsigned int source = mesh.Cell0DTotalNumber() > 3 ? 4 : 0;
    const auto fmm_mesh_graph = Gedim::Algorithms::FastMarchingMethod::create_mesh_graph(mesh);
    const auto fmm_config = Gedim::Algorithms::FastMarchingMethod::create_configuration(geometry_utilities,
                                                                                        fmm_mesh_graph,
                                                                                        { source },
                                                                                        R,
                                                                                        speed,
                                                                                        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes::Linear);
    const auto result = Gedim::Algorithms::FastMarchingMethod::fast_marching_method(geometry_utilities,
                                                                                    fmm_mesh_graph,
                                                                                    fmm_config);


    {
      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
        {
          const auto& vertex_config = fmm_config.cell0Ds_initial_condition.at(c);

          double cell0D_init_type = 0.0;
          switch (vertex_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Source:
              cell0D_init_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Initialized:
              cell0D_init_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unreachable:
              cell0D_init_type = std::numeric_limits<double>::infinity();
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unknown:
              cell0D_init_type = -1.0;
              break;
            default:
              throw std::runtime_error("status not managed");
          }

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
                {
                  {"init_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_init_type},
                  {"init_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.initial_time},
                  {"border_distance", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.border_distance},
                  {"status", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_status},
                  {"final_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &result.vertices_time.at(c)}
                });
        }

        exporter.Export(exportFolder + "/cell0Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
        {
          const auto& edge_config = fmm_config.cell1Ds_restriction.at(e);

          double cell1D_restriction_type = 0.0;
          switch (edge_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Free:
              cell1D_restriction_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Constrained:
              cell1D_restriction_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Unreachable:
              cell1D_restriction_type = std::numeric_limits<double>::infinity();
              break;
            default:
              throw std::runtime_error("restriction not managed");
          }

          exporter.AddSegment(
                mesh.Cell1DOriginCoordinates(e),
                mesh.Cell1DEndCoordinates(e),
                {
                  {"restriction_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell1D_restriction_type},
                });
        }

        exporter.Export(exportFolder + "/cell1Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                             mesh.Cell2DsVertices(),
                             {{"time", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(result.vertices_time.size()), result.vertices_time.data()}});

        exporter.Export(exportFolder + "/cell2Ds_time.vtu");
      }
    }
  }

  TEST(TestFastMarchingMethod, TestFastMarchingMethod_3D_Curves)
  {
    std::string exportFolder = "./Export/TestFastMarchingMethod/TestFastMarchingMethod_3D_Curves";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
      Gedim::MeshMatrices original_mesh_data;
      Gedim::MeshMatricesDAO original_mesh(original_mesh_data);

      {
        const double domain_volume = 1.0;
        const double max_relative_volume = 0.001;
        const double max_cell_edge = std::cbrt(domain_volume * max_relative_volume);

        const Gedim::GeometryUtilities::Polyhedron cube =
            geometry_utilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
        const Eigen::Vector3d parallelepipedsLengthTangent = cube.Vertices.col(1) - cube.Vertices.col(0);
        const Eigen::Vector3d parallelepipedsHeightTangent = cube.Vertices.col(3) - cube.Vertices.col(0);
        const Eigen::Vector3d parallelepipedsWidthTangent = cube.Vertices.col(4) - cube.Vertices.col(0);

        const std::vector<double> lengthMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(max_cell_edge, 0.0, 1.0, 1);
        const std::vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(max_cell_edge, 0.0, 1.0, 1);
        const std::vector<double> widthMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(max_cell_edge, 0.0, 1.0, 1);

        Gedim::MeshUtilities::CreateParallelepipedMesh(cube.Vertices.col(0),
                                                       parallelepipedsLengthTangent,
                                                       parallelepipedsHeightTangent,
                                                       parallelepipedsWidthTangent,
                                                       lengthMeshCurvilinearCoordinates,
                                                       heightMeshCurvilinearCoordinates,
                                                       widthMeshCurvilinearCoordinates,
                                                       mesh);
      }

      auto plane_x = [](const Eigen::Vector3d &point, const double &L, const double &c) {
          const Eigen::Vector3d p0 = Eigen::Vector3d(0.4, 0.4, 0.4);
          const Eigen::Vector3d n = Eigen::Vector3d(0.0, 1.0, 0.0).normalized();
          return (point - p0).dot(n);
      };

      auto generic_plane = [](const Eigen::Vector3d &point, const double &L, const double &c) {
          const Eigen::Vector3d p0 = Eigen::Vector3d(0.6, 0.4, 1.0);
          const Eigen::Vector3d n = Eigen::Vector3d(0.5, 1.0, 0.0).normalized();
          return (point - p0).dot(n);
      };

      auto sphere = [](const Eigen::Vector3d &point, const double &L, const double &c) {
          const double scale = 1.0;
          const double radius = 0.45;
          const double x = scale * (point.x() - c - 0.5), y = scale * (point.y() - c - 0.5), z = scale * (point.z() - c - 0.5);
          return x * x + y * y + z * z - radius * radius;
      };

      auto gyroid = [](const Eigen::Vector3d &point, const double &L, const double &c) {
          const double scale = (2.0 * std::numbers::pi) / L;
          const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();

          return std::sin(x) * std::cos(y) + std::sin(y) * std::cos(z) + std::sin(z) * std::cos(x) - c * scale;
      };

      auto gyroid_shifted = [&gyroid](const Eigen::Vector3d &point, const double &L, const double &c) {
          return gyroid(point, L, c + 0.1) * gyroid(point, L, c + 0.0);
      };

      auto diamond = [](const Eigen::Vector3d &point, const double &L, const double &c) {
          const double scale = (2.0 * std::numbers::pi) / L;
          const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
          return std::sin(x) * std::sin(y) * std::sin(z) + std::sin(x) * std::cos(y) * std::cos(z) +
                 std::cos(x) * std::sin(y) * std::cos(z) + std::cos(x) * std::cos(y) * std::sin(z) - c * scale;
      };

      auto diamond_shifted = [&diamond](const Eigen::Vector3d &point, const double &L, const double &c) {
          return diamond(point, L, c + 0.1) * diamond(point, L, c + 0.0);
      };

      auto split_p = [](const Eigen::Vector3d &point, const double &L, const double &c) {
          const double scale = (2.0 * std::numbers::pi) / L;
          const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
          return 1.1 * (std::sin(2.0 * x) * std::sin(z) * std::cos(y) + std::sin(2.0 * y) * std::sin(x) * std::cos(z) +
                        std::sin(2.0 * z) * std::sin(y) * std::cos(x)) -
                 0.2 * (std::cos(2.0 * x) * std::cos(2.0 * y) + std::cos(2.0 * y) * std::cos(2.0 * z) +
                        std::cos(2.0 * z) * std::cos(2.0 * x)) -
                 0.4 * (std::cos(2.0 * x) * std::cos(2.0 * y) * std::cos(2.0 * z)) - c * scale;
      };

      auto lidinoid = [](const Eigen::Vector3d &point, const double &L, const double &c) {
          const double scale = (2.0 * std::numbers::pi) / L;
          const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
          return std::sin(2.0 * x) * std::cos(y) * std::sin(z) + std::sin(2.0 * y) * std::cos(z) * std::sin(x) +
                 std::sin(2.0 * z) * std::cos(x) * std::sin(y) - std::cos(2.0 * x) * std::cos(2.0 * y) -
                 std::cos(2.0 * y) * std::cos(2.0 * z) - std::cos(2.0 * z) * std::cos(2.0 * x) - c * scale;
      };

      auto neovius = [](const Eigen::Vector3d &point, const double &L, const double &c) {
          const double scale = (2.0 * std::numbers::pi) / L;
          const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
          return 3.0 * (std::cos(x) + std::cos(y) + std::cos(z)) + 4.0 * std::cos(x) * std::cos(y) * std::cos(z) - c * scale;
      };

      auto surface_implicit_equation =
          [&plane_x, &generic_plane, &sphere, &gyroid, &diamond, &split_p, &lidinoid, &neovius, &gyroid_shifted, &diamond_shifted](
              const Eigen::MatrixXd &points) {
              float L = 1.0f;
              float c = 0.0f;

              Eigen::VectorXd surface_values(points.cols());
              for (unsigned int p = 0; p < surface_values.size(); p++)
                  surface_values[p] = sphere(points.col(p), L, c);

              return surface_values;
          };

      const auto marching_cubes_result = Gedim::MeshUtilities::MarchingCubes(geometry_utilities, surface_implicit_equation, mesh);

      const auto surface_intersection_graph =
          Gedim::MeshUtilities::ComputeMeshSurfaceIntersectionGraph(marching_cubes_result.surface_intersections,
                                                                    marching_cubes_result.surface_segments,
                                                                    marching_cubes_result.surface_triangles,
                                                                    mesh);

      Gedim::MeshMatrices conformed_mesh_data;
      Gedim::MeshMatricesDAO conformed_mesh(conformed_mesh_data);

      const auto conform_result = Gedim::MeshUtilities::ConformMeshToSurface(mesh,
                                                                             surface_intersection_graph,
                                                                             conformed_mesh);

      std::vector<unsigned int> cell3D_to_mantain;
      cell3D_to_mantain.reserve(conformed_mesh.Cell3DTotalNumber());
      for (unsigned int c = 0; c < conformed_mesh.Cell3DTotalNumber(); ++c)
      {
        switch (conform_result.new_cell3Ds_surface_position.at(c))
        {
          case Gedim::MeshUtilities::SurfacePositions::OnSurface:
          case Gedim::MeshUtilities::SurfacePositions::Positive:
            cell3D_to_mantain.push_back(c);
            break;
          case Gedim::MeshUtilities::SurfacePositions::Negative:
            break;
          default:
            throw std::runtime_error("Position not managed");
        }
      }
      cell3D_to_mantain.shrink_to_fit();

      const auto mesh_filter = Gedim::MeshUtilities::FilterMesh3D(cell3D_to_mantain,
                                                                  conformed_mesh);

      Gedim::MeshMatrices filtered_mesh_data;
      Gedim::MeshMatricesDAO filtered_mesh(filtered_mesh_data);

      Gedim::MeshUtilities::ExtractMesh3D(mesh_filter.Cell0Ds,
                                          mesh_filter.Cell1Ds,
                                          mesh_filter.Cell2Ds,
                                          mesh_filter.Cell3Ds,
                                          conformed_mesh,
                                          filtered_mesh);

      const auto cell2Ds_triangles = Gedim::MeshUtilities::ComputeMeshCell2DsTriangles(geometry_utilities,
                                                                                       Gedim::MeshUtilities::MeshCell2DsTrianglesData::TrianglesTypes::ByBarycenter,
                                                                                       filtered_mesh);

      const auto cell3Ds_tetra = Gedim::MeshUtilities::ComputeMeshCell3DsTetrahedra(geometry_utilities,
                                                                                    cell2Ds_triangles,
                                                                                    Gedim::MeshUtilities::MeshCell3DsTetrahedraData::TetrahedraTypes::ByBarycenter,
                                                                                    filtered_mesh);

      const auto extract_tetra = Gedim::MeshUtilities::ExtractTetrahedralMesh(cell2Ds_triangles,
                                                                              cell3Ds_tetra,
                                                                              filtered_mesh);

      Gedim::MeshUtilities::TetrahedralMeshDataToMeshDAO(extract_tetra,
                                                         filtered_mesh,
                                                         mesh);
      Gedim::MeshUtilities::ComputeCell0DCell1DNeighbours(mesh);
    }

    {
      const std::string mesh_folder = exportFolder + "/mesh";
      Gedim::Output::CreateFolder(mesh_folder);
      Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_folder, "mesh");
    }

    const double speed = 1.0;
    const unsigned int source = mesh.Cell0DTotalNumber() > 170 ? 170 : 0;
    const auto fmm_mesh_graph = Gedim::Algorithms::FastMarchingMethod::create_mesh_graph(mesh);
    const auto fmm_config = Gedim::Algorithms::FastMarchingMethod::create_configuration(geometry_utilities,
                                                                                        fmm_mesh_graph,
                                                                                        { source },
                                                                                        0.0,
                                                                                        speed,
                                                                                        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes::Linear);
    const auto result = Gedim::Algorithms::FastMarchingMethod::fast_marching_method(geometry_utilities,
                                                                                    fmm_mesh_graph,
                                                                                    fmm_config);
    {
      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
        {
          const auto& vertex_config = fmm_config.cell0Ds_initial_condition.at(c);

          double cell0D_init_type = 0.0;
          switch (vertex_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Source:
              cell0D_init_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Initialized:
              cell0D_init_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unreachable:
              cell0D_init_type = std::numeric_limits<double>::infinity();
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes::Unknown:
              cell0D_init_type = -1.0;
              break;
            default:
              throw std::runtime_error("status not managed");
          }

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
                {
                  {"init_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_init_type},
                  {"init_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.initial_time},
                  {"border_distance", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_config.border_distance},
                  {"status", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_status},
                  {"final_time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &result.vertices_time.at(c)}
                });
        }

        exporter.Export(exportFolder + "/cell0Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
        {
          const auto& edge_config = fmm_config.cell1Ds_restriction.at(e);

          double cell1D_restriction_type = 0.0;
          switch (edge_config.type)
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Free:
              cell1D_restriction_type = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Constrained:
              cell1D_restriction_type = 1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes::Unreachable:
              cell1D_restriction_type = std::numeric_limits<double>::infinity();
              break;
            default:
              throw std::runtime_error("restriction not managed");
          }

          exporter.AddSegment(
                mesh.Cell1DOriginCoordinates(e),
                mesh.Cell1DEndCoordinates(e),
                {
                  {"restriction_type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell1D_restriction_type},
                });
        }

        exporter.Export(exportFolder + "/cell1Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                             mesh.Cell2DsVertices(),
                             {{"time", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(result.vertices_time.size()), result.vertices_time.data()}});

        exporter.Export(exportFolder + "/cell2Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPolyhedrons(mesh.Cell0DsCoordinates(),
                             mesh.Cell3DsFacesVertices(),
                             {{"time", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(result.vertices_time.size()), result.vertices_time.data()}});

        exporter.Export(exportFolder + "/cell3Ds_time.vtu");
      }

    }
  }

} // namespace UnitTesting

#endif // __TEST_FastMarchingMethod_H
