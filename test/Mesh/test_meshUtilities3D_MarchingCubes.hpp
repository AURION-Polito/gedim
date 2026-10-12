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
    std::string exportFolder = "./Export/TestMeshUtilities/TestMeshUtilities_3D_MarchingCubes";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshUtilities mesh_utilities;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
}

} // namespace GedimUnitTesting

#endif // __TEST_meshUtilities3D_MarchingCubes_H
