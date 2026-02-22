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

#ifndef __TEST_meshUtilities2D_ImportTriangularMesh_H
#define __TEST_meshUtilities2D_ImportTriangularMesh_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshUtilities.hpp"

namespace GedimUnitTesting
{
TEST(TestMeshUtilities, TestMeshUtilities_2D_ImportTriangularMesh)
{
    std::string exportFolder = "./Export/TestMeshUtilities/TestMeshUtilities_2D_ImportTriangularMesh";
    Gedim::Output::CreateFolder(exportFolder);
}
} // namespace GedimUnitTesting

#endif // __TEST_MESH_H
