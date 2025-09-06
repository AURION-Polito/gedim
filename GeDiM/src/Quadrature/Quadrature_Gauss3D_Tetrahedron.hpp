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

#ifndef __Quadrature_Gauss3D_Tetrahedron_H
#define __Quadrature_Gauss3D_Tetrahedron_H

#include "QuadratureData.hpp"

namespace Gedim
{
namespace Quadrature
{
/// Gauss quadrature rule for Tetrahedrons
class Quadrature_Gauss3D_Tetrahedron
{
  public:
    static QuadratureData FillPointsAndWeights(const unsigned int &order);
};
} // namespace Quadrature
} // namespace Gedim

#endif // __Quadrature_Gauss3D_Tetrahedron_H
