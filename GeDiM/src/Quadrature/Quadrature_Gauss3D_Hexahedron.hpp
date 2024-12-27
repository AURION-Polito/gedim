#ifndef __Quadrature_Gauss3D_Hexahedron_H
#define __Quadrature_Gauss3D_Hexahedron_H

#include "QuadratureData.hpp"

namespace Gedim
{
namespace Quadrature
{
/// Gauss quadrature rule for Hexahedrons
class Quadrature_Gauss3D_Hexahedron
{
  public:
    static QuadratureData FillPointsAndWeights(const unsigned int &order);
};
} // namespace Quadrature
} // namespace Gedim

#endif // __Quadrature_Gauss3D_Hexahedron_H
