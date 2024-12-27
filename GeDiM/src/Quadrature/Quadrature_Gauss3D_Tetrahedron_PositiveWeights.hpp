#ifndef __Quadrature_Gauss3D_Tetrahedron_PositiveWeights_H
#define __Quadrature_Gauss3D_Tetrahedron_PositiveWeights_H

#include "QuadratureData.hpp"

namespace Gedim
{
namespace Quadrature
{
/// Gauss quadrature rule for Tetrahedrons
class Quadrature_Gauss3D_Tetrahedron_PositiveWeights
{
  public:
    static QuadratureData FillPointsAndWeights(const unsigned int &order);
};
} // namespace Quadrature
} // namespace Gedim

#endif // __Quadrature_Gauss3D_Tetrahedron_PositiveWeights_H
