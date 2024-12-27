#ifndef __Quadrature_Gauss2D_Triangle_H
#define __Quadrature_Gauss2D_Triangle_H

#include "QuadratureData.hpp"

namespace Gedim
{
namespace Quadrature
{
/// Gauss quadrature rule for triangles
class Quadrature_Gauss2D_Triangle
{
  public:
    static QuadratureData FillPointsAndWeights(const unsigned int &order);
};
} // namespace Quadrature
} // namespace Gedim

#endif // __Quadrature_Gauss2D_Triangle_H
