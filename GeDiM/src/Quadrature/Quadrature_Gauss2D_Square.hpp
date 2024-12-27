#ifndef __Quadrature_Gauss2D_Square_H
#define __Quadrature_Gauss2D_Square_H

#include "QuadratureData.hpp"

namespace Gedim
{
namespace Quadrature
{
/// Gauss quadrature rule for Squares
class Quadrature_Gauss2D_Square
{
  public:
    static QuadratureData FillPointsAndWeights(const unsigned int &order);
};
} // namespace Quadrature
} // namespace Gedim

#endif // __Quadrature_Gauss2D_Square_H
