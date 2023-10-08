#ifndef __GEDIM_UTILITIES_H
#define __GEDIM_UTILITIES_H

#include "IOUtilities.hpp"
#include <algorithm>
#include <numeric>
#include <random>

namespace Gedim
{
  /// \brief Static class fot generic functions
  /// \copyright See top level LICENSE file for details.
  class Utilities
  {
    public:
      /// \brief Tells the compiler the parameter is unused
      template<class T>
      static void Unused(const T&) { }

      /// \brief create Combination With Repetition with n elements in k subset
      /// \param n size of elements
      /// \param k size of subset
      static std::list<std::vector<int>> CombinationWithRepetition(int n,
                                                                   int k);


      /// \brief Shuffle an array
      template<class container>
      static void Shuffle(container& array,
                          const unsigned int& seed = time(nullptr))
      {
        std::shuffle(array.begin(),
                     array.end(),
                     std::default_random_engine(seed));
      }

      /// \param n size of elements
      /// \param seed the side
      /// \return random array
      static std::vector<unsigned int> RandomArrayNoRepetition(const unsigned int& n,
                                                               const unsigned int& maxNumber = RAND_MAX,
                                                               const unsigned int& seed = time(nullptr));


      /// \brief Gives you back the array sorted indices
      /// \param array the array to sort
      /// \return the indices of the array ordered
      template<typename T>
      static std::vector<unsigned int> SortArrayIndices(const std::vector<T>& array)
      {
        std::vector<unsigned int> indices(array.size());
        std::iota(begin(indices), end(indices), 0);

        std::sort(indices.begin(),
                  indices.end(),
                  [&array](int a, int b) { return (array.at(a) < array.at(b)); });

        return indices;
      }

      /// \brief  Compute the regression line slope of points (x,y).
      /// \param  x_begin  Start of x range.
      /// \param  y_begin  Start of y range.
      /// \param  n the number of points
      /// \return  The slope.
      template<typename _InputIterator>
      static double Slope(_InputIterator x_begin,
                          _InputIterator y_begin,
                          const unsigned int& n)
      {
        const auto s_x  = std::accumulate(x_begin, x_begin + n, 0.0);
        const auto s_y  = std::accumulate(y_begin, y_begin + n, 0.0);
        const auto s_xx = std::inner_product(x_begin, x_begin + n, x_begin, 0.0);
        const auto s_xy = std::inner_product(x_begin, x_begin + n, y_begin, 0.0);
        return (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
      }
  };
}

#endif // __GEDIM_UTILITIES_H
