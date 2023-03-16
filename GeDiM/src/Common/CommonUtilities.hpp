#ifndef __GEDIM_UTILITIES_H
#define __GEDIM_UTILITIES_H

#include "IOUtilities.hpp"

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
      static void Shuffle(std::vector<unsigned int>& array,
                          const unsigned int& seed = time(nullptr));

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
      static std::vector<unsigned int> SortArrayIndices(const std::vector<T>& array);
  };
}

#endif // __GEDIM_UTILITIES_H
