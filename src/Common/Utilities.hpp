#ifndef __UTILITIES_H
#define __UTILITIES_H

#include "Output.hpp"

using namespace std;

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
      static list<vector<int>> CombinationWithRepetition(int n,
                                                         int k);
  };
}

#endif // __UTILITIES_H
