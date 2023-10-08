#include "CommonUtilities.hpp"

using namespace std;

namespace Gedim
{
  // ***************************************************************************
  list<vector<int> > Utilities::CombinationWithRepetition(int n,
                                                          int k)
  {
    n--;
    vector<int> v(k+1, 0);
    list<vector<int> > result;

    while (true)
    {
      for (int i = 0; i < k; ++i){                //vai um
        if (v[i] > n){
          v[i + 1] += 1;
          for (int k = i; k >= 0; --k){
            v[k] = v[i + 1];
          }
          //v[i] = v[i + 1];
        }
      }

      if (v[k] > 0)
        break;

      result.push_back(vector<int>(v.begin(), v.end() - 1));
      v[0] += 1;
    }

    return result;
  }
  // ***************************************************************************
  std::vector<unsigned int> Utilities::RandomArrayNoRepetition(const unsigned int& n,
                                                               const unsigned int& maxNumber,
                                                               const unsigned int& seed)
  {
    if (n == 0)
      return { };

    Gedim::Output::Assert(n <= maxNumber);

    vector<unsigned int> randomNumbers(maxNumber);
    std::iota(begin(randomNumbers), end(randomNumbers), 0);
    Shuffle(randomNumbers, seed);
    randomNumbers.resize(n);

    return randomNumbers;
  }
  // ***************************************************************************
}
