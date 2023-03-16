#include "CommonUtilities.hpp"
#include <numeric>
#include <algorithm>

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
  void Utilities::Shuffle(std::vector<unsigned int>& array,
                          const unsigned int& seed)
  {
    srand(seed);
    const unsigned int n = array.size();

    for (unsigned int i = 0; i < n - 1; i++)
    {
      const unsigned int j = (i + rand() / (RAND_MAX / (n - i) + 1));
      const unsigned int temp = array[j];
      array[j] = array[i];
      array[i] = temp;
    }
  }
  // ***************************************************************************
  std::vector<unsigned int> Utilities::RandomArrayNoRepetition(const unsigned int& n,
                                                               const unsigned int& maxNumber,
                                                               const unsigned int& seed)
  {
    Output::Assert(n > 0);
    vector<unsigned int> randomNumbers(maxNumber);
    std::iota(begin(randomNumbers), end(randomNumbers), 0);
    Shuffle(randomNumbers, seed);
    randomNumbers.resize(n);

    return randomNumbers;
  }
  // ***************************************************************************
  template<typename T>
  std::vector<unsigned int> Utilities::SortArrayIndices(const std::vector<T>& array)
  {
    vector<unsigned int> indices(array.size());
    std::iota(begin(indices), end(indices), 0);

    std::sort(indices.begin(), indices.end(), [&array](int a, int b) { return (array.at(a) < array.at(b)); });

    return indices;
  }
  // ***************************************************************************
}
