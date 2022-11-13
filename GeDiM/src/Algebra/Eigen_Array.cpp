#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"

#include "IOUtilities.hpp"
#include "Eigen_Utilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  template class Eigen_Array<VectorXd, SparseMatrix<double>>;
  // ***************************************************************************
  template<typename Eigen_ArrayType, typename Eigen_SparseArrayType>
  void Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>::SetValues(const vector<int>& indices,
                                                                      const vector<double>& values)
  {
    Output::Assert(indices.size() == values.size());

    for (unsigned int i = 0; i < indices.size(); i++)
      SetValue(indices[i], values[i]);
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType, typename Eigen_SparseArrayType>
  void Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>::AddValues(const std::vector<int>& indices, const std::vector<double>& values)
  {
    Output::Assert(indices.size() == values.size());

    for (unsigned int i = 0; i < indices.size(); i++)
      AddValue(indices[i], values[i]);
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType, typename Eigen_SparseArrayType>
  void Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>::SumMultiplication(const ISparseArray& A, const IArray& w)
  {
    _vector += (const Eigen_SparseArrayType&)static_cast<const Eigen_SparseArray<Eigen_SparseArrayType>&>(A) *
               Cast(w);
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType, typename Eigen_SparseArrayType>
  void Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>::SubtractionMultiplication(const ISparseArray& A, const IArray& w)
  {
    _vector -= (const Eigen_SparseArrayType&)static_cast<const Eigen_SparseArray<Eigen_SparseArrayType>&>(A) *
               Cast(w);
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType, typename Eigen_SparseArrayType>
  void Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>::ToBinaryFile(const std::string& filePath,
                                                                         const unsigned int& dataSizeToWrite,
                                                                         const unsigned int& dataStartingPositionToWrite,
                                                                         const bool& append) const
  {
    Eigen_Utilities::WriteToBinaryFile(filePath,
                                       _vector,
                                       dataSizeToWrite,
                                       dataStartingPositionToWrite,
                                       append);
  }
  // ***************************************************************************
}
