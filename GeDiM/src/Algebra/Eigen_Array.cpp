#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"

#include "Eigen_Utilities.hpp"
#include "IOUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
template class Eigen_Array<VectorXd, SparseMatrix<double>>;
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>::SetValues(const vector<int> &indices, const vector<double> &values)
{
    Output::Assert(indices.size() == values.size());

    for (unsigned int i = 0; i < indices.size(); i++)
        SetValue(indices[i], values[i]);
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>::AddValues(const std::vector<int> &indices, const std::vector<double> &values)
{
    Output::Assert(indices.size() == values.size());

    for (unsigned int i = 0; i < indices.size(); i++)
        AddValue(indices[i], values[i]);
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
std::vector<double> Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>::GetValues(const std::vector<int> &indices) const
{
    std::vector<double> values;

    if (indices.size() == 0)
    {
        const unsigned int size = Size();
        values.resize(size);
        for (unsigned int i = 0; i < size; i++)
            values[i] = _vector[i];
    }
    else
    {
        values.resize(indices.size());
        for (unsigned int i = 0; i < indices.size(); i++)
            values[i] = _vector[indices.at(i)];
    }

    return values;
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>::SumMultiplication(const ISparseArray &A, const IArray &w)
{
    const Eigen_SparseArrayType &M =
        (const Eigen_SparseArrayType &)static_cast<const Eigen_SparseArray<Eigen_SparseArrayType> &>(A);

    switch (A.Type())
    {
    case ISparseArray::SparseArrayTypes::Symmetric:
        _vector += (Eigen::SparseMatrix<double>(M)).selfadjointView<Lower>() * Cast(w);
        break;
    case ISparseArray::SparseArrayTypes::None:
    case ISparseArray::SparseArrayTypes::Diagonal:
    case ISparseArray::SparseArrayTypes::Lower:
    case ISparseArray::SparseArrayTypes::Upper:
        _vector += M * Cast(w);
        break;
    default:
        throw std::runtime_error("Matrix type not managed");
    }
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>::SubtractionMultiplication(const ISparseArray &A, const IArray &w)
{

    const Eigen_SparseArrayType &M =
        (const Eigen_SparseArrayType &)static_cast<const Eigen_SparseArray<Eigen_SparseArrayType> &>(A);

    switch (A.Type())
    {
    case ISparseArray::SparseArrayTypes::Symmetric:
        _vector -= (Eigen::SparseMatrix<double>(M)).selfadjointView<Lower>() * Cast(w);
        break;
    case ISparseArray::SparseArrayTypes::None:
    case ISparseArray::SparseArrayTypes::Diagonal:
    case ISparseArray::SparseArrayTypes::Lower:
    case ISparseArray::SparseArrayTypes::Upper:
        _vector -= M * Cast(w);
        break;
    default:
        throw std::runtime_error("Matrix type not managed");
    }
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>::ToBinaryFile(const std::string &filePath,
                                                                       const unsigned int &dataSizeToWrite,
                                                                       const unsigned int &dataStartingPositionToWrite,
                                                                       const bool &append) const
{
    Eigen_Utilities::WriteToBinaryFile(filePath, _vector, dataSizeToWrite, dataStartingPositionToWrite, append);
}
// ***************************************************************************
} // namespace Gedim
