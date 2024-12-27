#include "Eigen_SparseArray.hpp"
#include "Eigen_Utilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
// ***************************************************************************
template class Eigen_SparseArray<VectorXd, SparseMatrix<double>>;
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_SparseArray<Eigen_ArrayType, Eigen_SparseArrayType>::Create()
{
    _matrix.setFromTriplets(_triplets.begin(), _triplets.end());
    _matrix.makeCompressed();
    _triplets.clear();
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_SparseArray<Eigen_ArrayType, Eigen_SparseArrayType>::Triplet(const unsigned int &i, const unsigned int &j, const double &value)
{
    switch (_matrixType)
    {
    case SparseArrayTypes::None:
        _triplets.push_back(Eigen::Triplet<double>(i, j, value));
        break;
    case SparseArrayTypes::Symmetric: // store only lower part
        if (j <= i)
            _triplets.push_back(Eigen::Triplet<double>(i, j, value));
        break;
    case SparseArrayTypes::Lower:
        if (j <= i)
            _triplets.push_back(Eigen::Triplet<double>(i, j, value));
        break;
    case SparseArrayTypes::Upper:
        if (j >= i)
            _triplets.push_back(Eigen::Triplet<double>(i, j, value));
        break;
    case SparseArrayTypes::Diagonal:
        if (j == i)
            _triplets.push_back(Eigen::Triplet<double>(i, j, value));
        break;
    default:
        throw runtime_error("Matrix type not supported");
    }
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_SparseArray<Eigen_ArrayType, Eigen_SparseArrayType>::Triplets(const vector<unsigned int> &i,
                                                                         const vector<unsigned int> &j,
                                                                         const vector<double> &values)
{
    if (i.size() != j.size() || i.size() != values.size())
        throw runtime_error("Invalid triplets size");

    unsigned int numTriplets = i.size();
    for (unsigned int t = 0; t < numTriplets; t++)
        Triplet(i[t], j[t], values[t]);
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_SparseArray<Eigen_ArrayType, Eigen_SparseArrayType>::ToBinaryFile(const std::string &filePath, const bool &append) const
{
    switch (_matrixType)
    {
    case SparseArrayTypes::None:
    case SparseArrayTypes::Diagonal:
        Eigen_Utilities::WriteToBinaryFile(filePath, _matrix, (Eigen::UpLoType)0, append);
        break;
    case SparseArrayTypes::Lower:
        Eigen_Utilities::WriteToBinaryFile(filePath, _matrix, Eigen::UpLoType::Lower, append);
        break;
    case SparseArrayTypes::Upper:
        Eigen_Utilities::WriteToBinaryFile(filePath, _matrix, Eigen::UpLoType::Upper, append);
        break;
    case SparseArrayTypes::Symmetric:
        Eigen_Utilities::WriteToBinaryFile(
            filePath,
            Eigen::SparseMatrix<double>(Eigen::SparseMatrix<double>(_matrix).selfadjointView<Lower>()),
            (Eigen::UpLoType)0,
            append);
        break;
    default:
        throw runtime_error("Matrix type not managed");
    }
}
// ***************************************************************************
} // namespace Gedim
