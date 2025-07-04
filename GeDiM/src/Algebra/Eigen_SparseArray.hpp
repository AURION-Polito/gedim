// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __Eigen_SparseArray_HPP
#define __Eigen_SparseArray_HPP

#include "Eigen/Eigen"
#include "IOUtilities.hpp"
#include "ISparseArray.hpp"
#include "LAPACK_utilities.hpp"

#if ENABLE_SUITESPARSE == 1
#include "SuiteSparse_Utilities.hpp"
#endif

namespace Gedim
{
/// \brief Eigen sparse array
template <typename Eigen_ArrayType = Eigen::VectorXd, typename Eigen_SparseArrayType = Eigen::SparseMatrix<double>>
class Eigen_SparseArray final : public ISparseArray
{
  private:
    Eigen_SparseArrayType _matrix; ///< internal matrix
    SparseArrayTypes _matrixType;  ///< matrix type
    std::list<Eigen::Triplet<double>> _triplets;

  public:
    Eigen_SparseArray()
    {
        _matrixType = SparseArrayTypes::None;
    }
    ~Eigen_SparseArray()
    {
        Reset();
    }

    Eigen_SparseArray(const Eigen_SparseArrayType &matrix, const SparseArrayTypes &type = SparseArrayTypes::None)
    {
        _matrix = matrix;
        _matrixType = type;
    }

    operator Eigen_SparseArrayType &()
    {
        return _matrix;
    }
    operator const Eigen_SparseArrayType &() const
    {
        return _matrix;
    }
    inline Eigen_SparseArrayType &Cast(ISparseArray &v)
    {
        return (Eigen_SparseArrayType &)static_cast<Eigen_SparseArray<Eigen_ArrayType, Eigen_SparseArrayType> &>(v);
    }
    inline const Eigen_SparseArrayType &Cast(const ISparseArray &v)
    {
        return (const Eigen_SparseArrayType &)static_cast<const Eigen_SparseArray<Eigen_ArrayType, Eigen_SparseArrayType> &>(v);
    }
    inline const Eigen_SparseArrayType &Cast(const ISparseArray &v) const
    {
        return (const Eigen_SparseArrayType &)static_cast<const Eigen_SparseArray<Eigen_ArrayType, Eigen_SparseArrayType> &>(v);
    }

    inline void SetSize(const unsigned int &numRows, const unsigned int &numCols, const SparseArrayTypes &type = SparseArrayTypes::None)
    {
        _matrix.resize(numRows, numCols);
        _matrixType = type;
    }

    inline SparseArrayTypes Type() const
    {
        return _matrixType;
    }

    void Create();
    inline void Destroy()
    {
        Reset();
    }
    inline void Flush()
    {
    }
    inline void Reset()
    {
        _triplets.clear();
        _matrix.prune(0.0);
    }

    void Triplet(const unsigned int &i, const unsigned int &j, const double &value);

    void Triplets(const std::vector<unsigned int> &i, const std::vector<unsigned int> &j, const std::vector<double> &values);

    inline std::ostream &Print(std::ostream &output) const
    {
        return output << _matrix;
    }

    inline ISparseArray &operator+=(const ISparseArray &A)
    {
        Gedim::Output::Assert(A.Type() == Type());
        _matrix += Cast(A);
        return *this;
    }
    inline ISparseArray &operator-=(const ISparseArray &A)
    {
        Gedim::Output::Assert(A.Type() == Type());
        _matrix -= Cast(A);
        return *this;
    }
    inline ISparseArray &operator*=(const double &c)
    {
        _matrix *= c;
        return *this;
    }
    inline ISparseArray &operator/=(const double &c)
    {
        _matrix /= c;
        return *this;
    }

    inline void Copy(const ISparseArray &A)
    {
        _matrix = Cast(A);
    }

    inline unsigned int rows() const
    {
        return _matrix.rows();
    }

    inline unsigned int cols() const
    {
        return _matrix.cols();
    }

    Eigen_SparseArray<Eigen_ArrayType, Eigen_SparseArrayType> &operator=(Eigen_SparseArrayType &&matrix)
    {
        _matrix = std::move(matrix);
        return *this;
    }

    void ToBinaryFile(const std::string &filePath, const bool &append = false) const;

    inline double Norm() const
    {
        return _matrix.norm();
    }

    inline double Cond(const ISparseArray::ConditionNumberAlgorithm &algorithm = ISparseArray::ConditionNumberAlgorithm::SVDLapack) const
    {
        switch (algorithm)
        {
        case ISparseArray::ConditionNumberAlgorithm::SVDLapack: // Lapack -> cond 2
            return LAPACK_utilities::cond(LAPACK_utilities::svd(Eigen::MatrixXd(_matrix)));
        case ISparseArray::ConditionNumberAlgorithm::CondestLapack: // Lapack -> cond 1
            return 1.0 / LAPACK_utilities::rcondest(_matrix);
#if ENABLE_SUITESPARSE == 1
        case ISparseArray::ConditionNumberAlgorithm::CondestSuiteSparse: // SuiteSparse
            return Gedim::SuiteSparse_Utilities::condest(_matrix);
#endif
        default:
            throw std::runtime_error("Not valid method to compute the condition number.");
        }
    }

    inline unsigned int NonZeros() const
    {
        return _matrix.nonZeros();
    }
};
} // namespace Gedim

#endif // __Eigen_SparseArray_HPP
