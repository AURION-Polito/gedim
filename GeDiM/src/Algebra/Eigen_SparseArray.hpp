#ifndef __Eigen_SparseArray_HPP
#define __Eigen_SparseArray_HPP

#include "ISparseArray.hpp"
#include "Eigen/Eigen"
#include "LAPACK_utilities.hpp"
#include "IOUtilities.hpp"

#if ENABLE_SUITESPARSE == 1
#include "SuiteSparse_Utilities.hpp"
#endif

namespace Gedim
{
/// \brief Eigen sparse array
template<typename Eigen_ArrayType = Eigen::VectorXd,
         typename Eigen_SparseArrayType = Eigen::SparseMatrix<double>>
class Eigen_SparseArray final : public ISparseArray
{
private:
    Eigen_SparseArrayType _matrix; ///< internal matrix
    SparseArrayTypes _matrixType; ///< matrix type
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

    operator Eigen_SparseArrayType&() { return _matrix; }
    operator const Eigen_SparseArrayType&() const { return _matrix; }
    inline Eigen_SparseArrayType& Cast(ISparseArray& v)
    { return (Eigen_SparseArrayType&)static_cast<Eigen_SparseArray<Eigen_ArrayType, Eigen_SparseArrayType>&>(v); }
    inline const Eigen_SparseArrayType& Cast(const ISparseArray& v)
    { return (const Eigen_SparseArrayType&)static_cast<const Eigen_SparseArray<Eigen_ArrayType, Eigen_SparseArrayType>&>(v); }

    inline void SetSize(const unsigned int& numRows,
                        const unsigned int& numCols,
                        const SparseArrayTypes& type = SparseArrayTypes::None)
    {
        _matrix.resize(numRows, numCols);
        _matrixType = type;
    }

    inline SparseArrayTypes Type() const
    { return _matrixType; }

    void Create();
    inline void Flush() { }
    inline void Reset()
    {
        _triplets.clear();
        _matrix.prune(0.0);
    }

    void Triplet(const unsigned int& i,
                 const unsigned int& j,
                 const double& value);

    void Triplets(const std::vector<unsigned int>& i,
                  const std::vector<unsigned int>& j,
                  const std::vector<double>& values);

    inline std::ostream& Print(std::ostream& output) const
    { return output<< _matrix; }

    inline ISparseArray& operator+=(const ISparseArray& A)
    {
        Gedim::Output::Assert(A.Type() == Type());
        _matrix += Cast(A);
        return *this;
    }
    inline ISparseArray& operator-=(const ISparseArray& A)
    {
        Gedim::Output::Assert(A.Type() == Type());
        _matrix -= Cast(A);
        return *this;
    }
    inline ISparseArray& operator*=(const double& c)
    {
        _matrix *= c;
        return *this;
    }
    inline ISparseArray& operator/=(const double& c)
    {
        _matrix /= c;
        return *this;
    }

    Eigen_SparseArray<Eigen_ArrayType, Eigen_SparseArrayType>& operator=(Eigen_SparseArrayType&& matrix)
    {
        _matrix = std::move(matrix);
        return *this;
    }

    void ToBinaryFile(const std::string& filePath,
                      const bool& append = false) const;

    inline double Norm() const
    { return _matrix.norm(); }

    inline double Cond(unsigned int condType = 0) const
    {
        switch(condType)
        {
        case 0: // Lapack -> cond 2
            return LAPACK_utilities::cond(LAPACK_utilities::svd(Eigen::MatrixXd(_matrix)));
        case 1: // Lapack -> cond 1
            return 1.0 / LAPACK_utilities::rcondest(_matrix);
#if ENABLE_SUITESPARSE == 1
        case 2: // SuiteSparse
            return Gedim::SuiteSparse_Utilities::condest(_matrix);
#endif
        default:
            throw std::runtime_error("Not valid condition type.");
        }
    }

    inline unsigned int NonZeros() const
    { return _matrix.nonZeros(); }
};
}

#endif // __Eigen_SparseArray_HPP
