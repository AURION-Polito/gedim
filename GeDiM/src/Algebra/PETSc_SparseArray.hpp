#ifndef __PETSc_SparseArray_HPP
#define __PETSc_SparseArray_HPP

#include "Macro.hpp"

#if ENABLE_PETSC == 1

#include "ISparseArray.hpp"

#include "IOUtilities.hpp"

#include "petscmat.h"

namespace Gedim
{
  /// \brief PETSc sparse array
  template<typename PETSc_ArrayType = Vec,
           typename PETSc_SparseArrayType = Mat>
  class PETSc_SparseArray final : public ISparseArray
  {
    private:
      PETSc_SparseArrayType _matrix; ///< internal matrix
      SparseArrayTypes _matrixType; ///< matrix type

    public:
      PETSc_SparseArray()
      {
        _matrixType = SparseArrayTypes::None;
      }
      ~PETSc_SparseArray()
      {
      }

      PETSc_SparseArray(const PETSc_SparseArrayType& matrix,
                        const SparseArrayTypes& type = SparseArrayTypes::None)
      {
        _matrix = matrix;
        _matrixType = type;
      }

      operator PETSc_SparseArrayType&() { return _matrix; }
      operator const PETSc_SparseArrayType&() const { return _matrix; }
      inline PETSc_SparseArrayType& Cast(ISparseArray& v)
      { return (PETSc_SparseArrayType&)static_cast<PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>&>(v); }
      inline const PETSc_SparseArrayType& Cast(const ISparseArray& v)
      { return (const PETSc_SparseArrayType&)static_cast<const PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>&>(v); }
      inline const PETSc_SparseArrayType& Cast(const ISparseArray& v) const
      { return (const PETSc_SparseArrayType&)static_cast<const PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>&>(v); }

      void SetSize(const unsigned int& numRows,
                   const unsigned int& numCols,
                   const SparseArrayTypes& type = SparseArrayTypes::None);

      inline SparseArrayTypes Type() const
      { return _matrixType; }

      void Create();
      void Destroy();
      void Flush();
      void Reset();

      void Triplet(const unsigned int& i,
                   const unsigned int& j,
                   const double& value);

      void Triplets(const std::vector<unsigned int>& i,
                    const std::vector<unsigned int>& j,
                    const std::vector<double>& values);

      std::ostream& Print(std::ostream& output) const;

      inline ISparseArray& operator+=(const ISparseArray& A)
      {
        MatAXPY(_matrix, 1.0, Cast(A), DIFFERENT_NONZERO_PATTERN);
        return *this;
      }
      inline ISparseArray& operator-=(const ISparseArray& A)
      {
        MatAXPY(_matrix, -1.0, Cast(A), DIFFERENT_NONZERO_PATTERN);
        return *this;
      }
      inline ISparseArray& operator*=(const double& c)
      {
        MatScale(_matrix, c);
        return *this;
      }
      inline ISparseArray& operator/=(const double& c)
      {
        MatScale(_matrix, 1.0 / c);
        return *this;
      }

      unsigned int rows() const;
      unsigned int cols() const;

      PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>& operator=(PETSc_SparseArrayType&&)
      { throw std::runtime_error("unimplemented method"); }

      void ToBinaryFile(const std::string&,
                        const bool&) const
      { throw std::runtime_error("unimplemented method"); }

      inline double Norm() const
      {
        double norm;
        MatNorm(_matrix, NORM_FROBENIUS, &norm);
        return norm;
      }

      inline double Cond(const ISparseArray::ConditionNumberAlgorithm& algorithm) const
      { throw std::runtime_error("unimplemented method"); }

      unsigned int NonZeros() const;
  };
}

#endif // ENABLE_PETSC

#endif // __PETSc_SparseArray_HPP
