#ifndef __PETSc_Array_HPP
#define __PETSc_Array_HPP

#include "Gedim_Macro.hpp"

#if ENABLE_PETSC == 1

#include "petscmat.h"
#include "petscvec.h"

#include "IArray.hpp"

namespace Gedim
{
/// \brief PETSc column vector
template <typename PETSc_ArrayType = Vec, typename PETSc_SparseArrayType = Mat> class PETSc_Array final : public IArray
{
  private:
    PETSc_ArrayType _vector;

  public:
    PETSc_Array()
    {
    }
    PETSc_Array(const PETSc_ArrayType &vector)
    {
        _vector = vector;
    }
    PETSc_Array(PETSc_ArrayType &&vector)
    {
        _vector = std::move(vector);
    }
    ~PETSc_Array()
    {
    }

    operator PETSc_ArrayType &()
    {
        return _vector;
    }
    operator const PETSc_ArrayType &() const
    {
        return _vector;
    }
    inline PETSc_ArrayType &Cast(IArray &v)
    {
        return (PETSc_ArrayType &)static_cast<PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType> &>(v);
    }
    inline const PETSc_ArrayType &Cast(const IArray &v) const
    {
        return (const PETSc_ArrayType &)static_cast<const PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType> &>(v);
    }

    void Create();
    inline void Destroy()
    {
        VecDestroy(&_vector);
    }
    inline void SetSize(const unsigned int &numCols)
    {
        SetSizes(numCols, 0);
    }
    void SetSizes(const unsigned int &numCols, const unsigned int &numLocalCols = 0);

    unsigned int Size() const;
    const double *Data() const
    {
        throw std::runtime_error("unimplemented method");
    }
    double *Data()
    {
        throw std::runtime_error("unimplemented method");
    }

    void SetValue(const int &i, const double &val);
    void SetValues(const std::vector<int> &indices, const std::vector<double> &values);
    void AddValue(const int &i, const double &val);
    void AddValues(const std::vector<int> &indices, const std::vector<double> &values);
    std::vector<double> GetValues(const std::vector<int> &indices = {}) const;

    inline void Zeros()
    {
        VecSet(_vector, 0.0);
    }
    inline void Ones()
    {
        VecSet(_vector, 1.0);
    }
    inline void Constant(const double &c)
    {
        VecSet(_vector, c);
    }

    void SumMultiplication(const ISparseArray &A, const IArray &w);
    void SubtractionMultiplication(const ISparseArray &A, const IArray &w);

    double &operator[](const int &)
    {
        throw std::runtime_error("unimplemented method");
    }

    const PetscScalar &operator[](const int &) const
    {
        throw std::runtime_error("unimplemented method");
    }

    inline void Concatenate(const IArray &, const IArray &)
    {
        throw std::runtime_error("unimplemented method");
    }

    inline double GetValue(const int &i) const
    {
        PetscScalar value;
        VecGetValues(_vector, 1, &i, &value);
        return static_cast<double>(value);
    }

    std::ostream &Print(std::ostream &output) const;

    inline IArray &operator+=(const IArray &v)
    {
        VecAXPY(_vector, 1.0, Cast(v));
        return *this;
    }
    inline IArray &operator-=(const IArray &v)
    {
        VecAXPY(_vector, -1.0, Cast(v));
        return *this;
    }
    inline IArray &operator*=(const double &c)
    {
        VecScale(_vector, c);
        return *this;
    }
    inline IArray &operator/=(const double &c)
    {
        VecScale(_vector, 1.0 / c);
        return *this;
    }

    PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType> operator+(const PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType> &) const
    {
        throw std::runtime_error("unimplemented method");
    }

    PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType> operator-(const PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType> &) const
    {
        throw std::runtime_error("unimplemented method");
    }

    PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType> operator*(const double &) const
    {
        throw std::runtime_error("unimplemented method");
    }

    PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType> &operator=(PETSc_ArrayType &&)
    {
        throw std::runtime_error("unimplemented method");
    }

    void ToBinaryFile(const std::string &, const unsigned int & = 0, const unsigned int & = 0, const bool & = false) const
    {
        throw std::runtime_error("unimplemented method");
    }

    double Norm() const;

    double Dot(const IArray &v) const;

    inline void Copy(const IArray &v)
    {
        VecDuplicate(Cast(v), &_vector);
    }
};
} // namespace Gedim

#endif // ENABLE_PETSC

#endif // __PETSc_Array_HPP
