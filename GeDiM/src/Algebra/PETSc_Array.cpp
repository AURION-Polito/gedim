#include "PETSc_Array.hpp"

#if ENABLE_PETSC == 1
#include "PETSc_SparseArray.hpp"
#include "IOUtilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  template class PETSc_Array<Vec, Mat>;
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  void PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::Create()
  {
    VecAssemblyBegin(_vector);
    VecAssemblyEnd(_vector);
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  void PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::SetSizes(const unsigned int& numCols,
                                                                     const unsigned int& numLocalCols)
  {
    if (numLocalCols == 0)
    {
      VecCreate(PETSC_COMM_WORLD,
                &_vector);
      VecSetType(_vector,
                 VECMPI);
      VecSetSizes(_vector,
                  PETSC_DECIDE,numCols);
    }
    else
    {
      VecCreate(PETSC_COMM_WORLD,
                &_vector);
      VecSetType(_vector,
                 VECMPI);
      VecSetSizes(_vector,
                  numLocalCols,
                  PETSC_DETERMINE);
    }

    VecSetFromOptions(_vector);
    VecSetUp(_vector);

    Zeros();
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  unsigned int PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::Size() const
  {
    PetscInt N;
    VecGetSize(_vector, &N);
    return static_cast<unsigned int>(N);
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  void PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::SetValue(const int& i,
                                                                     const double& val)
  {
    VecSetValues(_vector,
                 1,
                 (PetscInt*)(&i),
                 (const PetscScalar*)(&val),
                 INSERT_VALUES);
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  void PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::SetValues(const std::vector<int>& indices,
                                                                      const std::vector<double>& values)
  {
    Output::Assert(indices.size() == values.size());

    VecSetValues(_vector,
                 indices.size(),
                 (PetscInt*)(indices.data()),
                 (const PetscScalar*)(values.data()),
                 INSERT_VALUES);
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  void PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::AddValue(const int& i,
                                                                     const double& val)
  {
    VecSetValues(_vector,
                 1,
                 (PetscInt*)(&i),
                 (PetscScalar*)&val,
                 ADD_VALUES);
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  void PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::AddValues(const std::vector<int>& indices,
                                                                      const std::vector<double>& values)
  {
    Output::Assert(indices.size() == values.size());

    VecSetValues(_vector,
                 indices.size(),
                 (PetscInt*)(indices.data()),
                 (const PetscScalar*)(values.data()),
                 ADD_VALUES);
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  std::vector<double> PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::GetValues(const std::vector<int>& indices) const
  {
    std::vector<double> values;

    const PetscScalar *array; // Pointer to the vector data
    VecGetArrayRead(_vector,
                    &array);

    if (indices.size() == 0)
    {
      const unsigned int size = Size();
      values.resize(size);
      for (unsigned int i = 0; i < size; i++)
        values[i] = array[i];
    }
    else
    {
      values.resize(indices.size());
      for (unsigned int i = 0; i < indices.size(); i++)
        values[i] = array[indices.at(i)];
    }

    VecRestoreArrayRead(_vector,
                        &array);

    return values;
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType,
           typename PETSc_SparseArrayType>
  void PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::SumMultiplication(const ISparseArray& A,
                                                                              const IArray& w)
  {
    const auto& w_PETSc_Array = static_cast<const PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>&>(w);
    const PETSc_ArrayType& w_PETSc = static_cast<const PETSc_ArrayType&>(w_PETSc_Array);

    const auto& A_PETSc_Array = static_cast<const PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>&>(A);
    const PETSc_SparseArrayType& A_PETSc = static_cast<const PETSc_SparseArrayType&>(A_PETSc_Array);

    unsigned int size = Size();
    Vec temp;
    VecCreate(PETSC_COMM_WORLD, &temp);
    VecSetSizes(temp, PETSC_DECIDE, size);
    VecSetFromOptions(temp);

    switch (A.Type())
    {
      case ISparseArray::SparseArrayTypes::Symmetric:
      case ISparseArray::SparseArrayTypes::None:
      case ISparseArray::SparseArrayTypes::Diagonal:
      case ISparseArray::SparseArrayTypes::Lower:
      case ISparseArray::SparseArrayTypes::Upper:
        MatMult(A_PETSc, w_PETSc, temp);
        VecAXPY(_vector, 1.0, temp);
        break;
      default:
        throw std::runtime_error("Matrix type not managed");
    }
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  void PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::SubtractionMultiplication(const ISparseArray& A,
                                                                                      const IArray& w)
  {
    const auto& w_PETSc_Array = static_cast<const PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>&>(w);
    const PETSc_ArrayType& w_PETSc = static_cast<const PETSc_ArrayType&>(w_PETSc_Array);

    const auto& A_PETSc_Array = static_cast<const PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>&>(A);
    const PETSc_SparseArrayType& A_PETSc = static_cast<const PETSc_SparseArrayType&>(A_PETSc_Array);

    unsigned int size = Size();
    Vec temp;
    VecCreate(PETSC_COMM_WORLD, &temp);
    VecSetSizes(temp, PETSC_DECIDE, size);
    VecSetFromOptions(temp);

    switch (A.Type())
    {
      case ISparseArray::SparseArrayTypes::Symmetric:
      case ISparseArray::SparseArrayTypes::None:
      case ISparseArray::SparseArrayTypes::Diagonal:
      case ISparseArray::SparseArrayTypes::Lower:
      case ISparseArray::SparseArrayTypes::Upper:
        MatMult(A_PETSc, w_PETSc, temp);
        VecAXPY(_vector, -1.0, temp);
        break;
      default:
        throw std::runtime_error("Matrix type not managed");
    }
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  std::ostream& PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::Print(std::ostream& output) const
  {
    const unsigned int size = Size();

    const PetscScalar *array; // Pointer to the vector data
    VecGetArrayRead(_vector,
                    &array);

    output << "[";
    for (unsigned int i = 0; i < size; ++i)
      output<< (i == 0 ? "" : ",") << array[i];
    output<< "]" << std::endl;

    // Restore the array
    VecRestoreArrayRead(_vector,
                        &array);

    return output;
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  double PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::Norm() const
  {
    double norm2;
    VecNorm(_vector, NORM_2, &norm2);
    return norm2;
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  double PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::Dot(const IArray& v) const
  {
    double dot;
    VecDot(_vector, Cast(v), &dot);
    return dot;
  }
  // ***************************************************************************
}

#endif
