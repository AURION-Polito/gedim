#include "PETSc_Array.hpp"

#if ENABLE_PETSC

#include "IOUtilities.hpp"

using namespace std;

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
  void PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::SetValues(const vector<int>& indices,
                                                                      const vector<double>& values)
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
  void PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::SumMultiplication(const ISparseArray& A,
                                                                              const IArray& w)
  {
    throw std::runtime_error("unimplemented method");
    //    const PETSc_SparseArrayType& M = (const PETSc_SparseArrayType&)static_cast<const Mat&>(A);

    //    switch (A.Type())
    //    {
    //      case ISparseArray::SparseArrayTypes::Symmetric:
    //        _vector += (PETSc::SparseMatrix<double>(M)).selfadjointView<Lower>() *
    //                   Cast(w);
    //        break;
    //      case ISparseArray::SparseArrayTypes::None:
    //      case ISparseArray::SparseArrayTypes::Diagonal:
    //      case ISparseArray::SparseArrayTypes::Lower:
    //      case ISparseArray::SparseArrayTypes::Upper:
    //        _vector += M * Cast(w);
    //        break;
    //      default:
    //        throw std::runtime_error("Matrix type not managed");
    //    }
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  void PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::SubtractionMultiplication(const ISparseArray& A,
                                                                                      const IArray& w)
  {
    throw std::runtime_error("unimplemented method");
    //    const PETSc_SparseArrayType& M = (const PETSc_SparseArrayType&)static_cast<const PETSc_SparseArray<PETSc_SparseArrayType>&>(A);

    //    switch (A.Type())
    //    {
    //      case ISparseArray::SparseArrayTypes::Symmetric:
    //        _vector -= (PETSc::SparseMatrix<double>(M)).selfadjointView<Lower>() *
    //                   Cast(w);
    //        break;
    //      case ISparseArray::SparseArrayTypes::None:
    //      case ISparseArray::SparseArrayTypes::Diagonal:
    //      case ISparseArray::SparseArrayTypes::Lower:
    //      case ISparseArray::SparseArrayTypes::Upper:
    //        _vector -= M * Cast(w);
    //        break;
    //      default:
    //        throw std::runtime_error("Matrix type not managed");
    //    }
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType, typename PETSc_SparseArrayType>
  ostream& PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>::Print(std::ostream& output) const
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
