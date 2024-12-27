#include "PETSc_SparseArray.hpp"

#if ENABLE_PETSC == 1

namespace Gedim
{
// ***************************************************************************
template class PETSc_SparseArray<Vec, Mat>;
// ***************************************************************************
template <typename PETSc_ArrayType, typename PETSc_SparseArrayType>
void PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>::SetSize(const unsigned int &numRows,
                                                                        const unsigned int &numCols,
                                                                        const SparseArrayTypes &type)
{
    MatCreate(PETSC_COMM_WORLD, &_matrix);
    MatSetType(_matrix, MATAIJ);
    MatSetSizes(_matrix, PETSC_DECIDE, PETSC_DECIDE, numRows, numCols);
    MatSetFromOptions(_matrix);
    MatSetUp(_matrix);

    _matrixType = type;
}
// ***************************************************************************
template <typename PETSc_ArrayType, typename PETSc_SparseArrayType>
void PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>::Create()
{
    MatAssemblyBegin(_matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_matrix, MAT_FINAL_ASSEMBLY);
}
// ***************************************************************************
template <typename PETSc_ArrayType, typename PETSc_SparseArrayType>
void PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>::Destroy()
{
    MatDestroy(&_matrix);
}
// ***************************************************************************
template <typename PETSc_ArrayType, typename PETSc_SparseArrayType>
void PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>::Flush()
{
    MatAssemblyBegin(_matrix, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(_matrix, MAT_FLUSH_ASSEMBLY);
}
// ***************************************************************************
template <typename PETSc_ArrayType, typename PETSc_SparseArrayType>
void PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>::Reset()
{
    MatZeroEntries(_matrix);
}
// ***************************************************************************
template <typename PETSc_ArrayType, typename PETSc_SparseArrayType>
void PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>::Triplet(const unsigned int &i, const unsigned int &j, const double &value)
{
    switch (_matrixType)
    {
    case SparseArrayTypes::Symmetric:
    case SparseArrayTypes::Lower:
    case SparseArrayTypes::Upper:
    case SparseArrayTypes::Diagonal:
    case SparseArrayTypes::None: {
        MatSetValue(_matrix, i, j, value, ADD_VALUES);
    }
    break;
    default:
        throw std::runtime_error("Matrix type not supported");
    }
}
// ***************************************************************************
template <typename PETSc_ArrayType, typename PETSc_SparseArrayType>
void PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>::Triplets(const std::vector<unsigned int> &i,
                                                                         const std::vector<unsigned int> &j,
                                                                         const std::vector<double> &values)
{
    if (i.size() != j.size() || i.size() != values.size())
        throw std::runtime_error("Invalid triplets size");

    switch (_matrixType)
    {
    case SparseArrayTypes::Symmetric:
    case SparseArrayTypes::Lower:
    case SparseArrayTypes::Upper:
    case SparseArrayTypes::Diagonal:
    case SparseArrayTypes::None: {
        for (unsigned int s = 0; s < i.size(); ++s)
        {
            MatSetValue(_matrix, i[s], j[s], values[s], ADD_VALUES);
        }
    }
    break;
    default:
        throw std::runtime_error("Matrix type not supported");
    }
}
// ***************************************************************************
template <typename PETSc_ArrayType, typename PETSc_SparseArrayType>
std::ostream &PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>::Print(std::ostream &output) const
{
    PetscViewer viewer;
    PetscViewerCreate(PETSC_COMM_SELF, &viewer);

    // Create a temporary file for storing the matrix output
    PetscViewerASCIIOpen(PETSC_COMM_SELF, "temp_matrix_output.txt", &viewer);
    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE); // Optionally set format

    MatView(_matrix, viewer);

    PetscViewerDestroy(&viewer);

    // Read the file into a C++ ostringstream
    std::ifstream file("temp_matrix_output.txt");
    std::ostringstream oss;
    oss << file.rdbuf();
    file.close();

    // Write to the provided ostream
    output << oss.str();

    // Clean up the temporary file
    std::remove("temp_matrix_output.txt");

    return output;
}
// ***************************************************************************
template <typename PETSc_ArrayType, typename PETSc_SparseArrayType>
unsigned int PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>::rows() const
{
    PetscInt rows, cols;
    MatGetSize(_matrix, &rows, &cols);
    return static_cast<unsigned int>(rows);
}
// ***************************************************************************
template <typename PETSc_ArrayType, typename PETSc_SparseArrayType>
unsigned int PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>::cols() const
{
    PetscInt rows, cols;
    MatGetSize(_matrix, &rows, &cols);
    return static_cast<unsigned int>(cols);
}
// ***************************************************************************
template <typename PETSc_ArrayType, typename PETSc_SparseArrayType>
unsigned int PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>::NonZeros() const
{
    PetscObjectState non_zeros;
    MatGetNonzeroState(_matrix, &non_zeros);
    return static_cast<unsigned int>(non_zeros);
}
// ***************************************************************************
} // namespace Gedim

#endif
