#ifndef __EIGEN_UTILITIES_H
#define __EIGEN_UTILITIES_H

#include <string>
#include <iostream>

#include "Eigen/Eigen"
#include "IOUtilities.hpp"

namespace Gedim
{
std::ostream& operator<<(std::ostream& out, const Eigen::VectorXd& vec);
std::ostream& operator<<(std::ostream& out, const Eigen::Map<const Eigen::VectorXd>& vec);
std::ostream& operator<<(std::ostream& out, const Eigen::MatrixXd& mat);

class Eigen_Utilities
{
private:

public:
    template <typename Derived>
    static Eigen::Map<const Eigen::SparseMatrix<Derived>> SparseMatrixToMap(const unsigned int& rows,
                                                                            const unsigned int& cols,
                                                                            const unsigned int& nnz,
                                                                            const int* innerVectors,
                                                                            const int* innerIndeces,
                                                                            const Derived* values)
    {
        return Eigen::Map<const Eigen::SparseMatrix<Derived>>(rows,
                                                              cols,
                                                              nnz,
                                                              innerVectors,
                                                              innerIndeces,
                                                              values);
    }

    template <typename Derived>
    static Eigen::Map<const Eigen::SparseMatrix<Derived>> SparseMatrixToMap(const Eigen::Ref<const Eigen::SparseMatrix<Derived>> sparseMatrix)
    {
        return SparseMatrixToMap(sparseMatrix.rows(),
                                 sparseMatrix.cols(),
                                 sparseMatrix.nonZeros(),
                                 sparseMatrix.outerIndexPtr(),
                                 sparseMatrix.innerIndexPtr(),
                                 sparseMatrix.valuePtr());
    }

    template <typename Derived, int Rows, int Cols>
    static Eigen::Map<const Eigen::Array<Derived, Rows, Cols>> ArrayToMap(const int& rows,
                                                                          const int& cols,
                                                                          const Derived* values)
    {
        return Eigen::Map<const Eigen::Array<Derived, Rows, Cols>>(values, rows, cols);
    }

    /// \brief Convert a matrix expressions to map
    /// \param matrix A const Eigen::Matrix
    /// \param map Resulting map
    template <typename Derived, int Rows, int Cols>
    static Eigen::Map<const Eigen::Array<Derived, Rows, Cols>> ArrayToMap(const Eigen::Array<Derived, Rows, Cols>& matrix)
    {
        return ArrayToMap(matrix.rows(), matrix.cols(), matrix.data());
    }

    template <typename Derived, int Rows, int Cols = 1, int Options>
    static Eigen::Map<const Eigen::Matrix<Derived, Rows, Cols, Options>> MatrixToMap(const int& rows,
                                                                                     const int& cols,
                                                                                     const Derived* values)
    {
        return Eigen::Map<const Eigen::Matrix<Derived, Rows, Cols, Options>>(values, rows, cols);
    }

    /// \brief Convert a matrix expressions to map
    /// \param matrix A const Eigen::Matrix
    /// \param map Resulting map
    template <typename Derived, int Rows, int Cols, int Options>
    static Eigen::Map<const Eigen::Matrix<Derived, Rows, Cols, Options>> MatrixToMap(const Eigen::Matrix<Derived, Rows, Cols, Options>& matrix)
    {
        return MatrixToMap(matrix.rows(), matrix.cols(), matrix.data());
    }

    /// \brief Convert a std::vector of const Eigen::Matrix expressions to const map
    /// \param matrices A std::vector of const matrices
    /// \param maps Resulting std::vector of map
    template <typename Derived, int Rows, int Cols>
    static std::vector<Eigen::Map<const Eigen::Matrix<Derived, Rows, Cols>>> MatrixToMap(const std::vector<Eigen::Matrix<Derived, Rows, Cols>>& matrices)
    {
        std::vector<Eigen::Map<const Eigen::Matrix<Derived, Rows, Cols>>> maps;
        maps.reserve(matrices.size());

        for (const Eigen::Matrix<Derived, Rows, Cols>& matrix : matrices)
            maps.push_back(MatrixToMap(matrix));
    }

    template <typename Derived>
    static std::vector<Eigen::Triplet<Derived>> SparseMatrixToTriplets(const Eigen::SparseMatrix<Derived>& sparseMatrix)
    {
        std::vector<Eigen::Triplet<Derived>> triplets;
        for(int i = 0; i < sparseMatrix.outerSize(); i++)
            for(typename Eigen::SparseMatrix<Derived>::InnerIterator it(sparseMatrix,i); it; ++it)
                triplets.emplace_back(it.row(),it.col(),it.value());
        return triplets;
    }

    static void ReadFromBinaryFile(const std::string& nameFile,
                                   Eigen::VectorXd& dataToRead,
                                   const unsigned int& dataSizeToRead = 0,
                                   const unsigned int& startingPosition = 0);
    static void WriteToBinaryFile(const std::string& nameFile,
                                  const Eigen::VectorXd& vec,
                                  const unsigned int& dataSizeToWrite = 0,
                                  const unsigned int& dataStartingPositionToWrite = 0,
                                  const bool& append = false);
    static void ReadFromBinaryFile(const std::string& nameFile,
                                   Eigen::Map<Eigen::VectorXd>& dataToRead,
                                   const unsigned int& dataSizeToRead = 0,
                                   const unsigned int& startingPosition = 0);
    static void WriteToBinaryFile(const std::string& nameFile,
                                  const Eigen::Map<const Eigen::VectorXd>& vec,
                                  const unsigned int& dataSizeToWrite = 0,
                                  const unsigned int& dataStartingPositionToWrite = 0,
                                  const bool& append = false);
    static void ReadFromBinaryFile(const std::string& nameFile,
                                   unsigned int& rows,
                                   unsigned int& cols,
                                   std::vector<Eigen::Triplet<double> >& triplets,
                                   const unsigned int& startingPosition = 0);
    static void ReadFromBinaryFile(const std::string& nameFile,
                                   Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix,
                                   const unsigned int& startingPosition = 0);
    static void ReadFromBinaryFile(const std::string& nameFile,
                                   Eigen::SparseMatrix<double, Eigen::ColMajor>& matrix,
                                   const unsigned int& startingPosition = 0);
    static void WriteToBinaryFile(const std::string& nameFile,
                                  const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix,
                                  const Eigen::UpLoType& matrixType = (Eigen::UpLoType)0,
                                  const bool& append = false);
    static void WriteToBinaryFile(const std::string& nameFile,
                                  const Eigen::SparseMatrix<double, Eigen::ColMajor>& matrix,
                                  const Eigen::UpLoType& matrixType = (Eigen::UpLoType)0,
                                  const bool& append = false);

    static void RemoveRow(Eigen::MatrixXd& matrix,
                          const unsigned int& rowToRemove);
    static void RemoveColumn(Eigen::MatrixXd& matrix,
                             const unsigned int& colToRemove);
};
}
#endif // __EIGEN_UTILITIES_H
