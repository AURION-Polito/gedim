#ifndef __ISPARSEArray_HPP
#define __ISPARSEArray_HPP

#include "Eigen/Eigen"
#include <ostream>
#include <vector>

namespace Gedim
{
  class IArray;

  /// \brief Interface used for sparse Array of double
  class ISparseArray
  {
    public:
      enum class SparseArrayTypes
      {
        None = 0,
        Lower = 1,
        Upper = 2,
        Symmetric = 3,
        Diagonal = 4
      };

    public:
      virtual ~ISparseArray() {}

      /// \brief Resize the Matrix
      /// \param numRows the number of rows
      /// \param numCols the number of cols
      /// \param type the type of the matrix, default none
      virtual void SetSize(const unsigned int& numRows,
                           const unsigned int& numCols,
                           const SparseArrayTypes& type = SparseArrayTypes::None) = 0;

      /// \return The matrix type
      virtual SparseArrayTypes Type() const = 0;

      /// \brief Matrix Allocation is implemented in child classes
      /// \brief Matrix create call. Last call to complete the Array structure.
      virtual void Create() = 0;
      /// \brief Intermediate flush of the Matrix during creation.
      virtual void Flush() = 0;

      /// \brief Put zero-values in the Matrix
      virtual void Reset() = 0;
      /// \brief Given the row index i and the column index j the value val is put into the Matrix in ADD or INSERT mode
      virtual void Triplet(const unsigned int& i,
                           const unsigned int& j,
                           const double& val) = 0;
      /// \brief Given the row indices i and the column indices j the values val are put into the Matrix in ADD or INSERT mode
      virtual void Triplets(const std::vector<unsigned int>& i,
                            const std::vector<unsigned int>& j,
                            const std::vector<double>& val) = 0;

      /// \return The Frobenius norm of the sparse array
      virtual double Norm() const = 0;

      /// \brief Print the array
      virtual std::ostream& Print(std::ostream& output) const = 0;

      /// \brief Write sparse array to binary file
      /// \param filePath the file Path
      /// \param append true if append operation, default false
      virtual void ToBinaryFile(const std::string& filePath,
                                const bool& append = false) const = 0;

      /// \return the condition number of the sparse array in norm 2
      virtual double Cond(unsigned int condType = 0) const = 0;
      /// \return the number of non-zeros element
      virtual unsigned int NonZeros() const = 0;

      virtual ISparseArray& operator+=(const ISparseArray& A) = 0;
      virtual ISparseArray& operator-=(const ISparseArray& A) = 0;
      virtual ISparseArray& operator*=(const double& c) = 0;
      virtual ISparseArray& operator/=(const double& c) = 0;

      /** \addtogroup Operator overloadings
      *  @{ */

      /// \brief operator <<
      inline friend std::ostream& operator<<(std::ostream& output,
                                             const ISparseArray& b)
      { return b.Print(output); }
      /** @}*/
  };
}

#endif // __ISPARSEArray_HPP
