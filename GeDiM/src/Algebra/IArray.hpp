#ifndef __IArray_HPP
#define __IArray_HPP

#include <ostream>
#include <vector>
#include "ISparseArray.hpp"

namespace Gedim
{
  /// \brief Interface used for column double Array
  class IArray
  {
    public:
      virtual ~IArray() {}

      /// \brief Matrix create call. Last call to complete the matrix structure.
      virtual void Create() = 0;
      /// \brief Set the size the Array. Resize a Array for Eigen application.
      /// \param numCols the number of cols
      virtual void SetSize(const unsigned int& numCols) = 0;
      /// \brief Set the size the Array. Resize a Array for Eigen application.
      /// \param numCols the number of cols
      /// \param numLocalCols the local size for the Array (It can be set for each process. For PETSC applciation when not specified is set to PETSC_DECIDE)
      virtual void SetSizes(const unsigned int& numCols,
                            const unsigned int& numLocalCols) = 0;
      /// \return the size of the array
      virtual unsigned int Size() const = 0;
      /// \return the const pointer of the array
      virtual const double* Data() const = 0;
      /// \return the pointer of the array
      virtual double* Data() = 0;
      /// \brief Put zero-values in the Array
      virtual void Zeros() = 0;
      /// \brief Put one-values in the Array
      virtual void Ones() = 0;
      /// \brief operator []
      /// \param i the i_th position starting from 0
      /// \return the reference to i_th element
      virtual double& operator[] (const int& i) = 0;
      /// \brief const operator []
      /// \param i the i_th position starting from 0
      /// \return the const reference to i_th element
      virtual const double& operator[] (const int& i) const = 0;
      /// \brief Set the Array Value
      /// \param i the i_th global position starting from 0
      /// \param val element to set
      virtual void SetValue(const int& i,
                            const double& val) = 0;
      /// \brief Set the Array Values
      /// \param indices the indices
      /// \param values the values
      virtual void SetValues(const std::vector<int>& indices,
                             const std::vector<double>& values) = 0;
      /// \brief Add value to the Array
      /// \param i the i_th global position starting from 0
      /// \param val element to add
      virtual void AddValue(const int& i,
                            const double& val) = 0;
      /// \brief Add values in the indices to the Array
      /// \param indices the indices
      /// \param values the values
      virtual void AddValues(const std::vector<int>& indices,
                             const std::vector<double>& values) = 0;

      /// \brief Multiplication of matrix and vector, equivalent to v += A * w
      virtual void SumMultiplication(const ISparseArray& A,
                                     const IArray& w) = 0;

      /// \brief Multiplication of matrix and vector and sum to local vector, equivalent to v -= A * w
      virtual void SubtractionMultiplication(const ISparseArray& A,
                                             const IArray& w) = 0;

      /// \brief Concatenate vectors, equivalent u<< v, w
      virtual void Concatenate(const IArray& v,
                               const IArray& w) = 0;


      /// \return the \em l2 norm of array
      virtual double Norm() const = 0;
      /// \return the dot product
      virtual double Dot(const IArray& v) const = 0;

      /// \brief Copy the content of vector v inside this
      virtual void Copy(const IArray& v) = 0;

      /// \brief Print the array
      virtual std::ostream& Print(std::ostream& output) const = 0;

      /// \brief Write array to binary file
      /// \param filePath the file Path
      /// \param dataSizeToWrite the size of data to write, default 0 means all the vector
      /// \param dataStartingPositionToWrite the file starting position, default 0
      /// \param append true if append operation, default false
      virtual void ToBinaryFile(const std::string& filePath,
                                const unsigned int& dataSizeToWrite = 0,
                                const unsigned int& dataStartingPositionToWrite = 0,
                                const bool& append = false) const = 0;

      virtual IArray& operator+=(const IArray& v) = 0;
      virtual IArray& operator-=(const IArray& v) = 0;
      virtual IArray& operator*=(const double& c) = 0;
      virtual IArray& operator/=(const double& c) = 0;

      /** \addtogroup Operator overloadings
      *  @{ */

      /// \brief operator <<
      inline friend std::ostream& operator<<(std::ostream& output,
                                             const IArray& b)
      { return b.Print(output); }
      /** @}*/
  };
}

#endif // __IArray_HPP
