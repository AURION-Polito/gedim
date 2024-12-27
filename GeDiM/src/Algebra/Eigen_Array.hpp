#ifndef __Eigen_Array_HPP
#define __Eigen_Array_HPP

#include "Eigen/Eigen"
#include "IArray.hpp"

namespace Gedim
{
/// \brief Eigen column vector
template <typename Eigen_ArrayType = Eigen::VectorXd, typename Eigen_SparseArrayType = Eigen::SparseMatrix<double>>
class Eigen_Array final : public IArray
{
  private:
    Eigen_ArrayType _vector;

  public:
    Eigen_Array()
    {
    }
    Eigen_Array(const Eigen_ArrayType &vector)
    {
        _vector = vector;
    }
    Eigen_Array(Eigen_ArrayType &&vector)
    {
        _vector = std::move(vector);
    }
    ~Eigen_Array()
    {
    }

    operator Eigen_ArrayType &()
    {
        return _vector;
    }
    operator const Eigen_ArrayType &() const
    {
        return _vector;
    }
    inline Eigen_ArrayType &Cast(IArray &v)
    {
        return (Eigen_ArrayType &)static_cast<Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> &>(v);
    }
    inline const Eigen_ArrayType &Cast(const IArray &v) const
    {
        return (const Eigen_ArrayType &)static_cast<const Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> &>(v);
    }

    inline void Create()
    {
    }
    inline void Destroy()
    {
    }
    inline void SetSize(const unsigned int &numCols)
    {
        _vector.setZero(numCols);
    }
    inline void SetSizes(const unsigned int &numCols, const unsigned int & = 0)
    {
        return SetSize(numCols);
    }

    inline unsigned int Size() const
    {
        return _vector.size();
    }
    inline const double *Data() const
    {
        return _vector.data();
    }
    inline double *Data()
    {
        return _vector.data();
    }

    inline void SetValue(const int &i, const double &val)
    {
        _vector[i] = val;
    }
    void SetValues(const std::vector<int> &indices, const std::vector<double> &values);
    inline void AddValue(const int &i, const double &val)
    {
        _vector[i] += val;
    }
    void AddValues(const std::vector<int> &indices, const std::vector<double> &values);
    std::vector<double> GetValues(const std::vector<int> &indices = {}) const;

    inline double GetValue(const int &i) const
    {
        return _vector[i];
    }

    inline void Zeros()
    {
        _vector.setZero();
    }
    inline void Ones()
    {
        _vector.setOnes();
    }
    inline void Constant(const double &c)
    {
        _vector.setConstant(c);
    }

    void SumMultiplication(const ISparseArray &A, const IArray &w);
    void SubtractionMultiplication(const ISparseArray &A, const IArray &w);

    double &operator[](const int &i)
    {
        return _vector[i];
    }
    const double &operator[](const int &i) const
    {
        return _vector[i];
    }

    inline void Concatenate(const IArray &v, const IArray &w)
    {
        _vector << Cast(v), Cast(w);
    }

    inline std::ostream &Print(std::ostream &output) const
    {
        return output << _vector.transpose();
    }

    inline IArray &operator+=(const IArray &v)
    {
        _vector += Cast(v);
        return *this;
    }
    inline IArray &operator-=(const IArray &v)
    {
        _vector -= Cast(v);
        return *this;
    }
    inline IArray &operator*=(const double &c)
    {
        _vector *= c;
        return *this;
    }
    inline IArray &operator/=(const double &c)
    {
        _vector /= c;
        return *this;
    }

    Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> operator+(const Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> &other) const
    {
        Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> sumResult;
        sumResult._vector = _vector + other._vector;
        return sumResult;
    }

    Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> operator-(const Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> &other) const
    {
        Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> sumResult;
        sumResult._vector = _vector - other._vector;
        return sumResult;
    }

    Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> operator*(const double &other) const
    {
        Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> sumResult;
        sumResult._vector = other * _vector;
        return sumResult;
    }

    Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> &operator=(Eigen_ArrayType &&vector)
    {
        _vector = std::move(vector);
        return *this;
    }

    void ToBinaryFile(const std::string &filePath,
                      const unsigned int &dataSizeToWrite = 0,
                      const unsigned int &dataStartingPositionToWrite = 0,
                      const bool &append = false) const;

    inline double Norm() const
    {
        return _vector.norm();
    }

    inline double Dot(const IArray &v) const
    {
        return _vector.dot(Cast(v));
    }

    inline void Copy(const IArray &v)
    {
        _vector = Cast(v);
    }

    inline Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType> Segment(const unsigned int starter, const unsigned int numElements) const
    {
        return Eigen_Array<Eigen_ArrayType, Eigen_SparseArrayType>(_vector.segment(starter, numElements));
    }
};
} // namespace Gedim

#endif // __Eigen_Array_HPP
