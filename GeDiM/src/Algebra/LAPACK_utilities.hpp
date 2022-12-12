#ifndef LAPACK_UTILITIES_HPP
#define LAPACK_UTILITIES_HPP

#include "Eigen/Eigen"

namespace LAPACK_utilities
{
  /// Compute SVD: A = U * S * V'. It returns and computes only S and V'
  void svd(Eigen::MatrixXd A,
           Eigen::MatrixXd& V,
           Eigen::VectorXd& S);

  /// Compute SVD: A = U * S * V'. It returns and computes only S
  Eigen::VectorXd svd(Eigen::MatrixXd A);

  /// \brief Compute condition number in norm 2 given singular values
  /// \param s the singular values of matrix
  inline double cond(const Eigen::VectorXd& s)
  { return s[0] / s[s.size() - 1]; }

  /// Compute the modified Gram-Schmidt factorization of matrix X
  void MGS(const Eigen::MatrixXd& X,
           Eigen::MatrixXd& Q,
           Eigen::MatrixXd& R);

  /// Extract upper triangular part of matrix X
  Eigen::MatrixXd triu(const Eigen::MatrixXd& X,
                       const unsigned int& i);

  /// Compute eigenvalues and eigenvectors of A = R *  D * R'
  void eig(const Eigen::MatrixXd A,
           Eigen::VectorXd& D,
           Eigen::MatrixXd& R);

  /// Compute inverse of triangular matrix
  void inverseTri(const Eigen::MatrixXd A,
                  Eigen::MatrixXd& InvA,
                  const char& UPLO,
                  const char& DIAG);
}

#endif // LAPACK_UTILITIES_HPP
