// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef LAPACK_UTILITIES_HPP
#define LAPACK_UTILITIES_HPP

#include "Eigen/Eigen"

namespace LAPACK_utilities
{
struct QR_Factorization
{
    Eigen::MatrixXd Q;
    Eigen::MatrixXd R;
    Eigen::MatrixXd P;
    unsigned int Space_Dimension;
};

/// Compute SVD: A = U * S * V'. It returns U, S and V'
void svd(Eigen::MatrixXd A, Eigen::MatrixXd &U, Eigen::MatrixXd &V, Eigen::VectorXd &S);

/// Compute SVD: A = U * S * V'. It returns and computes only S and V'
void svd(Eigen::MatrixXd A, Eigen::MatrixXd &V, Eigen::VectorXd &S);

/// Compute SVD: A = U * S * V'. It returns and computes only S
Eigen::VectorXd svd(Eigen::MatrixXd A);

/// \brief Compute condition number in norm 2 given singular values
/// \param s the singular values of matrix
inline double cond(const Eigen::VectorXd &s)
{
    return s[0] / s[s.size() - 1];
}

/// \param s the singular values of matrix
unsigned int rank(const Eigen::VectorXd &s, const double &tolerance = std::numeric_limits<double>::epsilon());

/// Compute the modified Gram-Schmidt factorization of matrix X
void MGS(const Eigen::MatrixXd &X, Eigen::MatrixXd &Q, Eigen::MatrixXd &R);

/// Compute the modified Gram-Schmidt factorization of matrix X with no maximum rank
QR_Factorization MGS(const Eigen::MatrixXd &X, const double &tolerance = std::numeric_limits<double>::epsilon());

/// Compute the QR matrix based on Householder reflectors
QR_Factorization QR(const Eigen::MatrixXd &X, const double &tolerance = std::numeric_limits<double>::epsilon());

/// Compute the QR with pivoting matrix based on Householder reflectors
QR_Factorization QRP(const Eigen::MatrixXd &X, const double &tolerance = std::numeric_limits<double>::epsilon());

/// Extract upper triangular part of matrix X
Eigen::MatrixXd triu(const Eigen::MatrixXd &X, const unsigned int &i);

/// Compute eigenvalues and eigenvectors of A = R *  D * R'
void eig(const Eigen::MatrixXd A, Eigen::VectorXd &D, Eigen::MatrixXd &R);

/// Compute inverse of triangular matrix
void inverseTri(const Eigen::MatrixXd A, Eigen::MatrixXd &InvA, const char &UPLO, const char &DIAG);

double rcondest(const Eigen::SparseMatrix<double> &sparseA);
} // namespace LAPACK_utilities

#endif // LAPACK_UTILITIES_HPP
