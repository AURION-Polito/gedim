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

#include "LAPACK_utilities.hpp"

///  DGESVD computes the singular value decomposition (SVD) of a real M-by-N matrix A.
/// ( Double precision, simple driver)
extern "C" void dgesvd_(const char *JOBU,
                        const char *JOBVT,
                        const int *M,
                        const int *N,
                        double *A,
                        const int *LDA,
                        double *S,
                        double *U,
                        const int *LDU,
                        double *VT,
                        const int *LDVT,
                        double *WORK,
                        const int *LWORK,
                        int *INFO);

///  DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
/// the generalized eigenvalues, and optionally, the left and/or right
/// generalized eigenvectors.
extern "C" void dggev_(const char *JOBVL,
                       const char *JOBVR,
                       const int *N,
                       double *A,
                       const int *LDA,
                       double *B,
                       const int *LDB,
                       double *ALPHAR,
                       double *ALPHAI,
                       double *BETA,
                       double *VL,
                       const int *LDVL,
                       double *VR,
                       const int *LDVR,
                       double *WORK,
                       const int *LWORK,
                       int *INFO);

///  DGESVD computes the singular value decomposition (SVD) of a real M-by-N matrix A.
/// ( Double precision, Divide and conquer driver)
extern "C" void dgesdd_(const char *JOBZ,
                        const int *M,
                        const int *N,
                        double *A,
                        const int *LDA,
                        double *S,
                        double *U,
                        const int *LDU,
                        double *VT,
                        const int *LDVT,
                        double *WORK,
                        const int *LWORK,
                        int *IWORK,
                        int *INFO);

/// SSYEV computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A.
extern "C" void dsyev_(const char *JOBZ, const char *UPLO, const int *N, double *A, const int *LDA, double *D, double *WORK, const int *LWORK, int *INFO);

///  DTRTRI computes the inverse of a real upper or lower triangular matrix A.
extern "C" void dtrtri_(const char *UPLO, const char *DIAG, const int *N, double *A, const int *LDA, int *INFO);

extern "C" int
dgecon_(const char *norm, const int *n, double *a, const int *lda, const double *anorm, double *rcond, double *work, int *iwork, int *info, int len);

extern "C" int dgetrf_(const int *m, const int *n, double *a, const int *lda, int *lpiv, int *info);

extern "C" double dlange_(const char *norm, const int *m, const int *n, const double *a, const int *lda, double *work, const int norm_len);

///  DORGQR and DGEQRF computes the QR decomposition of a real M-by-N matrix A.
extern "C" void dgeqrf_(const int *m, const int *n, double *A, const int *lda, double *tau, double *work, int *lwork, int *info);
extern "C" void dorgqr_(const int *m, const int *n, const int *k, double *A, const int *lda, double *tau, double *work, int *lwork, int *info);

///  DGEQP3 computes the QR decomposition  with pivoting of a real M-by-N matrix A.
extern "C" void dgeqp3_(const int *m, const int *n, double *a, const int *lda, int *jpvt, double *tau, double *work, const int *lwork, int *info);

namespace LAPACK_utilities
{
// ***************************************************************************
void MGS(const Eigen::MatrixXd &X, Eigen::MatrixXd &Q, Eigen::MatrixXd &R)
{
    // Modified Gram-Schmidt.  [Q,R] = MGS(X);
    // G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
    const unsigned int m = X.rows();
    const unsigned int n = X.cols();

    Q = Eigen::MatrixXd::Zero(m, n);
    R = Eigen::MatrixXd::Zero(n, n);

    for (unsigned int i = 0; i < n; i++)
    {
        Q.col(i) = X.col(i);

        for (unsigned int j = 0; j < i; j++)
        {
            R(j, i) = Q.col(i).transpose() * Q.col(j);
            Q.col(i) = Q.col(i) - R(j, i) * Q.col(j);
        }

        R(i, i) = Q.col(i).norm();
        Q.col(i) = Q.col(i) / R(i, i);
    }
}
// ***************************************************************************
double rcondest(const Eigen::SparseMatrix<double> &sparseA)
{
    double rcond = 0.0;

    Eigen::MatrixXd A(sparseA);
    const int N = A.cols();
    const int LDA = A.rows();

    Eigen::VectorXi IWORK(N);
    Eigen::VectorXd WORK(4 * N);
    int INFO;

    double *Aptr = A.data();

    /* Computes the norm of x */
    double anorm = dlange_("1", &N, &N, Aptr, &LDA, WORK.data(), 1);

    /* Modifies x in place with a LU decomposition */
    dgetrf_(&N, &N, Aptr, &LDA, IWORK.data(), &INFO);
    if (INFO != 0)
        throw std::runtime_error("Error occurs in dgetrf");

    /* Computes the reciprocal norm */
    dgecon_("1", &N, Aptr, &LDA, &anorm, &rcond, WORK.data(), IWORK.data(), &INFO, 1);
    if (INFO != 0)
        throw std::runtime_error("Error occurs in dgecon");

    return rcond;
}
// ***************************************************************************
/// Extract the upper triangular matrix of matrix X. If i > 0, it returns the elements on and above the ith diagonal of
/// A.
Eigen::MatrixXd triu(const Eigen::MatrixXd &X, const unsigned int &i)
{
    unsigned int m = X.rows();
    unsigned int n = X.cols();

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, n);

    for (unsigned int yy = 0; yy < m - i; yy++)
    {
        for (unsigned int jj = i + yy; jj < n; jj++)
            A(yy, jj) = X(yy, jj);
    }

    return A;
}
// ***************************************************************************
/// Given A = U * S * V' returns only S and V'
void svd(Eigen::MatrixXd A, Eigen::MatrixXd &V, Eigen::VectorXd &S)
{
    char JOBU = 'N';  // all M columns of U are returned in array U
    char JOBVT = 'A'; // all N rows of V**T are returned in the array VT;

    int M = A.rows(); // The number of rows of the input matrix A.
    int N = A.cols(); // The number of columns of the input matrix A.

    V.resize(N, N);
    S.resize(std::min(M, N));

    const int LDA = M;
    const int LDU = 1;
    const int LDVT = N;

    double WORKDUMMY;
    int LWORK = -1; // Request optimum work size.
    int INFO = 0;

    dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), nullptr, &LDU, V.data(), &LDVT, &WORKDUMMY, &LWORK, &INFO);

    LWORK = int(WORKDUMMY) + 32;
    Eigen::VectorXd WORK(LWORK);

    dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), nullptr, &LDU, V.data(), &LDVT, WORK.data(), &LWORK, &INFO);

    if (INFO != 0)
        throw std::runtime_error("Error occurs in svd");

    WORK.resize(0);
}
// ***************************************************************************
void inverseTri(const Eigen::MatrixXd A, Eigen::MatrixXd &InvA, const char &UPLO, const char &DIAG)
{
    /// UPLO is CHARACTER*1 = 'U':  A is upper triangular or = 'L':  A is lower triangular.
    /// DIAG is CHARACTER*1 = 'N':  A is non-unit triangular; or = 'U':  A is unit triangular.

    int N = A.rows(); // The number of rows/columns of the input matrix A.
    int LDA = N;

    int INFO = 0;
    InvA = A;

    dtrtri_(&UPLO, &DIAG, &N, InvA.data(), &LDA, &INFO);

    if (INFO != 0)
        throw std::runtime_error("Error occurs in svd");
}
// ***************************************************************************
void eig(const Eigen::MatrixXd A, Eigen::VectorXd &D, Eigen::MatrixXd &R)
{
    char UPLO = 'U'; // Upper triangle of A is stored;
    char JOBZ = 'V'; // Compute eigenvalues and eigenvectors.

    int N = A.cols(); // The number of columns/rows of the input matrix A.

    R = LAPACK_utilities::triu(A, 0);

    Eigen::MatrixXd X = R;

    int LDA = N;

    D.resize(N);

    double WORKDUMMY;
    int LWORK = -1; // Request optimum work size.
    int INFO = 0;

    dsyev_(&JOBZ, &UPLO, &N, X.data(), &LDA, D.data(), &WORKDUMMY, &LWORK, &INFO);

    LWORK = int(WORKDUMMY) + 32;
    Eigen::VectorXd WORK(LWORK);

    dsyev_(&JOBZ, &UPLO, &N, R.data(), &LDA, D.data(), WORK.data(), &LWORK, &INFO);

    if (INFO != 0)
        throw std::runtime_error("Error occurs in eig");
}
// ***************************************************************************
/// Given A = U * S * V' returns only S
Eigen::VectorXd svd(Eigen::MatrixXd A)
{
    const char JOBU = 'N';  // no columns of U are returned in array U
    const char JOBVT = 'N'; // no rows of V**T are returned in the array VT;

    const int M = A.rows(); // The number of rows of the input matrix A.
    const int N = A.cols(); // The number of columns of the input matrix A.

    Eigen::VectorXd S(std::min(M, N));

    const int LDA = M;
    const int LDU = 1;
    const int LDVT = 1;

    double WORKDUMMY;
    int LWORK = -1; // Request optimum work size.
    int INFO = 0;

    dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), nullptr, &LDU, nullptr, &LDVT, &WORKDUMMY, &LWORK, &INFO);

    LWORK = int(WORKDUMMY) + 32;
    Eigen::VectorXd WORK(LWORK);

    dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), nullptr, &LDU, nullptr, &LDVT, WORK.data(), &LWORK, &INFO);

    if (INFO != 0)
        throw std::runtime_error("Error occurs in svd");

    return S;
}
// ***************************************************************************
/// Given A = U * S * V' returns U, S and V'
void svd(Eigen::MatrixXd A, Eigen::MatrixXd &U, Eigen::MatrixXd &V, Eigen::VectorXd &S)
{
    char JOBU = 'A';  // all M columns of U are returned in array U
    char JOBVT = 'A'; // all N rows of V**T are returned in the array VT;

    int M = A.rows(); // The number of rows of the input matrix A.
    int N = A.cols(); // The number of columns of the input matrix A.

    U.resize(M, M);
    V.resize(N, N);
    S.resize(std::min(M, N));

    int LDA = M;
    int LDU = M;
    int LDVT = N;

    double WORKDUMMY;
    int LWORK = -1; // Request optimum work size.
    int INFO = 0;

    dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, V.data(), &LDVT, &WORKDUMMY, &LWORK, &INFO);

    LWORK = int(WORKDUMMY) + 32;
    Eigen::VectorXd WORK(LWORK);

    dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, V.data(), &LDVT, WORK.data(), &LWORK, &INFO);

    if (INFO != 0)
        throw std::runtime_error("Error occurs in svd");
}
// ***************************************************************************
QR_Factorization MGS(const Eigen::MatrixXd &X, const double &tolerance)
{
    // Modified Gram-Schmidt.  [Q,R] = MGS(X);
    // X can have no max rank
    // G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
    const unsigned int m = X.rows();
    const unsigned int n = X.cols();

    QR_Factorization factorization = {Eigen::MatrixXd::Zero(m, n), Eigen::MatrixXd::Zero(n, n), Eigen::MatrixXd(), 0};

    Eigen::MatrixXd &Q = factorization.Q;
    Eigen::MatrixXd &R = factorization.R;

    for (unsigned int i = 0; i < n; i++)
    {
        Q.col(i) = X.col(i);

        for (unsigned int j = 0; j < i; j++)
        {
            R(j, i) = Q.col(i).transpose() * Q.col(j);
            Q.col(i) = Q.col(i) - R(j, i) * Q.col(j);
        }

        R(i, i) = Q.col(i).norm();

        if (std::abs(R(i, i)) <= tolerance)
            Q.col(i) = Eigen::VectorXd::Zero(m);
        else
        {
            Q.col(i) /= R(i, i);
            factorization.Space_Dimension++;
        }
    }

    return factorization;
}
// ***************************************************************************
QR_Factorization QR(const Eigen::MatrixXd &X, const double &tolerance)
{
    QR_Factorization qr_data;

    const unsigned int original_m = X.rows();
    const unsigned int original_n = X.cols();

    qr_data.R = X;

    // extend the matrix when rows > cols
    if (original_m > original_n)
    {
        qr_data.R.conservativeResize(original_m, original_m);
        qr_data.R.block(0, original_n, original_m, original_m - original_n).setZero();
    }

    const int m = qr_data.R.rows();
    const int n = qr_data.R.cols();

    std::vector<double> tau(std::min(m, n));
    std::vector<double> work(1);
    int lwork = -1;
    int info;

    // Compute workspace for dgeqrf
    dgeqrf_(&m, &n, qr_data.R.data(), &m, tau.data(), work.data(), &lwork, &info);
    lwork = static_cast<int>(work[0]);
    work.resize(lwork);

    // Compute QR factorization (X overwritten with R and reflectors)
    dgeqrf_(&m, &n, qr_data.R.data(), &m, tau.data(), work.data(), &lwork, &info);

    if (info != 0)
        throw std::runtime_error("Error occurs in dgeqrf");

    // Copy R to Q to build the Q matrix
    qr_data.Q = qr_data.R;
    qr_data.Q.conservativeResize(m, m);

    // Compute workspace for dorgqr
    work.resize(1);
    lwork = -1;
    dorgqr_(&m, &m, &m, qr_data.Q.data(), &m, tau.data(), work.data(), &lwork, &info);
    lwork = static_cast<int>(work[0]);
    work.resize(lwork);

    // Generate Q matrix
    dorgqr_(&m, &m, &m, qr_data.Q.data(), &m, tau.data(), work.data(), &lwork, &info);
    if (info != 0)
        throw std::runtime_error("Error occurs in dorgqr");

    qr_data.R.conservativeResize(original_m, original_n);
    qr_data.R = qr_data.R.triangularView<Eigen::Upper>();

    qr_data.Space_Dimension = rank(svd(qr_data.R), tolerance);

    return qr_data;
}
// ***************************************************************************
QR_Factorization QRP(const Eigen::MatrixXd &X, const double &tolerance)
{
    QR_Factorization qr_data;

    const unsigned int original_m = X.rows();
    const unsigned int original_n = X.cols();

    qr_data.R = X;

    // extend the matrix when rows > cols
    if (original_m > original_n)
    {
        qr_data.R.conservativeResize(original_m, original_m);
        qr_data.R.block(0, original_n, original_m, original_m - original_n).setZero();
    }

    const int m = qr_data.R.rows();
    const int n = qr_data.R.cols();

    std::vector<int> jpvt(n, 0);
    std::vector<double> tau(std::min(m, n));
    std::vector<double> work(1);
    int lwork = -1;
    int info;

    // Compute workspace for dgeqp3_
    dgeqp3_(&m, &n, qr_data.R.data(), &m, jpvt.data(), tau.data(), work.data(), &lwork, &info);
    lwork = static_cast<int>(work[0]);
    work.resize(lwork);

    // Compute QR factorization with pivoting (X overwritten with R and reflectors)
    dgeqp3_(&m, &n, qr_data.R.data(), &m, jpvt.data(), tau.data(), work.data(), &lwork, &info);

    if (info != 0)
        throw std::runtime_error("Error occurs in dgeqp3_");

    // Copy R to Q to build the Q matrix
    qr_data.Q = qr_data.R;
    qr_data.Q.conservativeResize(m, m);

    // Compute workspace for dorgqr
    work.resize(1);
    lwork = -1;
    dorgqr_(&m, &m, &m, qr_data.Q.data(), &m, tau.data(), work.data(), &lwork, &info);
    lwork = static_cast<int>(work[0]);
    work.resize(lwork);

    // Generate Q matrix
    dorgqr_(&m, &m, &m, qr_data.Q.data(), &m, tau.data(), work.data(), &lwork, &info);
    if (info != 0)
        throw std::runtime_error("Error occurs in dorgqr");

    qr_data.R.conservativeResize(original_m, original_n);
    qr_data.R = qr_data.R.triangularView<Eigen::Upper>();

    qr_data.Space_Dimension = 0;
    unsigned int min_dim = std::min(original_m, original_n);
    for (unsigned int i = 0; i < min_dim; ++i)
    {
        if (std::abs(qr_data.R(i, i)) <= tolerance)
            break;

        qr_data.Space_Dimension++;
    }

    qr_data.P.setZero(original_n, original_n);
    for (unsigned int i = 0; i < original_n; ++i)
        qr_data.P(jpvt[i] - 1, i) = 1.0;

    return qr_data;
}
// ***************************************************************************
unsigned int rank(const Eigen::VectorXd &s, const double &tolerance)
{
    unsigned int rank = 0;
    for (unsigned int i = 0; i < s.size(); i++)
    {
        if (std::abs(s[i]) <= tolerance)
            break;

        rank++;
    }

    return rank;
}
// ***************************************************************************
} // namespace LAPACK_utilities
