#include "Eigen_CholeskySolver.hpp"

#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"

using namespace Eigen;
using namespace std;

namespace Gedim
{
// ***************************************************************************
template class Eigen_CholeskySolver<VectorXd, SparseMatrix<double>, Eigen::SimplicialLLT<SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>>>;
template class Eigen_CholeskySolver<VectorXd, SparseMatrix<double>, Eigen::SimplicialLDLT<SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>>>;
template class Eigen_CholeskySolver<VectorXd, SparseMatrix<double>, Eigen::SimplicialLLT<SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>>>;
template class Eigen_CholeskySolver<VectorXd, SparseMatrix<double>, Eigen::SimplicialLDLT<SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>>>;
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType, typename Eigen_SolverType>
void Eigen_CholeskySolver<Eigen_ArrayType, Eigen_SparseArrayType, Eigen_SolverType>::Initialize(const ISparseArray &matrix,
                                                                                                const IArray &rightHandSide,
                                                                                                IArray &solution,
                                                                                                const Configuration &)
{
    _rightHandSide = &rightHandSide;
    _solution = &solution;

    Initialize(matrix);
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType, typename Eigen_SolverType>
void Eigen_CholeskySolver<Eigen_ArrayType, Eigen_SparseArrayType, Eigen_SolverType>::Initialize(const IArray &rightHandSide,
                                                                                                IArray &solution,
                                                                                                const Configuration &)
{
    _rightHandSide = &rightHandSide;
    _solution = &solution;
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType, typename Eigen_SolverType>
ILinearSolver::SolutionInfo Eigen_CholeskySolver<Eigen_ArrayType, Eigen_SparseArrayType, Eigen_SolverType>::Solve() const
{
    if (_rightHandSide == nullptr || _solution == nullptr)
        throw runtime_error("No initialization found");

    return Solve(*_rightHandSide, *_solution);
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType, typename Eigen_SolverType>
void Eigen_CholeskySolver<Eigen_ArrayType, Eigen_SparseArrayType, Eigen_SolverType>::Initialize(const ISparseArray &matrix,
                                                                                                const Configuration &)
{
    const SparseMatrix<double> &_matrix = static_cast<const Eigen_SparseArray<SparseMatrix<double>> &>(matrix);

    switch (matrix.Type())
    {
    case ISparseArray::SparseArrayTypes::Symmetric:
        linearSolver.compute(_matrix);
        break;
    default:
        throw std::runtime_error("Matrix type not supported");
    }

    if (linearSolver.info() != Eigen::ComputationInfo::Success)
        throw runtime_error("Cholesky Factorization computation failed");
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType, typename Eigen_SolverType>
ILinearSolver::SolutionInfo Eigen_CholeskySolver<Eigen_ArrayType, Eigen_SparseArrayType, Eigen_SolverType>::Solve(const IArray &rightHandSide,
                                                                                                                  IArray &solution) const
{
    const VectorXd &_rightHandSide = static_cast<const Eigen_Array<VectorXd> &>(rightHandSide);
    VectorXd &_solution = static_cast<Eigen_Array<VectorXd> &>(solution);

    _solution = linearSolver.solve(_rightHandSide);

    if (linearSolver.info() == Eigen::ComputationInfo::NumericalIssue || linearSolver.info() == Eigen::ComputationInfo::InvalidInput)
        throw runtime_error("Cholesky solver failed");

    return {};
}
// ***************************************************************************
} // namespace Gedim
