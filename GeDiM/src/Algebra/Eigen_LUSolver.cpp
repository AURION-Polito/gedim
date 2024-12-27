#include "Eigen_LUSolver.hpp"

#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"

using namespace Eigen;

namespace Gedim
{
// ***************************************************************************
template class Eigen_LUSolver<VectorXd, SparseMatrix<double>>;
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
Eigen_LUSolver<Eigen_ArrayType, Eigen_SparseArrayType>::Eigen_LUSolver()
{
    _rightHandSide = nullptr;
    _solution = nullptr;
}
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
Eigen_LUSolver<Eigen_ArrayType, Eigen_SparseArrayType>::~Eigen_LUSolver()
{
    _rightHandSide = nullptr;
    _solution = nullptr;
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_LUSolver<Eigen_ArrayType, Eigen_SparseArrayType>::Initialize(const ISparseArray &matrix,
                                                                        const IArray &rightHandSide,
                                                                        IArray &solution,
                                                                        const Configuration &)
{
    _rightHandSide = &rightHandSide;
    _solution = &solution;

    Initialize(matrix);
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_LUSolver<Eigen_ArrayType, Eigen_SparseArrayType>::Initialize(const IArray &rightHandSide, IArray &solution, const Configuration &)
{
    _rightHandSide = &rightHandSide;
    _solution = &solution;
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
ILinearSolver::SolutionInfo Eigen_LUSolver<Eigen_ArrayType, Eigen_SparseArrayType>::Solve() const
{
    if (_rightHandSide == nullptr || _solution == nullptr)
        throw std::runtime_error("No solver initialization found");

    return Solve(*_rightHandSide, *_solution);
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
void Eigen_LUSolver<Eigen_ArrayType, Eigen_SparseArrayType>::Initialize(const ISparseArray &matrix, const Configuration &)
{
    const SparseMatrix<double> &_matrix = static_cast<const Eigen_SparseArray<Eigen_SparseArrayType> &>(matrix);

    switch (matrix.Type())
    {
    case ISparseArray::SparseArrayTypes::Symmetric:
        linearSolver.compute(_matrix.selfadjointView<Lower>());
        break;
    case ISparseArray::SparseArrayTypes::None:
    case ISparseArray::SparseArrayTypes::Lower:
    case ISparseArray::SparseArrayTypes::Upper:
    case ISparseArray::SparseArrayTypes::Diagonal:
        linearSolver.compute(_matrix);
        break;
    default:
        throw std::runtime_error("Matrix type not supported");
    }

    if (linearSolver.info() != Eigen::ComputationInfo::Success)
        throw std::runtime_error("LU Factorization computation failed");
}
// ***************************************************************************
template <typename Eigen_ArrayType, typename Eigen_SparseArrayType>
ILinearSolver::SolutionInfo Eigen_LUSolver<Eigen_ArrayType, Eigen_SparseArrayType>::Solve(const IArray &rightHandSide,
                                                                                          IArray &solution) const
{
    const VectorXd &_rightHandSide = static_cast<const Eigen_Array<Eigen_ArrayType> &>(rightHandSide);
    VectorXd &_solution = static_cast<Eigen_Array<Eigen_ArrayType> &>(solution);

    _solution = linearSolver.solve(_rightHandSide);

    if (linearSolver.info() != Eigen::ComputationInfo::Success)
        throw std::runtime_error("LU solver failed");

    return {};
}
// ***************************************************************************
} // namespace Gedim
