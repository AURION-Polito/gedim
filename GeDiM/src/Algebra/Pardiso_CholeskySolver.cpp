#include "Pardiso_CholeskySolver.hpp"

#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"

#include "CommonUtilities.hpp"

using namespace Eigen;
using namespace std;

namespace Gedim
{
  // ***************************************************************************
  template class Pardiso_CholeskySolver<VectorXd, SparseMatrix<double>>;
  // ***************************************************************************
  template<typename Eigen_ArrayType, typename Eigen_SparseArrayType>
  Pardiso_CholeskySolver<Eigen_ArrayType, Eigen_SparseArrayType>::Pardiso_CholeskySolver()
  {
    _rightHandSide = nullptr;
    _solution = nullptr;
  }
  template<typename Eigen_ArrayType, typename Eigen_SparseArrayType>
  Pardiso_CholeskySolver<Eigen_ArrayType, Eigen_SparseArrayType>::~Pardiso_CholeskySolver()
  {
    _rightHandSide = nullptr;
    _solution = nullptr;
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType, typename Eigen_SparseArrayType>
  void Pardiso_CholeskySolver<Eigen_ArrayType, Eigen_SparseArrayType>::Initialize(const ISparseArray& matrix,
                                                                                  const IArray& rightHandSide,
                                                                                  IArray& solution,
                                                                                  const ILinearSolver::Configuration&)
  {
    _rightHandSide = &rightHandSide;
    _solution = &solution;

    Initialize(matrix);
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType, typename Eigen_SparseArrayType>
  ILinearSolver::SolutionInfo Pardiso_CholeskySolver<Eigen_ArrayType, Eigen_SparseArrayType>::Solve() const
  {
    if (_rightHandSide == nullptr ||
        _solution == nullptr)
      throw runtime_error("No initialization found");

    return Solve(*_rightHandSide,
          *_solution);
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType, typename Eigen_SparseArrayType>
  void Pardiso_CholeskySolver<Eigen_ArrayType, Eigen_SparseArrayType>::Initialize(const ISparseArray& matrix,
                                                                                  const ILinearSolver::Configuration&)
  {
#if ENABLE_MKL
    const SparseMatrix<double>& _matrix = static_cast<const Eigen_SparseArray<SparseMatrix<double>>&>(matrix);

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
#else
    Utilities::Unused(matrix);
#endif
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType, typename Eigen_SparseArrayType>
  ILinearSolver::SolutionInfo Pardiso_CholeskySolver<Eigen_ArrayType, Eigen_SparseArrayType>::Solve(const IArray& rightHandSide,
                                                                                                    IArray& solution) const
  {
#if ENABLE_MKL
    const VectorXd& _rightHandSide = static_cast<const Eigen_Array<VectorXd>&>(rightHandSide);
    VectorXd& _solution = static_cast<Eigen_Array<VectorXd>&>(solution);

    _solution = linearSolver.solve(_rightHandSide);

    if (linearSolver.info() == Eigen::ComputationInfo::NumericalIssue ||
        linearSolver.info() == Eigen::ComputationInfo::InvalidInput)
      throw runtime_error("Cholesky solver failed");
#else
    Utilities::Unused(rightHandSide);
    Utilities::Unused(solution);
#endif

    return {};
  }
  // ***************************************************************************
}
