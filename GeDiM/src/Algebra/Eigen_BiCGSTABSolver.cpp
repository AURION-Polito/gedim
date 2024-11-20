#include "Eigen_BiCGSTABSolver.hpp"

#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"

using namespace Eigen;
using namespace std;

namespace Gedim
{
  // ***************************************************************************
  template class Eigen_BiCGSTABSolver<VectorXd, SparseMatrix<double>, Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IdentityPreconditioner>>;
  template class Eigen_BiCGSTABSolver<VectorXd, SparseMatrix<double>, Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>>>;
  template class Eigen_BiCGSTABSolver<VectorXd, SparseMatrix<double>, Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>>>;
  // ***************************************************************************
  template<typename Eigen_ArrayType,
           typename Eigen_SparseArrayType,
           typename Eigen_SolverType>
  Eigen_BiCGSTABSolver<
  Eigen_ArrayType,
  Eigen_SparseArrayType,
  Eigen_SolverType>::Eigen_BiCGSTABSolver()
  {
    _rightHandSide = nullptr;
    _solution = nullptr;
  }
  template<typename Eigen_ArrayType,
           typename Eigen_SparseArrayType,
           typename Eigen_SolverType>
  Eigen_BiCGSTABSolver<
  Eigen_ArrayType,
  Eigen_SparseArrayType,
  Eigen_SolverType>::~Eigen_BiCGSTABSolver()
  {
    _rightHandSide = nullptr;
    _solution = nullptr;
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType,
           typename Eigen_SparseArrayType,
           typename Eigen_SolverType>
  void Eigen_BiCGSTABSolver<
  Eigen_ArrayType,
  Eigen_SparseArrayType,
  Eigen_SolverType>::Initialize(const ISparseArray& matrix,
                                const IArray& rightHandSide,
                                IArray& solution,
                                const Configuration& configuration)
  {
    _rightHandSide = &rightHandSide;
    _solution = &solution;

    Initialize(matrix,
               configuration);
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType,
           typename Eigen_SparseArrayType,
           typename Eigen_SolverType>
  void Eigen_BiCGSTABSolver<
  Eigen_ArrayType,
  Eigen_SparseArrayType,
  Eigen_SolverType>::Initialize(const IArray& rightHandSide,
                                IArray& solution,
                                const Configuration& configuration)
  {
    _rightHandSide = &rightHandSide;
    _solution = &solution;
    _config = configuration;
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType,
           typename Eigen_SparseArrayType,
           typename Eigen_SolverType>
  ILinearSolver::SolutionInfo Eigen_BiCGSTABSolver<
  Eigen_ArrayType,
  Eigen_SparseArrayType,
  Eigen_SolverType>::Solve() const
  {
    if (_rightHandSide == nullptr ||
        _solution == nullptr)
      throw runtime_error("No initialization found");

    return Solve(*_rightHandSide,
                 *_solution);
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType,
           typename Eigen_SparseArrayType,
           typename Eigen_SolverType>
  void Eigen_BiCGSTABSolver<
  Eigen_ArrayType,
  Eigen_SparseArrayType,
  Eigen_SolverType>::Initialize(const ISparseArray& matrix,
                                const Configuration& configuration)
  {
    _config = configuration;

    const SparseMatrix<double>& _matrix = static_cast<const Eigen_SparseArray<SparseMatrix<double>>&>(matrix);

    linearSolver.setTolerance(_config.Tolerance);
    linearSolver.setMaxIterations(_config.MaxIterations);

    switch (matrix.Type())
    {
      case ISparseArray::SparseArrayTypes::Symmetric:
        throw std::runtime_error("Matrix symmetric type not supported");
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
      throw runtime_error("BiCGSTAB initialization computation failed");
  }
  // ***************************************************************************
  template<typename Eigen_ArrayType,
           typename Eigen_SparseArrayType,
           typename Eigen_SolverType>
  ILinearSolver::SolutionInfo Eigen_BiCGSTABSolver<
  Eigen_ArrayType,
  Eigen_SparseArrayType,
  Eigen_SolverType>::Solve(const IArray& rightHandSide,
                           IArray& solution) const
  {
    const VectorXd& _rightHandSide = static_cast<const Eigen_Array<VectorXd>&>(rightHandSide);
    VectorXd& _solution = static_cast<Eigen_Array<VectorXd>&>(solution);

    _solution = linearSolver.solve(_rightHandSide);

    if (linearSolver.info() == Eigen::ComputationInfo::NumericalIssue ||
        linearSolver.info() == Eigen::ComputationInfo::InvalidInput)
      throw runtime_error("Cholesky solver failed");

    return
    {
      static_cast<unsigned int>(linearSolver.iterations()),
          linearSolver.error()
    };
  }
  // ***************************************************************************
}
