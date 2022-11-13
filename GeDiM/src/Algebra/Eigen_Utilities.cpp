#include "Eigen_Utilities.hpp"
#include <iostream>

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  std::ostream& operator<<(std::ostream& out, const Eigen::VectorXd& vec)
  {
    const unsigned int sizeVector = vec.size();
    const unsigned int startIndexVector = 0;

    out<< "[";
    for (unsigned int i = startIndexVector; i < startIndexVector + sizeVector; i++)
      out<< (i != startIndexVector ? ", " : "")<< vec(i);
    out<< "]";

    return out;
  }
  // ***************************************************************************
  std::ostream& operator<<(std::ostream& out, const Eigen::Map<const Eigen::VectorXd>& vec)
  {
    const unsigned int sizeVector = vec.size();
    const unsigned int startIndexVector = 0;

    out<< "[";
    for (unsigned int i = startIndexVector; i < startIndexVector + sizeVector; i++)
      out<< (i != startIndexVector ? ", " : "")<< vec(i);
    out<< "]";

    return out;
  }
  // ***************************************************************************
  std::ostream& operator<<(std::ostream& out, const Eigen::MatrixXd& mat)
  {
    out<< "[[";
    for (unsigned int i = 0; i < mat.rows(); i++)
    {
      out<< (i != 0 ? "\n" : "");
      for (unsigned int j = 0; j < mat.cols(); j++)
        out<< (j != 0 ? ", " : "")<< mat(i, j);
    }
    out<< "]]";

    return out;
  }
  // ***************************************************************************
  void Eigen_Utilities::ReadFromBinaryFile(const string& nameFile,
                                           VectorXd& dataToRead,
                                           const unsigned int& dataSizeToRead,
                                           const unsigned int& startingPosition)
  {
    /// <ul>

    /// <li> Get file dimension
    unsigned int fileSize = 0;

    if (dataSizeToRead == 0)
    {
      Output::GetBinaryFileSize(nameFile,
                                fileSize,
                                sizeof(double),
                                startingPosition);
      if (fileSize == 0)
        return;
    }

    int fileContentSize = dataSizeToRead == 0 ? fileSize / sizeof(double) : dataSizeToRead;

    if (fileContentSize == 0)
      throw runtime_error("Error in file '" +
                          nameFile +
                          "': uncorrect content size. Obtained " +
                          to_string(fileContentSize));

    dataToRead.resize(fileContentSize);
    Output::ReadBinaryFile(nameFile,
                           dataToRead.data(),
                           sizeof(double),
                           fileContentSize,
                           startingPosition);

    /// </ul>
  }
  // ***************************************************************************
  void Eigen_Utilities::WriteToBinaryFile(const string& nameFile,
                                          const VectorXd& vec,
                                          const unsigned int& dataSizeToWrite,
                                          const unsigned int& dataStartingPositionToWrite,
                                          const bool& append)
  {
    /// <ul>

    if (dataSizeToWrite > vec.size())
      throw runtime_error("Error in write file '" +
                          nameFile +
                          "': uncorrect size. Required " +
                          to_string(dataSizeToWrite) +
                          ", given " +
                          to_string(vec.size()));

    if (dataStartingPositionToWrite > vec.size())
      throw runtime_error("Error in write file '" +
                          nameFile +
                          "': starting position. Required " +
                          to_string(dataStartingPositionToWrite) +
                          ", given " +
                          to_string(vec.size()));

    unsigned int fileDataSize = (dataSizeToWrite == 0) ? vec.size() : dataSizeToWrite;

    if (fileDataSize == 0)
      return;

    /// <li> Exporting file
    Output::WriteBinaryFile(nameFile,
                            vec.data() + dataStartingPositionToWrite,
                            sizeof(double),
                            fileDataSize,
                            append);

    /// </ul>
  }
  // ***************************************************************************
  void Eigen_Utilities::ReadFromBinaryFile(const string& nameFile,
                                           Map<VectorXd>& dataToRead,
                                           const unsigned int& dataSizeToRead,
                                           const unsigned int& startingPosition)
  {
    /// <ul>

    /// <li> Get file dimension
    unsigned int fileSize = 0;

    if (dataSizeToRead == 0)
    {
      Output::GetBinaryFileSize(nameFile,
                                fileSize,
                                sizeof(double),
                                startingPosition);

      if (fileSize == 0)
        return;
    }

    int fileContentSize = dataSizeToRead == 0 ? fileSize / sizeof(double) : dataSizeToRead;

    if (fileContentSize == 0)
      throw runtime_error("Error in file '" +
                          nameFile +
                          "': uncorrect content size. Obtained " +
                          to_string(fileContentSize));

    dataToRead.resize(fileContentSize);
    Output::ReadBinaryFile(nameFile,
                           dataToRead.data(),
                           sizeof(double),
                           fileContentSize,
                           startingPosition);
    /// </ul>
  }
  // ***************************************************************************
  void Eigen_Utilities::WriteToBinaryFile(const string& nameFile,
                                          const Map<const VectorXd>& vec,
                                          const unsigned int& dataSizeToWrite,
                                          const unsigned int& dataStartingPositionToWrite,
                                          const bool& append)
  {
    /// <ul>

    if (dataSizeToWrite > vec.size())
      throw runtime_error("Error in write file '" +
                          nameFile +
                          "': uncorrect size. Required " +
                          to_string(dataSizeToWrite) +
                          ", given " +
                          to_string(vec.size()));

    if (dataStartingPositionToWrite > vec.size())
      throw runtime_error("Error in write file '" +
                          nameFile +
                          "': starting position. Required " +
                          to_string(dataStartingPositionToWrite) +
                          ", given " +
                          to_string(vec.size()));

    unsigned int fileDataSize = (dataSizeToWrite == 0) ? vec.size() : dataSizeToWrite;

    if (fileDataSize == 0)
      return;

    /// <li> Exporting file
    Output::WriteBinaryFile(nameFile,
                            vec.data() + dataStartingPositionToWrite,
                            sizeof(double),
                            fileDataSize,
                            append);

    /// </ul>
  }
  // ***************************************************************************
  void Eigen_Utilities::ReadFromBinaryFile(const string& nameFile,
                                           unsigned int& rows,
                                           unsigned int& cols,
                                           vector< Triplet<double> >& triplets,
                                           const unsigned int& startingPosition)
  {
    vector<double> vectorToImport;

    if (!Output::ReadBinaryFile(nameFile, vectorToImport, 0, startingPosition))
      throw runtime_error("Error in read matrix from file '" +
                          nameFile +
                          "'");

    if (vectorToImport.size() < 3)
      throw runtime_error("Error in read matrix from file '" +
                          nameFile +
                          "': no enough information get from file");

    rows = vectorToImport[0];
    cols = vectorToImport[1];
    unsigned int nonZeros = vectorToImport[2];

    if (rows == 0)
      throw runtime_error("Error in read matrix from file '" +
                          nameFile +
                          "': no rows value");

    if (cols == 0)
      throw runtime_error("Error in read matrix from file '" +
                          nameFile +
                          "': no cols value");

    if (nonZeros == 0)
      throw runtime_error("Error in read matrix from file '" +
                          nameFile +
                          "': no nonZeros value");

    if (vectorToImport.size() - 3 < 3 * nonZeros)
      throw runtime_error("Error in read matrix from file '" +
                          nameFile +
                          "': no enough non zeros values");

    triplets.clear();
    triplets.reserve(nonZeros);

    for (unsigned int counter = 0; counter < 3 * nonZeros; counter += 3)
      triplets.push_back(Triplet<double>(vectorToImport[3 + counter], vectorToImport[3 + counter + 1], vectorToImport[3 + counter + 2]));
  }
  // ***************************************************************************
  void Eigen_Utilities::ReadFromBinaryFile(const string& nameFile,
                                           SparseMatrix<double, RowMajor>& matrix,
                                           const unsigned int& startingPosition)
  {
    unsigned int rows, cols;
    vector< Triplet<double> > triplets;

    ReadFromBinaryFile(nameFile,
                       rows,
                       cols,
                       triplets,
                       startingPosition);

    matrix.resize(rows, cols);

    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();
  }
  // ***************************************************************************
  void Eigen_Utilities::ReadFromBinaryFile(const string& nameFile,
                                           SparseMatrix<double, ColMajor>& matrix,
                                           const unsigned int& startingPosition)
  {
    unsigned int rows, cols;
    vector< Triplet<double> > triplets;

    ReadFromBinaryFile(nameFile,
                       rows,
                       cols,
                       triplets,
                       startingPosition);

    matrix.resize(rows, cols);

    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();
  }
  // ***************************************************************************
  void Eigen_Utilities::WriteToBinaryFile(const string& nameFile,
                                          const SparseMatrix<double, RowMajor>& matrix,
                                          const UpLoType& matrixType,
                                          const bool& append)
  {
    unsigned int nonZeros = matrix.nonZeros();

    if (nonZeros == 0)
      throw runtime_error("Error in write matrix to file '" +
                          nameFile +
                          "': empty matrix");

    vector<double> vectorToExport;
    unsigned int dimension = 3 + nonZeros * 3;
    vectorToExport.reserve(dimension);

    vectorToExport.push_back(matrix.rows());
    vectorToExport.push_back(matrix.cols());
    vectorToExport.push_back(0);

    for (int k = 0; k < matrix.outerSize(); ++k)
    {
      for (SparseMatrix<double, RowMajor>::InnerIterator it(matrix, k); it; ++it)
      {
        switch (matrixType)
        {
          case Lower:
            if (it.row() < it.col())
              continue;
            break;
          case Upper:
            if (it.row() > it.col())
              continue;
            break;
          default:
            break;
        }

        vectorToExport.push_back(it.row());
        vectorToExport.push_back(it.col());
        vectorToExport.push_back(it.value());
      }
    }
    vectorToExport.shrink_to_fit();
    dimension = vectorToExport.size();
    vectorToExport[2] = (dimension - 3) / 3;

    if (!Output::WriteBinaryFile(nameFile, vectorToExport, dimension, 0, append))
      throw runtime_error("");
  }
  // ***************************************************************************
  void Eigen_Utilities::WriteToBinaryFile(const string& nameFile,
                                          const SparseMatrix<double, ColMajor>& matrix,
                                          const UpLoType& matrixType,
                                          const bool& append)
  {
    unsigned int nonZeros = matrix.nonZeros();

    if (nonZeros == 0)
      throw runtime_error("Error in write matrix to file '" +
                          nameFile +
                          "': empty matrix");

    vector<double> vectorToExport;
    unsigned int dimension = 3 + nonZeros * 3;
    vectorToExport.reserve(dimension);

    vectorToExport.push_back(matrix.rows());
    vectorToExport.push_back(matrix.cols());
    vectorToExport.push_back(0);

    for (int k = 0; k < matrix.outerSize(); ++k)
    {
      for (SparseMatrix<double, ColMajor>::InnerIterator it(matrix, k); it; ++it)
      {
        switch (matrixType) {
          case Lower:
            if (it.row() < it.col())
              continue;
            break;
          case Upper:
            if (it.row() > it.col())
              continue;
            break;
          default:
            break;
        }

        vectorToExport.push_back(it.row());
        vectorToExport.push_back(it.col());
        vectorToExport.push_back(it.value());
      }
    }
    vectorToExport.shrink_to_fit();
    dimension = vectorToExport.size();
    vectorToExport[2] = (dimension - 3) / 3;

    if (!Output::WriteBinaryFile(nameFile, vectorToExport, dimension, 0, append))
      throw runtime_error("");
  }
  // ***************************************************************************
}
