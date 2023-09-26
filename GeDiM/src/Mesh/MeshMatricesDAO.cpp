#include "MeshMatricesDAO.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  template<typename T>
  void MeshMatricesDAO::InitializeNuberVectorWithConstantElements(std::vector<unsigned int>& numberElementVector,
                                                                  std::vector<T>& elementVector,
                                                                  const unsigned int numberElementSize,
                                                                  const unsigned int numberElements,
                                                                  const T& elementInitialization)
  {
    numberElementVector.resize(numberElementSize + 1, 0);
    for (unsigned int e = 0; e < numberElementSize; e++)
      numberElementVector[e + 1] = numberElementVector[e] + numberElements;

    elementVector.resize(numberElementSize * numberElements,
                         elementInitialization);
  }
  // ***************************************************************************
  template<typename T>
  void MeshMatricesDAO::InitializeNumberVector(std::vector<unsigned int>& numberElementVector,
                                               std::vector<T>& elementVector,
                                               const std::vector<unsigned int>& numberElements,
                                               const T& elementInitialization)
  {
    const unsigned int numberElementSize = numberElements.size();
    unsigned int totalSize = 0;
    numberElementVector.resize(numberElementSize + 1, 0);
    for (unsigned int e = 0; e < numberElementSize; e++)
    {
      totalSize += numberElements[e];
      numberElementVector[e + 1] = numberElementVector[e] + numberElements[e];
    }

    elementVector.resize(totalSize,
                         elementInitialization);
  }
  // ***************************************************************************
  template<typename T>
  void MeshMatricesDAO::ResizeNumberVectorWithNewNumberElements(vector<unsigned int>& numberElementVector,
                                                                vector<T>& elementVector,
                                                                const unsigned int& numberElements,
                                                                const unsigned int& vectorIndex,
                                                                const unsigned int& newNumberElements,
                                                                const T& newElementInitialization)
  {
    int numOriginalElements = (numberElementVector[vectorIndex + 1] -
                              numberElementVector[vectorIndex]);
    int newElements = newNumberElements -
                      numOriginalElements;
    unsigned int numOldElements = numberElementVector[numberElements];

    for (unsigned int next = vectorIndex; next < numberElements; next++)
      numberElementVector[next + 1] = numberElementVector[next + 1] +
                                      newElements;

    if (newElements > 0)
    {
      elementVector.resize(numberElementVector[numberElements], newElementInitialization);
      std::rotate(elementVector.begin() +
                  numberElementVector[vectorIndex] +
                  numOriginalElements,
                  elementVector.begin() +
                  numOldElements,
                  elementVector.end());
    }
    else
      elementVector.erase(elementVector.begin() +
                          numberElementVector[vectorIndex] +
                          newNumberElements,
                          elementVector.begin() +
                          numberElementVector[vectorIndex] +
                          numOriginalElements);
  }
  // ***************************************************************************
  template<class Container, class T>
  void MeshMatricesDAO::AlignMapContainerHigherElements(unordered_map<unsigned int, Container>& elements,
                                                        const T& minElement,
                                                        const T& newElementInitialization)
  {
    unordered_map<unsigned int, Container> tempMap;
    list<T> elementsToRemove;
    for (auto it = elements.begin(); it != elements.end(); it++)
    {
      if (it->first == minElement)
        elementsToRemove.push_back(it->first);

      std::vector<T> v(it->second.begin(), it->second.end());
      AlignContainerHigherElements(v,
                                   minElement,
                                   newElementInitialization);
      it->second.clear();
      for (unsigned int i = 0; i < v.size(); i++)
      {
        if (v[i] != newElementInitialization)
          it->second.insert(v[i]);
      }

      if (it->second.empty())
        elementsToRemove.push_back(it->first);
      else if (it->first > minElement)
      {
        tempMap.insert(pair<unsigned int, Container>(it->first - 1, it->second));
        elementsToRemove.push_back(it->first);
      }
    }

    for (const unsigned int& index : elementsToRemove)
      elements.erase(index);

    for (auto it = tempMap.begin(); it != tempMap.end(); it++)
      elements.insert(pair<unsigned int, Container>(it->first, it->second));
  }
  // ***************************************************************************
  template<class Container, class T>
  void MeshMatricesDAO::AlignContainerHigherElements(Container& elements,
                                                     const T& minElement,
                                                     const T& newElementInitialization)
  {
    for (auto it = elements.begin(); it != elements.end(); it++)
    {
      if (*it > minElement)
        (*it)--;
      else if ((*it) == minElement)
        (*it) = newElementInitialization;
    }
  }
  // ***************************************************************************
  template<class Container, class T>
  void MeshMatricesDAO::AlignContainerElements(Container& elements,
                                               const T& element,
                                               const T& newElementInitialization)
  {
    for (auto it = elements.begin(); it != elements.end(); it++)
    {
      if ((*it) == element)
        (*it) = newElementInitialization;
    }
  }
  // ***************************************************************************
  template<typename T>
  void MeshMatricesDAO::AlignSparseMatrixHigherElements(Eigen::SparseMatrix<T>& matrix,
                                                        const T& minElement)
  {
    list<Eigen::Triplet<T>> triplets;
    for (unsigned int k = 0; k < matrix.outerSize(); k++)
    {
      for (typename SparseMatrix<T>::InnerIterator it(matrix, k); it; ++it)
      {
        if (it.value() > minElement)
          triplets.push_back(Eigen::Triplet<T>(it.row(), it.col(), it.value() - 1));
        else if (it.value() < minElement)
          triplets.push_back(Eigen::Triplet<T>(it.row(), it.col(), it.value()));
      }
    }

    matrix.setZero();
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();
    triplets.clear();
  }
  // ***************************************************************************
  MeshMatricesDAO::MeshMatricesDAO(MeshMatrices& mesh) :
    _mesh(mesh)
  {
  }
  MeshMatricesDAO::~MeshMatricesDAO()
  {
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DsInitialize(const unsigned int& numberCell0Ds)
  {
    _mesh.NumberCell0D = numberCell0Ds;
    _mesh.Cell0DCoordinates.resize(3 * _mesh.NumberCell0D, 0.0);
    _mesh.Cell0DMarkers.resize(_mesh.NumberCell0D, 0);
    _mesh.ActiveCell0D.resize(_mesh.NumberCell0D, false);
    _mesh.NumberCell0DNeighbourCell1D.resize(_mesh.NumberCell0D + 1, 0);
    _mesh.NumberCell0DNeighbourCell2D.resize(_mesh.NumberCell0D + 1, 0);
    for (unsigned int p = 0; p < Cell0DNumberDoubleProperties(); p++)
      _mesh.Cell0DDoublePropertySizes[p].resize(_mesh.NumberCell0D + 1, 0);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell0DAppend(const unsigned int& numberCell0Ds)
  {
    unsigned int oldNumberCell0Ds = _mesh.NumberCell0D;

    _mesh.NumberCell0D = oldNumberCell0Ds + numberCell0Ds;
    for (unsigned int nc = 0; nc < numberCell0Ds; nc++)
    {
      for (unsigned int v = 0; v < 3; v++)
        _mesh.Cell0DCoordinates.push_back(0.0);

      _mesh.Cell0DMarkers.push_back(0);
      _mesh.ActiveCell0D.push_back(false);
      _mesh.NumberCell0DNeighbourCell1D.push_back(0);
      _mesh.NumberCell0DNeighbourCell2D.push_back(0);
      for (unsigned int p = 0; p < Cell0DNumberDoubleProperties(); p++)
        _mesh.Cell0DDoublePropertySizes[p].push_back(0);
    }

    for (unsigned c = 0; c < numberCell0Ds; c++)
      _mesh.NumberCell0DNeighbourCell1D[oldNumberCell0Ds + c + 1] = _mesh.NumberCell0DNeighbourCell1D[oldNumberCell0Ds];

    for (unsigned c = 0; c < numberCell0Ds; c++)
      _mesh.NumberCell0DNeighbourCell2D[oldNumberCell0Ds + c + 1] = _mesh.NumberCell0DNeighbourCell2D[oldNumberCell0Ds];


    for (unsigned int p = 0; p < Cell0DNumberDoubleProperties(); p++)
    {
      for (unsigned c = 0; c < numberCell0Ds; c++)
        _mesh.Cell0DDoublePropertySizes[p][oldNumberCell0Ds + c + 1] = _mesh.Cell0DDoublePropertySizes[p][oldNumberCell0Ds];
    }

    return oldNumberCell0Ds;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DRemove(const unsigned int& cell0DIndex)
  {
    Output::Assert(cell0DIndex < Cell0DTotalNumber());

    for (unsigned int p = 0; p < Cell0DNumberDoubleProperties(); p++)
    {
      ResizeNumberVectorWithNewNumberElements(_mesh.Cell0DDoublePropertySizes[p],
                                              _mesh.Cell0DDoublePropertyValues[p],
                                              _mesh.NumberCell0D,
                                              cell0DIndex,
                                              0);
      _mesh.Cell0DDoublePropertySizes[p].erase(std::next(_mesh.Cell0DDoublePropertySizes[p].begin(),
                                                         cell0DIndex));
    }

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell0DNeighbourCell1D,
                                            _mesh.Cell0DNeighbourCell1Ds,
                                            _mesh.NumberCell0D,
                                            cell0DIndex,
                                            0);
    _mesh.NumberCell0DNeighbourCell1D.erase(std::next(_mesh.NumberCell0DNeighbourCell1D.begin(),
                                                      cell0DIndex));

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell0DNeighbourCell2D,
                                            _mesh.Cell0DNeighbourCell2Ds,
                                            _mesh.NumberCell0D,
                                            cell0DIndex,
                                            0);
    _mesh.NumberCell0DNeighbourCell2D.erase(std::next(_mesh.NumberCell0DNeighbourCell2D.begin(),
                                                      cell0DIndex));



    _mesh.UpdatedCell0Ds.erase(cell0DIndex);

    _mesh.Cell0DCoordinates.erase(std::next(_mesh.Cell0DCoordinates.begin(), 3 * cell0DIndex),
                                  std::next(_mesh.Cell0DCoordinates.begin(), 3 * cell0DIndex + 3));
    _mesh.Cell0DMarkers.erase(std::next(_mesh.Cell0DMarkers.begin(), cell0DIndex));
    _mesh.ActiveCell0D.erase(std::next(_mesh.ActiveCell0D.begin(), cell0DIndex));
    _mesh.NumberCell0D--;

    AlignContainerHigherElements(_mesh.Cell1DVertices,
                                 cell0DIndex,
                                 _mesh.NumberCell0D);
    AlignContainerHigherElements(_mesh.Cell2DVertices,
                                 cell0DIndex,
                                 _mesh.NumberCell0D);
    AlignContainerHigherElements(_mesh.Cell3DVertices,
                                 cell0DIndex,
                                 _mesh.NumberCell0D);

    AlignMapContainerHigherElements(_mesh.UpdatedCell0Ds,
                                    cell0DIndex,
                                    _mesh.NumberCell0D);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DInsertCoordinates(const unsigned int& cell0DIndex,
                                                const Vector3d& coordinates)
  {
    Output::Assert(cell0DIndex < Cell0DTotalNumber());

    _mesh.Cell0DCoordinates[3 * cell0DIndex] = coordinates.x();
    _mesh.Cell0DCoordinates[3 * cell0DIndex + 1] = coordinates.y();
    _mesh.Cell0DCoordinates[3 * cell0DIndex + 2] = coordinates.z();
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DsInsertCoordinates(const MatrixXd& coordinates)
  {
    Output::Assert(coordinates.rows() == 3 && coordinates.cols() == Cell0DTotalNumber());

    for (unsigned int c = 0; c < Cell0DTotalNumber(); c++)
      Cell0DInsertCoordinates(c, coordinates.col(c));
  }
  // ***************************************************************************
  MatrixXd MeshMatricesDAO::Cell0DsCoordinates() const
  {
    MatrixXd coordinates(3, Cell0DTotalNumber());
    for (unsigned int v = 0; v < Cell0DTotalNumber(); v++)
      coordinates.col(v) << Cell0DCoordinates(v);
    return coordinates;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DInsertUpdatedCell0D(const unsigned int& cell0DIndex,
                                                  const unsigned int& updatedCell0DIdex)
  {
    Output::Assert(cell0DIndex < Cell0DTotalNumber());
    Output::Assert(updatedCell0DIdex < Cell0DTotalNumber());

    if (!Cell0DHasUpdatedCell0Ds(cell0DIndex))
      _mesh.UpdatedCell0Ds.insert(pair<unsigned int, unordered_set<unsigned int>>(cell0DIndex, {}));
    _mesh.UpdatedCell0Ds.at(cell0DIndex).insert(updatedCell0DIdex);
  }
  // ***************************************************************************
  bool MeshMatricesDAO::Cell0DUpdatedCell0Ds(const unsigned int& cell0DIndex,
                                             list<unsigned int>& updatedCell0DIds) const
  {
    Output::Assert(cell0DIndex < Cell0DTotalNumber());
    unordered_map<unsigned int, unordered_set<unsigned int>>::const_iterator iter = _mesh.UpdatedCell0Ds.find(cell0DIndex);

    if (iter == _mesh.UpdatedCell0Ds.end())
    {
      updatedCell0DIds.push_back(cell0DIndex);
      return true;
    }

    for (const unsigned int& childId : iter->second)
    {
      Cell0DUpdatedCell0Ds(childId,
                           updatedCell0DIds);
    }

    return false;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DsInitializeNeighbourCell1Ds(const std::vector<unsigned int>& numberNeighbourCell1Ds)
  {
    Output::Assert(numberNeighbourCell1Ds.size() == Cell0DTotalNumber());
    InitializeNumberVector(_mesh.NumberCell0DNeighbourCell1D,
                           _mesh.Cell0DNeighbourCell1Ds,
                           numberNeighbourCell1Ds,
                           std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DInitializeNeighbourCell1Ds(const unsigned int& cell0DIndex,
                                                         const unsigned int& numberNeighbourCell1Ds)
  {
    Output::Assert(cell0DIndex < Cell0DTotalNumber());

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell0DNeighbourCell1D,
                                            _mesh.Cell0DNeighbourCell1Ds,
                                            _mesh.NumberCell0D,
                                            cell0DIndex,
                                            numberNeighbourCell1Ds,
                                            std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DsInitializeNeighbourCell2Ds(const std::vector<unsigned int>& numberNeighbourCell2Ds)
  {
    Output::Assert(numberNeighbourCell2Ds.size() == Cell0DTotalNumber());
    InitializeNumberVector(_mesh.NumberCell0DNeighbourCell2D,
                           _mesh.Cell0DNeighbourCell2Ds,
                           numberNeighbourCell2Ds,
                           std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DInitializeNeighbourCell2Ds(const unsigned int& cell0DIndex, const unsigned int& numberNeighbourCell2Ds)
  {
    Output::Assert(cell0DIndex < Cell0DTotalNumber());

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell0DNeighbourCell2D,
                                            _mesh.Cell0DNeighbourCell2Ds,
                                            _mesh.NumberCell0D,
                                            cell0DIndex,
                                            numberNeighbourCell2Ds,
                                            std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DInitializeDoubleProperties(const unsigned int& numberDoubleProperties)
  {
    _mesh.Cell0DDoublePropertyIds.reserve(numberDoubleProperties);
    _mesh.Cell0DDoublePropertySizes.reserve(numberDoubleProperties);
    _mesh.Cell0DDoublePropertyValues.reserve(numberDoubleProperties);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell0DAddDoubleProperty(const string& propertyId)
  {
    Output::Assert(!Cell0DDoublePropertyExists(propertyId));
    unsigned int propertyIndex = _mesh.Cell0DDoublePropertySizes.size();

    _mesh.Cell0DDoublePropertyIndices.insert(pair<string, unsigned int>(propertyId,
                                                                        propertyIndex));
    _mesh.Cell0DDoublePropertyIds.push_back(propertyId);
    _mesh.Cell0DDoublePropertySizes.push_back(vector<unsigned int>(_mesh.NumberCell0D + 1, 0));
    _mesh.Cell0DDoublePropertyValues.push_back(vector<double>());
    return propertyIndex;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DsInitializeDoublePropertyValues(const unsigned int& propertyIndex,
                                                              const std::vector<unsigned int>& propertySizes)
  {
    Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());
    Output::Assert(propertySizes.size() == Cell0DTotalNumber());

    InitializeNumberVector(_mesh.Cell0DDoublePropertySizes[propertyIndex],
                           _mesh.Cell0DDoublePropertyValues[propertyIndex],
                           propertySizes,
                           0.0);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell0DInitializeDoublePropertyValues(const unsigned int& cell0DIndex,
                                                             const unsigned int& propertyIndex,
                                                             const unsigned int& propertySize)
  {
    Output::Assert(cell0DIndex < Cell0DTotalNumber());
    Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());

    ResizeNumberVectorWithNewNumberElements(_mesh.Cell0DDoublePropertySizes[propertyIndex],
                                            _mesh.Cell0DDoublePropertyValues[propertyIndex],
                                            _mesh.NumberCell0D,
                                            cell0DIndex,
                                            propertySize,
                                            0.0);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DsInitialize(const unsigned int& numberCell1Ds)
  {
    _mesh.NumberCell1D = numberCell1Ds;
    _mesh.Cell1DVertices.resize(2 * _mesh.NumberCell1D, 0);
    _mesh.Cell1DMarkers.resize(_mesh.NumberCell1D, 0);
    _mesh.ActiveCell1D.resize(_mesh.NumberCell1D, false);
    _mesh.Cell1DOriginalCell1Ds.resize(_mesh.NumberCell1D, std::numeric_limits<unsigned int>::max());
    _mesh.NumberCell1DNeighbourCell2D.resize(_mesh.NumberCell1D + 1, 0);
    for (unsigned int p = 0; p < Cell1DNumberDoubleProperties(); p++)
      _mesh.Cell1DDoublePropertySizes[p].resize(_mesh.NumberCell1D + 1, 0);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell1DAppend(const unsigned int& numberCell1Ds)
  {
    unsigned int oldNumberCell1Ds = _mesh.NumberCell1D;

    _mesh.NumberCell1D = oldNumberCell1Ds + numberCell1Ds;
    for (unsigned int nc = 0; nc < numberCell1Ds; nc++)
    {
      for (unsigned int v = 0; v < 2; v++)
        _mesh.Cell1DVertices.push_back(0.0);
      _mesh.Cell1DMarkers.push_back(0);
      _mesh.ActiveCell1D.push_back(false);
      _mesh.Cell1DOriginalCell1Ds.push_back(std::numeric_limits<unsigned int>::max());
      _mesh.NumberCell1DNeighbourCell2D.push_back(0);
      for (unsigned int p = 0; p < Cell1DNumberDoubleProperties(); p++)
        _mesh.Cell1DDoublePropertySizes[p].push_back(0);
    }

    for (unsigned c = 0; c < numberCell1Ds; c++)
      _mesh.NumberCell1DNeighbourCell2D[oldNumberCell1Ds + c + 1] = _mesh.NumberCell1DNeighbourCell2D[oldNumberCell1Ds];

    for (unsigned int p = 0; p < Cell1DNumberDoubleProperties(); p++)
    {
      for (unsigned c = 0; c < numberCell1Ds; c++)
        _mesh.Cell1DDoublePropertySizes[p][oldNumberCell1Ds + c + 1] = _mesh.Cell1DDoublePropertySizes[p][oldNumberCell1Ds];
    }

    return oldNumberCell1Ds;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DRemove(const unsigned int& cell1DIndex)
  {
    Output::Assert(cell1DIndex < Cell1DTotalNumber());

    for (unsigned int p = 0; p < Cell1DNumberDoubleProperties(); p++)
    {
      ResizeNumberVectorWithNewNumberElements(_mesh.Cell1DDoublePropertySizes[p],
                                              _mesh.Cell1DDoublePropertyValues[p],
                                              _mesh.NumberCell1D,
                                              cell1DIndex,
                                              0);
      _mesh.Cell1DDoublePropertySizes[p].erase(std::next(_mesh.Cell1DDoublePropertySizes[p].begin(),
                                                         cell1DIndex));
    }

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell1DNeighbourCell2D,
                                            _mesh.Cell1DNeighbourCell2Ds,
                                            _mesh.NumberCell1D,
                                            cell1DIndex,
                                            0);
    _mesh.NumberCell1DNeighbourCell2D.erase(std::next(_mesh.NumberCell1DNeighbourCell2D.begin(),
                                                      cell1DIndex));

    _mesh.UpdatedCell1Ds.erase(cell1DIndex);

    _mesh.Cell1DVertices.erase(std::next(_mesh.Cell1DVertices.begin(), 2 * cell1DIndex),
                               std::next(_mesh.Cell1DVertices.begin(), 2 * cell1DIndex + 2));
    _mesh.Cell1DOriginalCell1Ds.erase(std::next(_mesh.Cell1DOriginalCell1Ds.begin(), cell1DIndex));
    _mesh.Cell1DMarkers.erase(std::next(_mesh.Cell1DMarkers.begin(), cell1DIndex));
    _mesh.ActiveCell1D.erase(std::next(_mesh.ActiveCell1D.begin(), cell1DIndex));
    _mesh.NumberCell1D--;

    AlignContainerHigherElements(_mesh.Cell0DNeighbourCell1Ds,
                                 cell1DIndex,
                                 std::numeric_limits<unsigned int>::max());
    AlignContainerHigherElements(_mesh.Cell2DEdges,
                                 cell1DIndex,
                                 _mesh.NumberCell1D);
    AlignContainerHigherElements(_mesh.Cell3DEdges,
                                 cell1DIndex,
                                 _mesh.NumberCell1D);

    AlignMapContainerHigherElements(_mesh.UpdatedCell1Ds,
                                    cell1DIndex,
                                    _mesh.NumberCell1D);
    AlignContainerHigherElements(_mesh.Cell1DOriginalCell1Ds,
                                 cell1DIndex,
                                 std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DInsertExtremes(const unsigned int& cell1DIndex,
                                             const unsigned int& originCell0DIndex,
                                             const unsigned int& endCell0DIndex)
  {
    Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
    _mesh.Cell1DVertices[2 * cell1DIndex] = originCell0DIndex;
    _mesh.Cell1DVertices[2 * cell1DIndex + 1] = endCell0DIndex;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DsInsertExtremes(const Eigen::MatrixXi& cell1DExtremes)
  {
    Gedim::Output::Assert(cell1DExtremes.rows() == 2 &&
                          cell1DExtremes.cols() == Cell1DTotalNumber());

    std::list<Eigen::Triplet<unsigned int>> triplets;
    for (unsigned int cell1DIndex = 0; cell1DIndex < Cell1DTotalNumber(); cell1DIndex++)
    {
      const unsigned int& originCell0DIndex = cell1DExtremes(0, cell1DIndex);
      const unsigned int& endCell0DIndex = cell1DExtremes(1, cell1DIndex);

      Gedim::Output::Assert(originCell0DIndex < Cell0DTotalNumber());
      Gedim::Output::Assert(endCell0DIndex < Cell0DTotalNumber());
      _mesh.Cell1DVertices[2 * cell1DIndex] = originCell0DIndex;
      _mesh.Cell1DVertices[2 * cell1DIndex + 1] = endCell0DIndex;
      triplets.push_back(Eigen::Triplet<unsigned int>(originCell0DIndex,
                                                      endCell0DIndex,
                                                      cell1DIndex + 1));
    }
  }
  // ***************************************************************************
  Eigen::MatrixXi MeshMatricesDAO::Cell1DsExtremes() const
  {
    Eigen::MatrixXi extremes(2, _mesh.NumberCell1D);

    for (unsigned int e = 0; e < _mesh.NumberCell1D; e++)
      extremes.col(e)<< Cell1DOrigin(e), Cell1DEnd(e);

    return extremes;
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell1DByExtremes(const unsigned int& originCell0DIndex,
                                                 const unsigned int& endCell0DIndex) const
  {
    Gedim::Output::Assert(originCell0DIndex < Cell0DTotalNumber());
    Gedim::Output::Assert(endCell0DIndex < Cell0DTotalNumber());

    for (unsigned int e = 0; e < Cell1DTotalNumber(); e++)
    {
      if (_mesh.Cell1DVertices[2 * e] == originCell0DIndex &&
          _mesh.Cell1DVertices[2 * e + 1] == endCell0DIndex)
        return e;
    }

    return Cell1DTotalNumber();
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DsInitializeNeighbourCell2Ds(const unsigned int& numberNeighbourCell2Ds)
  {
    InitializeNuberVectorWithConstantElements(_mesh.NumberCell1DNeighbourCell2D,
                                              _mesh.Cell1DNeighbourCell2Ds,
                                              _mesh.NumberCell1D,
                                              numberNeighbourCell2Ds,
                                              std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DsInitializeNeighbourCell2Ds(const std::vector<unsigned int>& numberNeighbourCell2Ds)
  {
    Output::Assert(numberNeighbourCell2Ds.size() == Cell1DTotalNumber());
    InitializeNumberVector(_mesh.NumberCell1DNeighbourCell2D,
                           _mesh.Cell1DNeighbourCell2Ds,
                           numberNeighbourCell2Ds,
                           std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DInitializeNeighbourCell2Ds(const unsigned int& cell1DIndex,
                                                         const unsigned int& numberNeighbourCell2Ds)
  {
    Output::Assert(cell1DIndex < Cell1DTotalNumber());

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell1DNeighbourCell2D,
                                            _mesh.Cell1DNeighbourCell2Ds,
                                            _mesh.NumberCell1D,
                                            cell1DIndex,
                                            numberNeighbourCell2Ds,
                                            std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DInsertUpdatedCell1D(const unsigned int& cell1DIndex,
                                                  const unsigned int& updatedCell1DIdex)
  {
    Output::Assert(cell1DIndex < Cell1DTotalNumber());
    Output::Assert(updatedCell1DIdex < Cell1DTotalNumber());

    if (!Cell1DHasUpdatedCell1Ds(cell1DIndex))
      _mesh.UpdatedCell1Ds.insert(pair<unsigned int, unordered_set<unsigned int>>(cell1DIndex, {}));
    _mesh.UpdatedCell1Ds.at(cell1DIndex).insert(updatedCell1DIdex);
    _mesh.Cell1DOriginalCell1Ds[updatedCell1DIdex] = cell1DIndex;
  }
  // ***************************************************************************
  bool MeshMatricesDAO::Cell1DUpdatedCell1Ds(const unsigned int& cell1DIndex,
                                             list<unsigned int>& updatedCell1DIds) const
  {
    Output::Assert(cell1DIndex < Cell1DTotalNumber());
    unordered_map<unsigned int, unordered_set<unsigned int>>::const_iterator iter = _mesh.UpdatedCell1Ds.find(cell1DIndex);

    if (iter == _mesh.UpdatedCell1Ds.end())
    {
      updatedCell1DIds.push_back(cell1DIndex);
      return true;
    }

    for (const unsigned int& childId : iter->second)
    {
      Cell1DUpdatedCell1Ds(childId,
                           updatedCell1DIds);
    }

    return false;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DInitializeDoubleProperties(const unsigned int& numberDoubleProperties)
  {
    _mesh.Cell1DDoublePropertyIds.reserve(numberDoubleProperties);
    _mesh.Cell1DDoublePropertySizes.reserve(numberDoubleProperties);
    _mesh.Cell1DDoublePropertyValues.reserve(numberDoubleProperties);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell1DAddDoubleProperty(const string& propertyId)
  {
    Output::Assert(!Cell1DDoublePropertyExists(propertyId));
    unsigned int propertyIndex = _mesh.Cell1DDoublePropertySizes.size();

    _mesh.Cell1DDoublePropertyIndices.insert(pair<string, unsigned int>(propertyId,
                                                                        propertyIndex));
    _mesh.Cell1DDoublePropertyIds.push_back(propertyId);
    _mesh.Cell1DDoublePropertySizes.push_back(vector<unsigned int>(_mesh.NumberCell1D + 1, 0));
    _mesh.Cell1DDoublePropertyValues.push_back(vector<double>());
    return propertyIndex;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DsInitializeDoublePropertyValues(const unsigned int& propertyIndex,
                                                              const std::vector<unsigned int>& propertySizes)
  {
    Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());
    Output::Assert(propertySizes.size() == Cell1DTotalNumber());

    InitializeNumberVector(_mesh.Cell1DDoublePropertySizes[propertyIndex],
                           _mesh.Cell1DDoublePropertyValues[propertyIndex],
                           propertySizes,
                           0.0);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DInitializeDoublePropertyValues(const unsigned int& cell1DIndex,
                                                             const unsigned int& propertyIndex,
                                                             const unsigned int& propertySize)
  {
    Output::Assert(cell1DIndex < Cell1DTotalNumber());
    Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());

    ResizeNumberVectorWithNewNumberElements(_mesh.Cell1DDoublePropertySizes[propertyIndex],
                                            _mesh.Cell1DDoublePropertyValues[propertyIndex],
                                            _mesh.NumberCell1D,
                                            cell1DIndex,
                                            propertySize,
                                            0.0);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DsInitialize(const unsigned int& numberCell2Ds)
  {
    _mesh.NumberCell2D = numberCell2Ds;
    _mesh.NumberCell2DVertices.resize(_mesh.NumberCell2D + 1, 0);
    _mesh.NumberCell2DEdges.resize(_mesh.NumberCell2D + 1, 0);
    _mesh.NumberCell2DSubdivision.resize(_mesh.NumberCell2D + 1, 0);
    _mesh.Cell2DMarkers.resize(_mesh.NumberCell2D, 0);
    _mesh.ActiveCell2D.resize(_mesh.NumberCell2D, false);
    _mesh.Cell2DOriginalCell2Ds.resize(_mesh.NumberCell2D, std::numeric_limits<unsigned int>::max());
    _mesh.NumberCell2DNeighbourCell3D.resize(_mesh.NumberCell2D + 1, 0);
    for (unsigned int p = 0; p < Cell2DNumberDoubleProperties(); p++)
      _mesh.Cell2DDoublePropertySizes[p].resize(_mesh.NumberCell2D + 1, 0);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell2DAppend(const unsigned int& numberCell2Ds)
  {
    unsigned int oldNumberCell2Ds = _mesh.NumberCell2D;

    _mesh.NumberCell2D = oldNumberCell2Ds + numberCell2Ds;
    for (unsigned int c = 0; c < numberCell2Ds; c++)
    {
      _mesh.NumberCell2DVertices.push_back(0);
      _mesh.NumberCell2DEdges.push_back(0);
      _mesh.NumberCell2DSubdivision.push_back(0);
      _mesh.Cell2DMarkers.push_back(0);
      _mesh.ActiveCell2D.push_back(false);
      _mesh.Cell2DOriginalCell2Ds.push_back(std::numeric_limits<unsigned int>::max());
      _mesh.NumberCell2DNeighbourCell3D.push_back(0);
      for (unsigned int p = 0; p < Cell2DNumberDoubleProperties(); p++)
        _mesh.Cell2DDoublePropertySizes[p].push_back(0);
    }

    for (unsigned c = 0; c < numberCell2Ds; c++)
    {
      _mesh.NumberCell2DVertices[oldNumberCell2Ds + c + 1] = _mesh.NumberCell2DVertices[oldNumberCell2Ds];
      _mesh.NumberCell2DEdges[oldNumberCell2Ds + c + 1] = _mesh.NumberCell2DEdges[oldNumberCell2Ds];
      _mesh.NumberCell2DSubdivision[oldNumberCell2Ds + c + 1] = _mesh.NumberCell2DSubdivision[oldNumberCell2Ds];
    }

    for (unsigned c = 0; c < numberCell2Ds; c++)
      _mesh.NumberCell2DNeighbourCell3D[oldNumberCell2Ds + c + 1] = _mesh.NumberCell2DNeighbourCell3D[oldNumberCell2Ds];

    for (unsigned int p = 0; p < Cell2DNumberDoubleProperties(); p++)
    {
      for (unsigned c = 0; c < numberCell2Ds; c++)
        _mesh.Cell2DDoublePropertySizes[p][oldNumberCell2Ds + c + 1] = _mesh.Cell2DDoublePropertySizes[p][oldNumberCell2Ds];
    }

    return oldNumberCell2Ds;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DRemove(const unsigned int& cell2DIndex)
  {
    Output::Assert(cell2DIndex < Cell2DTotalNumber());

    for (unsigned int p = 0; p < Cell2DNumberDoubleProperties(); p++)
    {
      ResizeNumberVectorWithNewNumberElements(_mesh.Cell2DDoublePropertySizes[p],
                                              _mesh.Cell2DDoublePropertyValues[p],
                                              _mesh.NumberCell2D,
                                              cell2DIndex,
                                              0);
      _mesh.Cell2DDoublePropertySizes[p].erase(std::next(_mesh.Cell2DDoublePropertySizes[p].begin(),
                                                         cell2DIndex));
    }

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell2DNeighbourCell3D,
                                            _mesh.Cell2DNeighbourCell3Ds,
                                            _mesh.NumberCell2D,
                                            cell2DIndex,
                                            0);
    _mesh.NumberCell2DNeighbourCell3D.erase(std::next(_mesh.NumberCell2DNeighbourCell3D.begin(),
                                                      cell2DIndex));

    _mesh.UpdatedCell2Ds.erase(cell2DIndex);

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell2DVertices,
                                            _mesh.Cell2DVertices,
                                            _mesh.NumberCell2D,
                                            cell2DIndex,
                                            0);
    _mesh.NumberCell2DVertices.erase(std::next(_mesh.NumberCell2DVertices.begin(),
                                               cell2DIndex));

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell2DEdges,
                                            _mesh.Cell2DEdges,
                                            _mesh.NumberCell2D,
                                            cell2DIndex,
                                            0);
    _mesh.NumberCell2DEdges.erase(std::next(_mesh.NumberCell2DEdges.begin(),
                                            cell2DIndex));

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell2DSubdivision,
                                            _mesh.Cell2DSubdivision,
                                            _mesh.NumberCell2D,
                                            cell2DIndex,
                                            0);
    _mesh.NumberCell2DSubdivision.erase(std::next(_mesh.NumberCell2DSubdivision.begin(),
                                                  cell2DIndex));

    _mesh.Cell2DOriginalCell2Ds.erase(std::next(_mesh.Cell2DOriginalCell2Ds.begin(),
                                                cell2DIndex));
    _mesh.Cell2DMarkers.erase(std::next(_mesh.Cell2DMarkers.begin(), cell2DIndex));
    _mesh.ActiveCell2D.erase(std::next(_mesh.ActiveCell2D.begin(), cell2DIndex));
    _mesh.NumberCell2D--;

    AlignContainerHigherElements(_mesh.Cell0DNeighbourCell2Ds,
                                 cell2DIndex,
                                 std::numeric_limits<unsigned int>::max());
    AlignContainerHigherElements(_mesh.Cell1DNeighbourCell2Ds,
                                 cell2DIndex,
                                 std::numeric_limits<unsigned int>::max());
    AlignContainerHigherElements(_mesh.Cell3DFaces,
                                 cell2DIndex,
                                 _mesh.NumberCell2D);
    AlignMapContainerHigherElements(_mesh.UpdatedCell2Ds,
                                    cell2DIndex,
                                    _mesh.NumberCell2D);
    AlignContainerHigherElements(_mesh.Cell2DOriginalCell2Ds,
                                 cell2DIndex,
                                 std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DsInitializeVertices(const unsigned int& numberCell2DVertices)
  {
    InitializeNuberVectorWithConstantElements(_mesh.NumberCell2DVertices,
                                              _mesh.Cell2DVertices,
                                              _mesh.NumberCell2D,
                                              numberCell2DVertices,
                                              std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DsInitializeVertices(const std::vector<unsigned int>& numberCell2DsVertices)
  {
    Gedim::Output::Assert(numberCell2DsVertices.size() ==
                          _mesh.NumberCell2D);
    InitializeNumberVector(_mesh.NumberCell2DVertices,
                           _mesh.Cell2DVertices,
                           numberCell2DsVertices,
                           std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DInitializeVertices(const unsigned int& cell2DIndex,
                                                 const unsigned int& numberCell2DVertices)
  {
    Output::Assert(cell2DIndex < Cell2DTotalNumber());
    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell2DVertices,
                                            _mesh.Cell2DVertices,
                                            _mesh.NumberCell2D,
                                            cell2DIndex,
                                            numberCell2DVertices);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DsInitializeEdges(const unsigned int& numberCell2DEdges)
  {
    InitializeNuberVectorWithConstantElements(_mesh.NumberCell2DEdges,
                                              _mesh.Cell2DEdges,
                                              _mesh.NumberCell2D,
                                              numberCell2DEdges,
                                              std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DsInitializeEdges(const std::vector<unsigned int>& numberCell2DsEdges)
  {
    Output::Assert(numberCell2DsEdges.size() ==
                   Cell2DTotalNumber());
    InitializeNumberVector(_mesh.NumberCell2DEdges,
                           _mesh.Cell2DEdges,
                           numberCell2DsEdges,
                           std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DInitializeEdges(const unsigned int& cell2DIndex,
                                              const unsigned int& numberCell2DEdges)
  {
    Output::Assert(cell2DIndex < Cell2DTotalNumber());
    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell2DEdges,
                                            _mesh.Cell2DEdges,
                                            _mesh.NumberCell2D,
                                            cell2DIndex,
                                            numberCell2DEdges);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DAddVertices(const unsigned int& cell2DIndex,
                                          const vector<unsigned int>& verticesCell0DIndices)
  {
    Cell2DInitializeVertices(cell2DIndex,
                             verticesCell0DIndices.size());
    for (unsigned int v = 0; v < verticesCell0DIndices.size(); v++)
      Cell2DInsertVertex(cell2DIndex,
                         v,
                         verticesCell0DIndices[v]);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DAddEdges(const unsigned int& cell2DIndex,
                                       const vector<unsigned int>& edgesCell1DIndices)
  {
    Cell2DInitializeEdges(cell2DIndex,
                          edgesCell1DIndices.size());
    for (unsigned int e = 0; e < edgesCell1DIndices.size(); e++)
      Cell2DInsertEdge(cell2DIndex,
                       e,
                       edgesCell1DIndices[e]);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DAddVerticesAndEdges(const unsigned int& cell2DIndex,
                                                  const Eigen::MatrixXi& verticesAndEdgesIndices)
  {
    Cell2DInitializeVertices(cell2DIndex,
                             verticesAndEdgesIndices.cols());
    Cell2DInitializeEdges(cell2DIndex,
                          verticesAndEdgesIndices.cols());

    for (unsigned int v = 0; v < verticesAndEdgesIndices.cols(); v++)
    {
      Cell2DInsertVertex(cell2DIndex,
                         v,
                         verticesAndEdgesIndices(0, v));
      Cell2DInsertEdge(cell2DIndex,
                       v,
                       verticesAndEdgesIndices(1, v));
    }
  }
  // ***************************************************************************
  std::vector<std::vector<unsigned int>> MeshMatricesDAO::Cell2DsVertices() const
  {
    vector<std::vector<unsigned int>> polygonVertices(Cell2DTotalNumber());

    for (unsigned int p = 0; p < Cell2DTotalNumber(); p++)
      polygonVertices[p] = Cell2DVertices(p);

    return polygonVertices;
  }
  // ***************************************************************************
  MatrixXd MeshMatricesDAO::Cell2DVerticesCoordinates(const unsigned int& cell2DIndex) const
  {
    MatrixXd coordinates(3, Cell2DNumberVertices(cell2DIndex));
    for (unsigned int v = 0; v < Cell2DNumberVertices(cell2DIndex); v++)
      coordinates.col(v) << Cell2DVertexCoordinates(cell2DIndex, v);
    return coordinates;
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell2DFindVertex(const unsigned int& cell2DIndex,
                                                 const unsigned int& cell0DIndex) const
  {
    const vector<unsigned int> vertices = Cell2DVertices(cell2DIndex);
    const vector<unsigned int>::const_iterator it = std::find(vertices.begin(),
                                                              vertices.end(),
                                                              cell0DIndex);
    if (it != vertices.end())
      return std::distance(vertices.begin(), it);
    else
      throw runtime_error("Vertex not found");
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell2DFindEdge(const unsigned int& cell2DIndex,
                                               const unsigned int& cell1DIndex) const
  {
    const vector<unsigned int> edges = Cell2DEdges(cell2DIndex);
    const vector<unsigned int>::const_iterator it = std::find(edges.begin(),
                                                              edges.end(),
                                                              cell1DIndex);
    if (it != edges.end())
      return std::distance(edges.begin(), it);
    else
      throw runtime_error("Edge not found");
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell2DFindEdgeByExtremes(const unsigned int& cell2DIndex,
                                                         const unsigned int& originCell0DIndex,
                                                         const unsigned int& endCell0DIndex) const
  {
    const vector<unsigned int> edges = Cell2DEdges(cell2DIndex);

    for (unsigned int e = 0; e < edges.size(); e++)
    {
      if (_mesh.Cell1DVertices[2 * edges[e]] == originCell0DIndex &&
          _mesh.Cell1DVertices[2 * edges[e] + 1] == endCell0DIndex)
        return e;
    }

    return edges.size();
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DInsertUpdatedCell2D(const unsigned int& cell2DIndex,
                                                  const unsigned int& updatedCell2DIdex)
  {
    Output::Assert(cell2DIndex < Cell2DTotalNumber());
    Output::Assert(updatedCell2DIdex < Cell2DTotalNumber());

    if (!Cell2DHasUpdatedCell2Ds(cell2DIndex))
      _mesh.UpdatedCell2Ds.insert(pair<unsigned int, unordered_set<unsigned int>>(cell2DIndex, {}));
    _mesh.UpdatedCell2Ds.at(cell2DIndex).insert(updatedCell2DIdex);
    _mesh.Cell2DOriginalCell2Ds[updatedCell2DIdex] = cell2DIndex;
  }
  // ***************************************************************************
  bool MeshMatricesDAO::Cell2DUpdatedCell2Ds(const unsigned int& cell2DIndex,
                                             list<unsigned int>& updatedCell2DIds) const
  {
    Output::Assert(cell2DIndex < Cell2DTotalNumber());
    unordered_map<unsigned int, unordered_set<unsigned int>>::const_iterator iter = _mesh.UpdatedCell2Ds.find(cell2DIndex);

    if (iter == _mesh.UpdatedCell2Ds.end())
    {
      updatedCell2DIds.push_back(cell2DIndex);
      return true;
    }

    for (const unsigned int& childId : iter->second)
    {
      Cell2DUpdatedCell2Ds(childId,
                           updatedCell2DIds);
    }

    return false;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DsInitializeNeighbourCell3Ds(const std::vector<unsigned int>& numberNeighbourCell3Ds)
  {
    Output::Assert(numberNeighbourCell3Ds.size() == Cell2DTotalNumber());
    InitializeNumberVector(_mesh.NumberCell2DNeighbourCell3D,
                           _mesh.Cell2DNeighbourCell3Ds,
                           numberNeighbourCell3Ds,
                           std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DInitializeNeighbourCell3Ds(const unsigned int& cell2DIndex,
                                                         const unsigned int& numberNeighbourCell3Ds)
  {
    Output::Assert(cell2DIndex < Cell2DTotalNumber());

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell2DNeighbourCell3D,
                                            _mesh.Cell2DNeighbourCell3Ds,
                                            _mesh.NumberCell2D,
                                            cell2DIndex,
                                            numberNeighbourCell3Ds,
                                            std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DInitializeDoubleProperties(const unsigned int& numberDoubleProperties)
  {
    _mesh.Cell2DDoublePropertyIds.reserve(numberDoubleProperties);
    _mesh.Cell2DDoublePropertySizes.reserve(numberDoubleProperties);
    _mesh.Cell2DDoublePropertyValues.reserve(numberDoubleProperties);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell2DAddDoubleProperty(const string& propertyId)
  {
    Output::Assert(!Cell2DDoublePropertyExists(propertyId));
    unsigned int propertyIndex = _mesh.Cell2DDoublePropertySizes.size();

    _mesh.Cell2DDoublePropertyIndices.insert(pair<string, unsigned int>(propertyId,
                                                                        propertyIndex));
    _mesh.Cell2DDoublePropertyIds.push_back(propertyId);
    _mesh.Cell2DDoublePropertySizes.push_back(vector<unsigned int>(_mesh.NumberCell2D + 1, 0));
    _mesh.Cell2DDoublePropertyValues.push_back(vector<double>());
    return propertyIndex;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DsInitializeDoublePropertyValues(const unsigned int& propertyIndex,
                                                              const std::vector<unsigned int>& propertySizes)
  {
    Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());
    Output::Assert(propertySizes.size() == Cell2DTotalNumber());

    InitializeNumberVector(_mesh.Cell2DDoublePropertySizes[propertyIndex],
                           _mesh.Cell2DDoublePropertyValues[propertyIndex],
                           propertySizes,
                           0.0);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DInitializeDoublePropertyValues(const unsigned int& cell2DIndex,
                                                             const unsigned int& propertyIndex,
                                                             const unsigned int& propertySize)
  {
    Output::Assert(cell2DIndex < Cell2DTotalNumber());
    Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());

    ResizeNumberVectorWithNewNumberElements(_mesh.Cell2DDoublePropertySizes[propertyIndex],
                                            _mesh.Cell2DDoublePropertyValues[propertyIndex],
                                            _mesh.NumberCell2D,
                                            cell2DIndex,
                                            propertySize,
                                            0.0);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DsInitializeSubDivision(const std::vector<unsigned int>& numberSubDivisions)
  {
    Output::Assert(numberSubDivisions.size() == Cell2DTotalNumber());

    InitializeNumberVector(_mesh.NumberCell2DSubdivision,
                           _mesh.Cell2DSubdivision,
                           numberSubDivisions,
                           static_cast<unsigned int>(0));
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DInitializeSubDivision(const unsigned int& cell2DIndex,
                                                    const unsigned int& numberSubDivision)
  {
    Output::Assert(cell2DIndex < Cell2DTotalNumber());
    Output::Assert(numberSubDivision % 3 == 0);
    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell2DSubdivision,
                                            _mesh.Cell2DSubdivision,
                                            _mesh.NumberCell2D,
                                            cell2DIndex,
                                            numberSubDivision);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DsInitialize(const unsigned int& numberCell3Ds)
  {
    _mesh.NumberCell3D = numberCell3Ds;
    _mesh.NumberCell3DVertices.resize(_mesh.NumberCell3D + 1, 0);
    _mesh.NumberCell3DEdges.resize(_mesh.NumberCell3D + 1, 0);
    _mesh.NumberCell3DFaces.resize(_mesh.NumberCell3D + 1, 0);
    _mesh.Cell3DMarkers.resize(_mesh.NumberCell3D, 0);
    _mesh.ActiveCell3D.resize(_mesh.NumberCell3D, false);
    for (unsigned int p = 0; p < Cell3DNumberDoubleProperties(); p++)
      _mesh.Cell3DDoublePropertySizes[p].resize(_mesh.NumberCell3D + 1, 0);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell3DAppend(const unsigned int& numberCell3Ds)
  {
    unsigned int oldNumberCell3Ds = _mesh.NumberCell3D;

    _mesh.NumberCell3D = oldNumberCell3Ds + numberCell3Ds;

    for (unsigned c = 0; c < numberCell3Ds; c++)
    {
      _mesh.NumberCell3DVertices.push_back(0);
      _mesh.NumberCell3DEdges.push_back(0);
      _mesh.NumberCell3DFaces.push_back(0);
      _mesh.Cell3DMarkers.push_back(0);
      _mesh.ActiveCell3D.push_back(false);
      for (unsigned int p = 0; p < Cell3DNumberDoubleProperties(); p++)
        _mesh.Cell3DDoublePropertySizes[p].push_back(0);
    }

    for (unsigned c = 0; c < numberCell3Ds; c++)
    {
      _mesh.NumberCell3DVertices[oldNumberCell3Ds + c + 1] = _mesh.NumberCell3DVertices[oldNumberCell3Ds];
      _mesh.NumberCell3DEdges[oldNumberCell3Ds + c + 1] = _mesh.NumberCell3DEdges[oldNumberCell3Ds];
      _mesh.NumberCell3DFaces[oldNumberCell3Ds + c + 1] = _mesh.NumberCell3DFaces[oldNumberCell3Ds];
    }

    for (unsigned int p = 0; p < Cell3DNumberDoubleProperties(); p++)
    {
      for (unsigned c = 0; c < numberCell3Ds; c++)
        _mesh.Cell3DDoublePropertySizes[p][oldNumberCell3Ds + c + 1] = _mesh.Cell3DDoublePropertySizes[p][oldNumberCell3Ds];
    }

    return oldNumberCell3Ds;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DRemove(const unsigned int& cell3DIndex)
  {
    Output::Assert(cell3DIndex < Cell3DTotalNumber());

    for (unsigned int p = 0; p < Cell3DNumberDoubleProperties(); p++)
    {
      ResizeNumberVectorWithNewNumberElements(_mesh.Cell3DDoublePropertySizes[p],
                                              _mesh.Cell3DDoublePropertyValues[p],
                                              _mesh.NumberCell3D,
                                              cell3DIndex,
                                              0);
      _mesh.Cell3DDoublePropertySizes[p].erase(std::next(_mesh.Cell3DDoublePropertySizes[p].begin(),
                                                         cell3DIndex));
    }

    _mesh.UpdatedCell3Ds.erase(cell3DIndex);

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell3DVertices,
                                            _mesh.Cell3DVertices,
                                            _mesh.NumberCell3D,
                                            cell3DIndex,
                                            0);
    _mesh.NumberCell3DVertices.erase(std::next(_mesh.NumberCell3DVertices.begin(),
                                               cell3DIndex));

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell3DEdges,
                                            _mesh.Cell3DEdges,
                                            _mesh.NumberCell3D,
                                            cell3DIndex,
                                            0);
    _mesh.NumberCell3DEdges.erase(std::next(_mesh.NumberCell3DEdges.begin(),
                                            cell3DIndex));

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell3DFaces,
                                            _mesh.Cell3DFaces,
                                            _mesh.NumberCell3D,
                                            cell3DIndex,
                                            0);
    _mesh.NumberCell3DFaces.erase(std::next(_mesh.NumberCell3DFaces.begin(),
                                            cell3DIndex));

    _mesh.Cell3DMarkers.erase(std::next(_mesh.Cell3DMarkers.begin(), cell3DIndex));
    _mesh.ActiveCell3D.erase(std::next(_mesh.ActiveCell3D.begin(), cell3DIndex));
    _mesh.NumberCell3D--;

    AlignContainerHigherElements(_mesh.Cell2DNeighbourCell3Ds,
                                 cell3DIndex,
                                 std::numeric_limits<unsigned int>::max());

    AlignMapContainerHigherElements(_mesh.UpdatedCell3Ds,
                                    cell3DIndex,
                                    _mesh.NumberCell3D);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DsInitializeVertices(const std::vector<unsigned int>& numberCell3DsVertices)
  {
    Gedim::Output::Assert(numberCell3DsVertices.size() ==
                          _mesh.NumberCell3D);
    InitializeNumberVector(_mesh.NumberCell3DVertices,
                           _mesh.Cell3DVertices,
                           numberCell3DsVertices,
                           std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DsInitializeEdges(const std::vector<unsigned int>& numberCell3DsEdges)
  {
    Gedim::Output::Assert(numberCell3DsEdges.size() ==
                          _mesh.NumberCell3D);
    InitializeNumberVector(_mesh.NumberCell3DEdges,
                           _mesh.Cell3DEdges,
                           numberCell3DsEdges,
                           std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DsInitializeFaces(const std::vector<unsigned int>& numberCell3DsFaces)
  {
    Gedim::Output::Assert(numberCell3DsFaces.size() ==
                          _mesh.NumberCell3D);
    InitializeNumberVector(_mesh.NumberCell3DFaces,
                           _mesh.Cell3DFaces,
                           numberCell3DsFaces,
                           std::numeric_limits<unsigned int>::max());
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DInitializeVertices(const unsigned int& cell3DIndex,
                                                 const unsigned int& numberCell3DVertices)
  {
    Output::Assert(cell3DIndex < Cell3DTotalNumber());
    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell3DVertices,
                                            _mesh.Cell3DVertices,
                                            _mesh.NumberCell3D,
                                            cell3DIndex,
                                            numberCell3DVertices);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DInitializeEdges(const unsigned int& cell3DIndex,
                                              const unsigned int& numberCell3DEdges)
  {
    Output::Assert(cell3DIndex < Cell3DTotalNumber());
    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell3DEdges,
                                            _mesh.Cell3DEdges,
                                            _mesh.NumberCell3D,
                                            cell3DIndex,
                                            numberCell3DEdges);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DInitializeFaces(const unsigned int& cell3DIndex,
                                              const unsigned int& numberCell3DFaces)
  {
    Output::Assert(cell3DIndex < Cell3DTotalNumber());
    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell3DFaces,
                                            _mesh.Cell3DFaces,
                                            _mesh.NumberCell3D,
                                            cell3DIndex,
                                            numberCell3DFaces);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DAddVertices(const unsigned int& cell3DIndex,
                                          const vector<unsigned int>& verticesCell0DIndices)
  {
    Cell3DInitializeVertices(cell3DIndex,
                             verticesCell0DIndices.size());
    for (unsigned int v = 0; v < verticesCell0DIndices.size(); v++)
      Cell3DInsertVertex(cell3DIndex,
                         v,
                         verticesCell0DIndices[v]);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DAddEdges(const unsigned int& cell3DIndex,
                                       const vector<unsigned int>& edgesCell0DIndices)
  {
    Cell3DInitializeEdges(cell3DIndex,
                          edgesCell0DIndices.size());
    for (unsigned int e = 0; e < edgesCell0DIndices.size(); e++)
      Cell3DInsertEdge(cell3DIndex,
                       e,
                       edgesCell0DIndices[e]);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DAddFaces(const unsigned int& cell3DIndex,
                                       const vector<unsigned int>& facesCell0DIndices)
  {
    Cell3DInitializeFaces(cell3DIndex,
                          facesCell0DIndices.size());
    for (unsigned int e = 0; e < facesCell0DIndices.size(); e++)
      Cell3DInsertFace(cell3DIndex,
                       e,
                       facesCell0DIndices[e]);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell3DFindEdgeByExtremes(const unsigned int& cell3DIndex,
                                                         const unsigned int& originCell0DIndex,
                                                         const unsigned int& endCell0DIndex) const
  {
    const vector<unsigned int> edges = Cell3DEdges(cell3DIndex);

    for (unsigned int e = 0; e < edges.size(); e++)
    {
      if (_mesh.Cell1DVertices[2 * edges[e]] == originCell0DIndex &&
          _mesh.Cell1DVertices[2 * edges[e] + 1] == endCell0DIndex)
        return e;
    }

    return edges.size();
  }
  // ***************************************************************************
  std::vector<unsigned int> MeshMatricesDAO::Cell3DVertices(const unsigned int& cell3DIndex) const
  {
    const unsigned int numVertices = Cell3DNumberVertices(cell3DIndex);

    vector<unsigned int> vertices(numVertices);
    for (unsigned int v = 0; v < numVertices; v++)
      vertices[v] = _mesh.Cell3DVertices[_mesh.NumberCell3DVertices[cell3DIndex] + v];

    return vertices;
  }
  // ***************************************************************************
  MatrixXd MeshMatricesDAO::Cell3DVerticesCoordinates(const unsigned int& cell3DIndex) const
  {
    MatrixXd coordinates(3, Cell3DNumberVertices(cell3DIndex));
    for (unsigned int v = 0; v < Cell3DNumberVertices(cell3DIndex); v++)
      coordinates.col(v) << Cell3DVertexCoordinates(cell3DIndex, v);
    return coordinates;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DInsertUpdatedCell3D(const unsigned int& cell3DIndex,
                                                  const unsigned int& updatedCell3DIdex)
  {
    Output::Assert(cell3DIndex < Cell3DTotalNumber());
    Output::Assert(updatedCell3DIdex < Cell3DTotalNumber());

    if (!Cell3DHasUpdatedCell3Ds(cell3DIndex))
      _mesh.UpdatedCell3Ds.insert(pair<unsigned int, unordered_set<unsigned int>>(cell3DIndex, {}));
    _mesh.UpdatedCell3Ds.at(cell3DIndex).insert(updatedCell3DIdex);
  }
  // ***************************************************************************
  bool MeshMatricesDAO::Cell3DUpdatedCell3Ds(const unsigned int& cell3DIndex,
                                             list<unsigned int>& updatedCell3DIds) const
  {
    Output::Assert(cell3DIndex < Cell3DTotalNumber());
    unordered_map<unsigned int, unordered_set<unsigned int>>::const_iterator iter = _mesh.UpdatedCell3Ds.find(cell3DIndex);

    if (iter == _mesh.UpdatedCell3Ds.end())
    {
      updatedCell3DIds.push_back(cell3DIndex);
      return true;
    }

    for (const unsigned int& childId : iter->second)
    {
      Cell3DUpdatedCell3Ds(childId,
                           updatedCell3DIds);
    }

    return false;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DInitializeDoubleProperties(const unsigned int& numberDoubleProperties)
  {
    _mesh.Cell3DDoublePropertyIds.reserve(numberDoubleProperties);
    _mesh.Cell3DDoublePropertySizes.reserve(numberDoubleProperties);
    _mesh.Cell3DDoublePropertyValues.reserve(numberDoubleProperties);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell3DAddDoubleProperty(const string& propertyId)
  {
    Output::Assert(!Cell3DDoublePropertyExists(propertyId));
    unsigned int propertyIndex = _mesh.Cell3DDoublePropertySizes.size();

    _mesh.Cell3DDoublePropertyIndices.insert(pair<string, unsigned int>(propertyId,
                                                                        propertyIndex));
    _mesh.Cell3DDoublePropertyIds.push_back(propertyId);
    _mesh.Cell3DDoublePropertySizes.push_back(vector<unsigned int>(_mesh.NumberCell3D + 1, 0));
    _mesh.Cell3DDoublePropertyValues.push_back(vector<double>());
    return propertyIndex;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DsInitializeDoublePropertyValues(const unsigned int& propertyIndex,
                                                              const std::vector<unsigned int>& propertySizes)
  {
    Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());
    Output::Assert(propertySizes.size() == Cell3DTotalNumber());

    InitializeNumberVector(_mesh.Cell3DDoublePropertySizes[propertyIndex],
                           _mesh.Cell3DDoublePropertyValues[propertyIndex],
                           propertySizes,
                           0.0);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell3DInitializeDoublePropertyValues(const unsigned int& cell3DIndex,
                                                             const unsigned int& propertyIndex,
                                                             const unsigned int& propertySize)
  {
    Output::Assert(cell3DIndex < Cell3DTotalNumber());
    Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());

    ResizeNumberVectorWithNewNumberElements(_mesh.Cell3DDoublePropertySizes[propertyIndex],
                                            _mesh.Cell3DDoublePropertyValues[propertyIndex],
                                            _mesh.NumberCell3D,
                                            cell3DIndex,
                                            propertySize,
                                            0.0);
  }
  // ***************************************************************************
  vector<unsigned int> MeshMatricesDAO::Cell3DEdges(const unsigned int& cell3DIndex) const
  {
    const unsigned int numEdges = Cell3DNumberEdges(cell3DIndex);

    vector<unsigned int> edges(numEdges);
    for (unsigned int e = 0; e < numEdges; e++)
      edges[e] = _mesh.Cell3DEdges[_mesh.NumberCell3DEdges[cell3DIndex] + e];

    return edges;
  }
  // ***************************************************************************
  vector<unsigned int> MeshMatricesDAO::Cell3DFaces(const unsigned int& cell3DIndex) const
  {
    const unsigned int numFaces = Cell3DNumberFaces(cell3DIndex);

    vector<unsigned int> faces(numFaces);
    for (unsigned int f = 0; f < numFaces; f++)
      faces[f] = _mesh.Cell3DFaces[_mesh.NumberCell3DFaces[cell3DIndex] + f];

    return faces;
  }
  // ***************************************************************************
  std::vector<std::vector<std::vector<unsigned int>>> MeshMatricesDAO::Cell3DsFacesVertices() const
  {
    vector<std::vector<std::vector<unsigned int>>> cell3DsFacesVertices(Cell3DTotalNumber());

    for (unsigned int p = 0; p < Cell3DTotalNumber(); p++)
    {
      const vector<unsigned int> cell3DFaces = Cell3DFaces(p);
      cell3DsFacesVertices[p].resize(cell3DFaces.size());

      for (unsigned int f = 0; f < cell3DFaces.size(); f++)
        cell3DsFacesVertices[p][f] = Cell2DVertices(cell3DFaces[f]);
    }

    return cell3DsFacesVertices;
  }
  // ***************************************************************************
  void MeshMatricesDAO::Compress()
  {
    _mesh.Cell0DCoordinates.shrink_to_fit();
    _mesh.Cell0DMarkers.shrink_to_fit();
    _mesh.ActiveCell0D.shrink_to_fit();
    _mesh.NumberCell0DNeighbourCell1D.shrink_to_fit();
    _mesh.Cell0DNeighbourCell1Ds.shrink_to_fit();
    _mesh.NumberCell0DNeighbourCell2D.shrink_to_fit();
    _mesh.Cell0DNeighbourCell2Ds.shrink_to_fit();
    _mesh.Cell0DDoublePropertyIds.shrink_to_fit();
    _mesh.Cell0DDoublePropertySizes.shrink_to_fit();
    for (unsigned int s = 0; s < _mesh.Cell0DDoublePropertySizes.size(); s++)
      _mesh.Cell0DDoublePropertySizes[s].shrink_to_fit();
    _mesh.Cell0DDoublePropertyValues.shrink_to_fit();
    for (unsigned int s = 0; s < _mesh.Cell0DDoublePropertyValues.size(); s++)
      _mesh.Cell0DDoublePropertyValues[s].shrink_to_fit();

    _mesh.Cell1DVertices.shrink_to_fit();
    _mesh.NumberCell1DNeighbourCell2D.shrink_to_fit();
    _mesh.Cell1DNeighbourCell2Ds.shrink_to_fit();
    _mesh.Cell1DMarkers.shrink_to_fit();
    _mesh.Cell1DOriginalCell1Ds.shrink_to_fit();
    _mesh.ActiveCell1D.shrink_to_fit();
    _mesh.Cell1DDoublePropertyIds.shrink_to_fit();
    _mesh.Cell1DDoublePropertySizes.shrink_to_fit();
    for (unsigned int s = 0; s < _mesh.Cell1DDoublePropertySizes.size(); s++)
      _mesh.Cell1DDoublePropertySizes[s].shrink_to_fit();
    _mesh.Cell1DDoublePropertyValues.shrink_to_fit();
    for (unsigned int s = 0; s < _mesh.Cell1DDoublePropertyValues.size(); s++)
      _mesh.Cell1DDoublePropertyValues[s].shrink_to_fit();

    _mesh.NumberCell2DVertices.shrink_to_fit();
    _mesh.NumberCell2DEdges.shrink_to_fit();
    _mesh.NumberCell2DSubdivision.shrink_to_fit();
    _mesh.Cell2DVertices.shrink_to_fit();
    _mesh.Cell2DEdges.shrink_to_fit();
    _mesh.Cell2DSubdivision.shrink_to_fit();
    _mesh.NumberCell2DNeighbourCell3D.shrink_to_fit();
    _mesh.Cell2DNeighbourCell3Ds.shrink_to_fit();
    _mesh.Cell2DMarkers.shrink_to_fit();
    _mesh.Cell2DOriginalCell2Ds.shrink_to_fit();
    _mesh.ActiveCell2D.shrink_to_fit();
    _mesh.Cell2DDoublePropertyIds.shrink_to_fit();
    _mesh.Cell2DDoublePropertySizes.shrink_to_fit();
    for (unsigned int s = 0; s < _mesh.Cell2DDoublePropertySizes.size(); s++)
      _mesh.Cell2DDoublePropertySizes[s].shrink_to_fit();
    _mesh.Cell2DDoublePropertyValues.shrink_to_fit();
    for (unsigned int s = 0; s < _mesh.Cell2DDoublePropertyValues.size(); s++)
      _mesh.Cell2DDoublePropertyValues[s].shrink_to_fit();

    _mesh.NumberCell3DVertices.shrink_to_fit();
    _mesh.NumberCell3DEdges.shrink_to_fit();
    _mesh.NumberCell3DFaces.shrink_to_fit();
    _mesh.Cell3DVertices.shrink_to_fit();
    _mesh.Cell3DEdges.shrink_to_fit();
    _mesh.Cell3DFaces.shrink_to_fit();
    _mesh.Cell3DMarkers.shrink_to_fit();
    _mesh.ActiveCell3D.shrink_to_fit();
    _mesh.Cell3DDoublePropertyIds.shrink_to_fit();
    _mesh.Cell3DDoublePropertySizes.shrink_to_fit();
    for (unsigned int s = 0; s < _mesh.Cell3DDoublePropertySizes.size(); s++)
      _mesh.Cell3DDoublePropertySizes[s].shrink_to_fit();
    _mesh.Cell3DDoublePropertyValues.shrink_to_fit();
    for (unsigned int s = 0; s < _mesh.Cell3DDoublePropertyValues.size(); s++)
      _mesh.Cell3DDoublePropertyValues[s].shrink_to_fit();
  }
  // ***************************************************************************
  string MeshMatricesDAO::ToString()
  {
    ostringstream converter;
    converter.precision(16);

    converter<< scientific<< "Dimension = "<< _mesh.Dimension<< ";"<< endl;
    converter<< scientific<< "NumberCell0D = "<< _mesh.NumberCell0D<< ";"<< endl;
    converter<< scientific<< "Cell0DCoordinates = "<< _mesh.Cell0DCoordinates<< ";"<< endl;
    converter<< scientific<< "Cell0DMarkers = "<< _mesh.Cell0DMarkers<< ";"<< endl;
    converter<< scientific<< "ActiveCell0D = "<< _mesh.ActiveCell0D<< ";"<< endl;
    converter<< scientific<< "UpdatedCell0Ds = "<< _mesh.UpdatedCell0Ds<< ";"<< endl;
    converter<< scientific<< "NumberCell0DNeighbourCell1D = "<< _mesh.NumberCell0DNeighbourCell1D<< ";"<< endl;
    converter<< scientific<< "Cell0DNeighbourCell1Ds = "<< _mesh.Cell0DNeighbourCell1Ds<< ";"<< endl;
    converter<< scientific<< "NumberCell0DNeighbourCell2D = "<< _mesh.NumberCell0DNeighbourCell2D<< ";"<< endl;
    converter<< scientific<< "Cell0DNeighbourCell2Ds = "<< _mesh.Cell0DNeighbourCell2Ds<< ";"<< endl;
    converter<< scientific<< "Cell0DDoublePropertyIds = "<< _mesh.Cell0DDoublePropertyIds<< ";"<< endl;
    converter<< scientific<< "Cell0DDoublePropertyIndices = "<< _mesh.Cell0DDoublePropertyIndices<< ";"<< endl;
    converter<< scientific<< "Cell0DDoublePropertySizes = "<< _mesh.Cell0DDoublePropertySizes<< ";"<< endl;
    converter<< scientific<< "Cell0DDoublePropertyValues = "<< _mesh.Cell0DDoublePropertyValues<< ";"<< endl;
    converter<< scientific<< "NumberCell1D = "<< _mesh.NumberCell1D<< ";"<< endl;
    converter<< scientific<< "Cell1DVertices = "<< _mesh.Cell1DVertices<< ";"<< endl;
    converter<< scientific<< "Cell1DMarkers = "<< _mesh.Cell1DMarkers<< ";"<< endl;
    converter<< scientific<< "Cell1DOriginalCell1Ds = "<< _mesh.Cell1DOriginalCell1Ds<< ";"<< endl;
    converter<< scientific<< "ActiveCell1D = "<< _mesh.ActiveCell1D<< ";"<< endl;
    converter<< scientific<< "UpdatedCell1Ds = "<< _mesh.UpdatedCell1Ds<< ";"<< endl;
    converter<< scientific<< "NumberCell1DNeighbourCell2D = "<< _mesh.NumberCell1DNeighbourCell2D<< ";"<< endl;
    converter<< scientific<< "Cell1DNeighbourCell2Ds = "<< _mesh.Cell1DNeighbourCell2Ds<< ";"<< endl;
    converter<< scientific<< "Cell1DDoublePropertyIds = "<< _mesh.Cell1DDoublePropertyIds<< ";"<< endl;
    converter<< scientific<< "Cell1DDoublePropertyIndices = "<< _mesh.Cell1DDoublePropertyIndices<< ";"<< endl;
    converter<< scientific<< "Cell1DDoublePropertySizes = "<< _mesh.Cell1DDoublePropertySizes<< ";"<< endl;
    converter<< scientific<< "Cell1DDoublePropertyValues = "<< _mesh.Cell1DDoublePropertyValues<< ";"<< endl;
    converter<< scientific<< "NumberCell2D = "<< _mesh.NumberCell2D<< ";"<< endl;
    converter<< scientific<< "NumberCell2DVertices = "<< _mesh.NumberCell2DVertices<< ";"<< endl;
    converter<< scientific<< "Cell2DVertices = "<< _mesh.Cell2DVertices<< ";"<< endl;
    converter<< scientific<< "NumberCell2DEdges = "<< _mesh.NumberCell2DEdges<< ";"<< endl;
    converter<< scientific<< "Cell2DEdges = "<< _mesh.Cell2DEdges<< ";"<< endl;
    converter<< scientific<< "NumberCell2DSubdivision = "<< _mesh.NumberCell2DSubdivision<< ";"<< endl;
    converter<< scientific<< "Cell2DSubdivision = "<< _mesh.Cell2DSubdivision<< ";"<< endl;
    converter<< scientific<< "Cell2DMarkers = "<< _mesh.Cell2DMarkers<< ";"<< endl;
    converter<< scientific<< "ActiveCell2D = "<< _mesh.ActiveCell2D<< ";"<< endl;
    converter<< scientific<< "Cell2DOriginalCell2Ds = "<< _mesh.Cell2DOriginalCell2Ds<< ";"<< endl;
    converter<< scientific<< "UpdatedCell2Ds = "<< _mesh.UpdatedCell2Ds<< ";"<< endl;
    converter<< scientific<< "NumberCell2DNeighbourCell3D = "<< _mesh.NumberCell2DNeighbourCell3D<< ";"<< endl;
    converter<< scientific<< "Cell2DNeighbourCell3Ds = "<< _mesh.Cell2DNeighbourCell3Ds<< ";"<< endl;
    converter<< scientific<< "Cell2DDoublePropertyIds = "<< _mesh.Cell2DDoublePropertyIds<< ";"<< endl;
    converter<< scientific<< "Cell2DDoublePropertyIndices = "<< _mesh.Cell2DDoublePropertyIndices<< ";"<< endl;
    converter<< scientific<< "Cell2DDoublePropertySizes = "<< _mesh.Cell2DDoublePropertySizes<< ";"<< endl;
    converter<< scientific<< "Cell2DDoublePropertyValues = "<< _mesh.Cell2DDoublePropertyValues<< ";"<< endl;
    converter<< scientific<< "NumberCell3D = "<< _mesh.NumberCell3D<< ";"<< endl;
    converter<< scientific<< "NumberCell3DVertices = "<< _mesh.NumberCell3DVertices<< ";"<< endl;
    converter<< scientific<< "Cell3DVertices = "<< _mesh.Cell3DVertices<< ";"<< endl;
    converter<< scientific<< "NumberCell3DEdges = "<< _mesh.NumberCell3DEdges<< ";"<< endl;
    converter<< scientific<< "Cell3DEdges = "<< _mesh.Cell3DEdges<< ";"<< endl;
    converter<< scientific<< "NumberCell3DFaces = "<< _mesh.NumberCell3DFaces<< ";"<< endl;
    converter<< scientific<< "Cell3DFaces = "<< _mesh.Cell3DFaces<< ";"<< endl;
    converter<< scientific<< "Cell3DMarkers = "<< _mesh.Cell3DMarkers<< ";"<< endl;
    converter<< scientific<< "ActiveCell3D = "<< _mesh.ActiveCell3D<< ";"<< endl;
    converter<< scientific<< "UpdatedCell3Ds = "<< _mesh.UpdatedCell3Ds<< ";"<< endl;
    converter<< scientific<< "Cell3DDoublePropertyIds = "<< _mesh.Cell3DDoublePropertyIds<< ";"<< endl;
    converter<< scientific<< "Cell3DDoublePropertyIndices = "<< _mesh.Cell3DDoublePropertyIndices<< ";"<< endl;
    converter<< scientific<< "Cell3DDoublePropertySizes = "<< _mesh.Cell3DDoublePropertySizes<< ";"<< endl;
    converter<< scientific<< "Cell3DDoublePropertyValues = "<< _mesh.Cell3DDoublePropertyValues<< ";"<< endl;

    return converter.str();
  }
  // ***************************************************************************
}
