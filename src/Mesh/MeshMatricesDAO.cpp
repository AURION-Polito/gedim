#include "MeshMatricesDAO.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
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
  void MeshMatricesDAO::AlignMapContainerHigherElements(map<unsigned int, Container>& elements,
                                                        const T& minElement,
                                                        const T& newElementInitialization)
  {
    map<unsigned int, Container> tempMap;
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
    _mesh.Cell1DAdjacency.conservativeResize(_mesh.NumberCell0D,
                                             _mesh.NumberCell0D);
    _mesh.NumberCell0DNeighbourCell1D.resize(_mesh.NumberCell0D + 1, 0);
    _mesh.NumberCell0DNeighbourCell2D.resize(_mesh.NumberCell0D + 1, 0);
    for (unsigned int p = 0; p < Cell0DNumberDoubleProperties(); p++)
      _mesh.Cell0DDoublePropertySizes[p].resize(_mesh.NumberCell0D + 1, 0);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell0DAppend(const unsigned int& numberCell0Ds)
  {
    unsigned int oldNumberCell0Ds = _mesh.NumberCell0D;
    Cell0DsInitialize(oldNumberCell0Ds + numberCell0Ds);

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

    _mesh.Cell1DAdjacency.conservativeResize(_mesh.NumberCell0D,
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
  MatrixXd MeshMatricesDAO::Cell0DCoordinates() const
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
      _mesh.UpdatedCell0Ds.insert(pair<unsigned int, set<unsigned int>>(cell0DIndex, {}));
    _mesh.UpdatedCell0Ds.at(cell0DIndex).insert(updatedCell0DIdex);
  }
  // ***************************************************************************
  bool MeshMatricesDAO::Cell0DUpdatedCell0Ds(const unsigned int& cell0DIndex,
                                             list<unsigned int>& updatedCell0DIds) const
  {
    Output::Assert(cell0DIndex < Cell0DTotalNumber());
    map<unsigned int, set<unsigned int>>::const_iterator iter = _mesh.UpdatedCell0Ds.find(cell0DIndex);

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
  void MeshMatricesDAO::Cell0DInitializeNeighbourCell1Ds(const unsigned int& cell0DIndex,
                                                         const unsigned int& numberNeighbourCell1Ds)
  {
    Output::Assert(cell0DIndex < Cell0DTotalNumber());

    ResizeNumberVectorWithNewNumberElements(_mesh.NumberCell0DNeighbourCell1D,
                                            _mesh.Cell0DNeighbourCell1Ds,
                                            _mesh.NumberCell0D,
                                            cell0DIndex,
                                            numberNeighbourCell1Ds,
                                            _mesh.NumberCell1D);
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
                                            _mesh.NumberCell2D);
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
  void MeshMatricesDAO::Cell0DInitializeDoublePropertyValues(const unsigned int& cell0DIndex,
                                                             const unsigned int& propertyIndex,
                                                             const unsigned int& porpertySize)
  {
    Output::Assert(cell0DIndex < Cell0DTotalNumber());
    Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());

    ResizeNumberVectorWithNewNumberElements(_mesh.Cell0DDoublePropertySizes[propertyIndex],
                                            _mesh.Cell0DDoublePropertyValues[propertyIndex],
                                            _mesh.NumberCell0D,
                                            cell0DIndex,
                                            porpertySize,
                                            0.0);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DsInitialize(const unsigned int& numberCell1Ds)
  {
    unsigned int oldNumberCell1D = _mesh.NumberCell1D;

    _mesh.NumberCell1D = numberCell1Ds;
    _mesh.Cell1DVertices.resize(2 * _mesh.NumberCell1D, 0);
    _mesh.Cell1DMarkers.resize(_mesh.NumberCell1D, 0);
    _mesh.ActiveCell1D.resize(_mesh.NumberCell1D, false);
    _mesh.NumberCell1DNeighbourCell2D.resize(_mesh.NumberCell1D + 1, 0);
    for (unsigned int p = 0; p < Cell1DNumberDoubleProperties(); p++)
      _mesh.Cell1DDoublePropertySizes[p].resize(_mesh.NumberCell1D + 1, 0);

    // update neighbours
    AlignContainerElements(_mesh.Cell0DNeighbourCell1Ds,
                           oldNumberCell1D,
                           _mesh.NumberCell1D);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell1DAppend(const unsigned int& numberCell1Ds)
  {
    unsigned int oldNumberCell1Ds = _mesh.NumberCell1D;
    Cell1DsInitialize(oldNumberCell1Ds + numberCell1Ds);

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
    _mesh.Cell1DMarkers.erase(std::next(_mesh.Cell1DMarkers.begin(), cell1DIndex));
    _mesh.ActiveCell1D.erase(std::next(_mesh.ActiveCell1D.begin(), cell1DIndex));
    _mesh.NumberCell1D--;

    AlignContainerHigherElements(_mesh.Cell0DNeighbourCell1Ds,
                                 cell1DIndex,
                                 _mesh.NumberCell1D);
    AlignContainerHigherElements(_mesh.Cell2DEdges,
                                 cell1DIndex,
                                 _mesh.NumberCell1D);
    AlignContainerHigherElements(_mesh.Cell3DEdges,
                                 cell1DIndex,
                                 _mesh.NumberCell1D);

    AlignMapContainerHigherElements(_mesh.UpdatedCell1Ds,
                                    cell1DIndex,
                                    _mesh.NumberCell1D);

    AlignSparseMatrixHigherElements(_mesh.Cell1DAdjacency,
                                    cell1DIndex + 1);
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
                                            _mesh.NumberCell2D);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell1DInsertUpdatedCell1D(const unsigned int& cell1DIndex,
                                                  const unsigned int& updatedCell1DIdex)
  {
    Output::Assert(cell1DIndex < Cell1DTotalNumber());
    Output::Assert(updatedCell1DIdex < Cell1DTotalNumber());

    if (!Cell1DHasUpdatedCell1Ds(cell1DIndex))
      _mesh.UpdatedCell1Ds.insert(pair<unsigned int, set<unsigned int>>(cell1DIndex, {}));
    _mesh.UpdatedCell1Ds.at(cell1DIndex).insert(updatedCell1DIdex);
  }
  // ***************************************************************************
  bool MeshMatricesDAO::Cell1DUpdatedCell1Ds(const unsigned int& cell1DIndex,
                                             list<unsigned int>& updatedCell1DIds) const
  {
    Output::Assert(cell1DIndex < Cell1DTotalNumber());
    map<unsigned int, set<unsigned int>>::const_iterator iter = _mesh.UpdatedCell1Ds.find(cell1DIndex);

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
  void MeshMatricesDAO::Cell1DInitializeDoublePropertyValues(const unsigned int& cell1DIndex,
                                                             const unsigned int& propertyIndex,
                                                             const unsigned int& porpertySize)
  {
    Output::Assert(cell1DIndex < Cell1DTotalNumber());
    Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());

    ResizeNumberVectorWithNewNumberElements(_mesh.Cell1DDoublePropertySizes[propertyIndex],
                                            _mesh.Cell1DDoublePropertyValues[propertyIndex],
                                            _mesh.NumberCell1D,
                                            cell1DIndex,
                                            porpertySize,
                                            0.0);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Cell2DsInitialize(const unsigned int& numberCell2Ds)
  {
    unsigned int oldNumberCell2D = _mesh.NumberCell2D;

    _mesh.NumberCell2D = numberCell2Ds;
    _mesh.NumberCell2DVertices.resize(_mesh.NumberCell2D + 1, 0);
    _mesh.NumberCell2DEdges.resize(_mesh.NumberCell2D + 1, 0);
    _mesh.NumberCell2DSubdivision.resize(_mesh.NumberCell2D + 1, 0);
    _mesh.Cell2DMarkers.resize(_mesh.NumberCell2D, 0);
    _mesh.ActiveCell2D.resize(_mesh.NumberCell2D, false);
    _mesh.Cell2DOriginalCell2Ds.resize(_mesh.NumberCell2D, _mesh.NumberCell2D);
    for (unsigned int p = 0; p < Cell2DNumberDoubleProperties(); p++)
      _mesh.Cell2DDoublePropertySizes[p].resize(_mesh.NumberCell2D + 1, 0);

    // update original cells
    AlignContainerElements(_mesh.Cell2DOriginalCell2Ds,
                           oldNumberCell2D,
                           _mesh.NumberCell2D);

    // update neighbours
    AlignContainerElements(_mesh.Cell0DNeighbourCell2Ds,
                           oldNumberCell2D,
                           _mesh.NumberCell2D);
    AlignContainerElements(_mesh.Cell1DNeighbourCell2Ds,
                           oldNumberCell2D,
                           _mesh.NumberCell2D);
  }
  // ***************************************************************************
  unsigned int MeshMatricesDAO::Cell2DAppend(const unsigned int& numberCell2Ds)
  {
    unsigned int oldNumberCell2Ds = _mesh.NumberCell2D;
    Cell2DsInitialize(oldNumberCell2Ds + numberCell2Ds);

    for (unsigned c = 0; c < numberCell2Ds; c++)
    {
      _mesh.NumberCell2DVertices[oldNumberCell2Ds + c + 1] = _mesh.NumberCell2DVertices[oldNumberCell2Ds];
      _mesh.NumberCell2DEdges[oldNumberCell2Ds + c + 1] = _mesh.NumberCell2DEdges[oldNumberCell2Ds];
      _mesh.NumberCell2DSubdivision[oldNumberCell2Ds + c + 1] = _mesh.NumberCell2DSubdivision[oldNumberCell2Ds];
    }

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

    _mesh.Cell2DOriginalCell2Ds.erase(std::next(_mesh.Cell2DOriginalCell2Ds.begin(), cell2DIndex));
    _mesh.Cell2DMarkers.erase(std::next(_mesh.Cell2DMarkers.begin(), cell2DIndex));
    _mesh.ActiveCell2D.erase(std::next(_mesh.ActiveCell2D.begin(), cell2DIndex));
    _mesh.NumberCell2D--;

    AlignContainerHigherElements(_mesh.Cell0DNeighbourCell2Ds,
                                 cell2DIndex,
                                 _mesh.NumberCell2D);
    AlignContainerHigherElements(_mesh.Cell1DNeighbourCell2Ds,
                                 cell2DIndex,
                                 _mesh.NumberCell2D);
    AlignContainerHigherElements(_mesh.Cell3DFaces,
                                 cell2DIndex,
                                 _mesh.NumberCell2D);
    AlignMapContainerHigherElements(_mesh.UpdatedCell2Ds,
                                    cell2DIndex,
                                    _mesh.NumberCell2D);
    AlignContainerHigherElements(_mesh.Cell2DOriginalCell2Ds,
                                 cell2DIndex,
                                 _mesh.NumberCell2D);
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
                                       const vector<unsigned int>& edgesCell0DIndices)
  {
    Cell2DInitializeEdges(cell2DIndex,
                          edgesCell0DIndices.size());
    for (unsigned int e = 0; e < edgesCell0DIndices.size(); e++)
      Cell2DInsertEdge(cell2DIndex,
                       e,
                       edgesCell0DIndices[e]);
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
  void MeshMatricesDAO::Cell2DInsertUpdatedCell2D(const unsigned int& cell2DIndex,
                                                  const unsigned int& updatedCell2DIdex)
  {
    Output::Assert(cell2DIndex < Cell2DTotalNumber());
    Output::Assert(updatedCell2DIdex < Cell2DTotalNumber());

    if (!Cell2DHasUpdatedCell2Ds(cell2DIndex))
      _mesh.UpdatedCell2Ds.insert(pair<unsigned int, set<unsigned int>>(cell2DIndex, {}));
    _mesh.UpdatedCell2Ds.at(cell2DIndex).insert(updatedCell2DIdex);
    _mesh.Cell2DOriginalCell2Ds[updatedCell2DIdex] = cell2DIndex;
  }
  // ***************************************************************************
  bool MeshMatricesDAO::Cell2DUpdatedCell2Ds(const unsigned int& cell2DIndex,
                                             list<unsigned int>& updatedCell2DIds) const
  {
    Output::Assert(cell2DIndex < Cell2DTotalNumber());
    map<unsigned int, set<unsigned int>>::const_iterator iter = _mesh.UpdatedCell2Ds.find(cell2DIndex);

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
  void MeshMatricesDAO::Cell2DInitializeDoublePropertyValues(const unsigned int& cell2DIndex,
                                                             const unsigned int& propertyIndex,
                                                             const unsigned int& porpertySize)
  {
    Output::Assert(cell2DIndex < Cell2DTotalNumber());
    Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());

    ResizeNumberVectorWithNewNumberElements(_mesh.Cell2DDoublePropertySizes[propertyIndex],
                                            _mesh.Cell2DDoublePropertyValues[propertyIndex],
                                            _mesh.NumberCell2D,
                                            cell2DIndex,
                                            porpertySize,
                                            0.0);
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
    Cell3DsInitialize(oldNumberCell3Ds + numberCell3Ds);

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

    AlignMapContainerHigherElements(_mesh.UpdatedCell3Ds,
                                    cell3DIndex,
                                    _mesh.NumberCell3D);
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
  void MeshMatricesDAO::Cell3DInsertUpdatedCell3D(const unsigned int& cell3DIndex,
                                                  const unsigned int& updatedCell3DIdex)
  {
    Output::Assert(cell3DIndex < Cell3DTotalNumber());
    Output::Assert(updatedCell3DIdex < Cell3DTotalNumber());

    if (!Cell3DHasUpdatedCell3Ds(cell3DIndex))
      _mesh.UpdatedCell3Ds.insert(pair<unsigned int, set<unsigned int>>(cell3DIndex, {}));
    _mesh.UpdatedCell3Ds.at(cell3DIndex).insert(updatedCell3DIdex);
  }
  // ***************************************************************************
  bool MeshMatricesDAO::Cell3DUpdatedCell3Ds(const unsigned int& cell3DIndex,
                                             list<unsigned int>& updatedCell3DIds) const
  {
    Output::Assert(cell3DIndex < Cell3DTotalNumber());
    map<unsigned int, set<unsigned int>>::const_iterator iter = _mesh.UpdatedCell3Ds.find(cell3DIndex);

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
  void MeshMatricesDAO::Cell3DInitializeDoublePropertyValues(const unsigned int& cell3DIndex,
                                                             const unsigned int& propertyIndex,
                                                             const unsigned int& porpertySize)
  {
    Output::Assert(cell3DIndex < Cell3DTotalNumber());
    Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());

    ResizeNumberVectorWithNewNumberElements(_mesh.Cell3DDoublePropertySizes[propertyIndex],
                                            _mesh.Cell3DDoublePropertyValues[propertyIndex],
                                            _mesh.NumberCell3D,
                                            cell3DIndex,
                                            porpertySize,
                                            0.0);
  }
  // ***************************************************************************
  void MeshMatricesDAO::Compress()
  {
    _mesh.Cell0DCoordinates.shrink_to_fit();
    _mesh.Cell0DMarkers.shrink_to_fit();
    _mesh.ActiveCell0D.shrink_to_fit();
    _mesh.Cell0DDoublePropertyIds.shrink_to_fit();
    _mesh.Cell0DDoublePropertySizes.shrink_to_fit();
    for (unsigned int s = 0; s < _mesh.Cell0DDoublePropertySizes.size(); s++)
      _mesh.Cell0DDoublePropertySizes[s].shrink_to_fit();
    _mesh.Cell0DDoublePropertyValues.shrink_to_fit();
    for (unsigned int s = 0; s < _mesh.Cell0DDoublePropertyValues.size(); s++)
      _mesh.Cell0DDoublePropertyValues[s].shrink_to_fit();

    _mesh.Cell1DVertices.shrink_to_fit();
    _mesh.Cell1DAdjacency.makeCompressed();
    _mesh.NumberCell1DNeighbourCell2D.shrink_to_fit();
    _mesh.Cell1DNeighbourCell2Ds.shrink_to_fit();
    _mesh.Cell1DMarkers.shrink_to_fit();
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

    converter<< scientific<< "Cell1DAdjacency.resize("<< _mesh.Cell1DAdjacency.rows()<< ", "<< _mesh.Cell1DAdjacency.rows()<< ");"<< endl;
    vector<unsigned int> tripletI(_mesh.Cell1DAdjacency.nonZeros());
    vector<unsigned int> tripletJ(_mesh.Cell1DAdjacency.nonZeros());
    vector<unsigned int> tripletValue(_mesh.Cell1DAdjacency.nonZeros());
    unsigned int element = 0;
    for (int k=0; k< _mesh.Cell1DAdjacency.outerSize(); ++k)
    {
      for (SparseMatrix<unsigned int>::InnerIterator it(_mesh.Cell1DAdjacency, k); it; ++it)
      {
        tripletValue[element] = it.value();
        tripletI[element] = it.row();
        tripletJ[element] = it.col();
        element++;
      }
    }

    converter<< scientific<< "Cell1DAdjacency.reserve("<< _mesh.Cell1DAdjacency.nonZeros()<< ");"<< endl;
    for (unsigned int e = 0; e < _mesh.Cell1DAdjacency.nonZeros(); e++)
    {
      converter<< scientific<< "Cell1DAdjacency.insert("<< tripletI[e];
      converter<< scientific<< ", "<<tripletJ[e]<< ") = ";
      converter<< scientific<< tripletValue[e]<< ";"<< endl;
    }
    converter<< scientific<< "Cell1DAdjacency.makeCompressed();"<< endl;

    converter<< scientific<< "Cell1DMarkers = "<< _mesh.Cell1DMarkers<< ";"<< endl;
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
