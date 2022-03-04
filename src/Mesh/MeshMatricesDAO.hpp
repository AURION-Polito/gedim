#ifndef __MeshMatricesWrapper_H
#define __MeshMatricesWrapper_H

#include "IOUtilities.hpp"
#include "Eigen"
#include "MeshMatrices.hpp"
#include "IMeshDAO.hpp"

using namespace std;

namespace Gedim
{
  class MeshMatricesDAO final : public IMeshDAO
  {
    private:
      MeshMatrices& _mesh;

      /// \brief for each i,j element of sparse matrix A if A[i,j] > minElement then A[i, j]--
      /// if A[i, j] == minElement the A[i, j] = 0
      /// \param matrix the sparse matrix A
      /// \param minElement the minElement
      /// \param newElementInitialization the new element initialization
      template<typename T>
      void AlignSparseMatrixHigherElements(Eigen::SparseMatrix<T>& matrix,
                                           const T& minElement);

      /// \brief for each i element of container on each map key v if v[i] > minElement then v[i]--
      /// if v[i] == minElement the v[i] = newElementInitialization
      /// \param elements the container map v
      /// \param minElement the minElement
      /// \param newElementInitialization the new element initialization
      template<class Container, class T>
      void AlignMapContainerHigherElements(map<unsigned int, Container>& elements,
                                           const T& minElement,
                                           const T& newElementInitialization);

      /// \brief for each i element of container v if v[i] > minElement then v[i]--
      /// if v[i] == minElement the v[i] = newElementInitialization
      /// \param elements the container v
      /// \param minElement the minElement
      /// \param newElementInitialization the new element initialization
      template<class Container, class T>
      void AlignContainerHigherElements(Container& elements,
                                        const T& minElement,
                                        const T& newElementInitialization);

      /// \brief for each i element of container v if v[i] == element then v[i] = newElementInitialization
      /// \param elements the container v
      /// \param element the element
      /// \param newElementInitialization the new element initialization
      template<class Container, class T>
      void AlignContainerElements(Container& elements,
                                  const T& element,
                                  const T& newElementInitialization);

      template<typename T>
      void ResizeNumberVectorWithNewNumberElements(vector<unsigned int>& numberElementVector,
                                                   vector<T>& elementVector,
                                                   const unsigned int& numberElements,
                                                   const unsigned int& vectorIndex,
                                                   const unsigned int& newNumberElements,
                                                   const T& newElementInitialization = T());

    public:
      MeshMatricesDAO(MeshMatrices& mesh);
      ~MeshMatricesDAO();

      inline void InitializeDimension(const unsigned int& dimension)
      { _mesh.Dimension = dimension; }
      inline unsigned int Dimension() const
      { return _mesh.Dimension; }

      void Cell0DsInitialize(const unsigned int& numberCell0Ds);
      unsigned int Cell0DAppend(const unsigned int& numberCell0Ds);
      void Cell0DRemove(const unsigned int& cell0DIndex);

      void Cell0DInsertCoordinates(const unsigned int& cell0DIndex,
                                   const Eigen::Vector3d& coordinates);
      void Cell0DsInsertCoordinates(const Eigen::MatrixXd& coordinates);
      inline void Cell0DSetMarker(const unsigned int& cell0DIndex,
                                  const unsigned int& marker)
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        _mesh.Cell0DMarkers[cell0DIndex] = marker;
      }
      inline void Cell0DSetState(const unsigned int& cell0DIndex,
                                 const bool& state)
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        _mesh.ActiveCell0D[cell0DIndex] = state;
      }

      inline unsigned int Cell0DTotalNumber() const
      { return _mesh.NumberCell0D; }
      inline double Cell0DCoordinateX(const unsigned int& cell0DIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DCoordinates[3 * cell0DIndex];
      }
      inline double Cell0DCoordinateY(const unsigned int& cell0DIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DCoordinates[3 * cell0DIndex + 1];
      }
      inline double Cell0DCoordinateZ(const unsigned int& cell0DIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DCoordinates[3 * cell0DIndex + 2];
      }
      inline Eigen::Vector3d Cell0DCoordinates(const unsigned int& cell0DIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return Eigen::Vector3d(Cell0DCoordinateX(cell0DIndex),
                               Cell0DCoordinateY(cell0DIndex),
                               Cell0DCoordinateZ(cell0DIndex));
      }
      Eigen::MatrixXd Cell0DCoordinates() const;
      inline unsigned int Cell0DMarker(const unsigned int& cell0DIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DMarkers[cell0DIndex];
      }
      inline bool Cell0DIsActive(const unsigned int& cell0DIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.ActiveCell0D[cell0DIndex];
      }

      inline bool Cell0DHasUpdatedCell0Ds(const unsigned int& cell0DIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.UpdatedCell0Ds.find(cell0DIndex) != _mesh.UpdatedCell0Ds.end();
      }
      inline unsigned int Cell0DNumberUpdatedCell0Ds(const unsigned int& cell0DIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.UpdatedCell0Ds.at(cell0DIndex).size();
      }
      inline bool Cell0DHasUpdatedCell0D(const unsigned int& cell0DIndex,
                                         const unsigned int& updatedCell0DIdex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Output::Assert(updatedCell0DIdex < Cell0DTotalNumber());
        return _mesh.UpdatedCell0Ds.at(cell0DIndex).find(updatedCell0DIdex) != _mesh.UpdatedCell0Ds.at(cell0DIndex).end();
      }
      void Cell0DInsertUpdatedCell0D(const unsigned int& cell0DIndex,
                                     const unsigned int& updatedCell0DIdex);
      bool Cell0DUpdatedCell0Ds(const unsigned int& cell0DIndex,
                                list<unsigned int>& updatedCell0DIds) const;

      inline void Cell0DSetId(const unsigned int& cell0DIndex,
                              const unsigned int& id) { return; }
      inline unsigned int Cell0DId(const unsigned int& cell0DIndex) const { return cell0DIndex; }
      inline void Cell0DInitializeNeighbourCell1Ds(const unsigned int& cell0DIndex,
                                                   const unsigned int& numberNeighbourCell1Ds);
      inline void Cell0DInsertNeighbourCell1D(const unsigned int& cell0DIndex,
                                              const unsigned int& neighbourIndex,
                                              const unsigned int& neigbourCell1DIndex)
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell1D(cell0DIndex));
        Output::Assert(neigbourCell1DIndex < Cell1DTotalNumber());

        _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] +
            neighbourIndex] = neigbourCell1DIndex;
      }
      inline unsigned int Cell0DNumberNeighbourCell1D(const unsigned int& cell0DIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.NumberCell0DNeighbourCell1D[cell0DIndex + 1] -
            _mesh.NumberCell0DNeighbourCell1D[cell0DIndex];
      }
      inline unsigned int Cell0DNeighbourCell1D(const unsigned int& cell0DIndex,
                                                const unsigned int& neighbourIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell1D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] +
            neighbourIndex];
      }
      inline bool Cell0DHasNeighbourCell1D(const unsigned int& cell0DIndex,
                                           const unsigned int& neighbourIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell1D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] +
            neighbourIndex] < _mesh.NumberCell1D;
      }
      inline void Cell0DResetNeighbourCell1D(const unsigned int& cell0DIndex,
                                             const unsigned int& neighbourIndex)
      {
        _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] +
            neighbourIndex] = _mesh.NumberCell1D;
      }

      inline void Cell0DInitializeNeighbourCell2Ds(const unsigned int& cell0DIndex,
                                                   const unsigned int& numberNeighbourCell2Ds);
      inline void Cell0DInsertNeighbourCell2D(const unsigned int& cell0DIndex,
                                              const unsigned int& neighbourIndex,
                                              const unsigned int& neigbourCell2DIndex)
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell2D(cell0DIndex));
        Output::Assert(neigbourCell2DIndex < Cell2DTotalNumber());

        _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] +
            neighbourIndex] = neigbourCell2DIndex;
      }
      inline unsigned int Cell0DNumberNeighbourCell2D(const unsigned int& cell0DIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.NumberCell0DNeighbourCell2D[cell0DIndex + 1] -
            _mesh.NumberCell0DNeighbourCell2D[cell0DIndex];
      }
      inline unsigned int Cell0DNeighbourCell2D(const unsigned int& cell0DIndex,
                                                const unsigned int& neighbourIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell2D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] +
            neighbourIndex];
      }
      inline bool Cell0DHasNeighbourCell2D(const unsigned int& cell0DIndex,
                                           const unsigned int& neighbourIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell2D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] +
            neighbourIndex] < _mesh.NumberCell2D;
      }
      inline void Cell0DResetNeighbourCell2D(const unsigned int& cell0DIndex,
                                             const unsigned int& neighbourIndex)
      {
        _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] +
            neighbourIndex] = _mesh.NumberCell2D;
      }
      inline void Cell0DInitializeNeighbourCell3Ds(const unsigned int& cell0DIndex,
                                                   const unsigned int& numberNeighbourCell3Ds) { return; }
      inline void Cell0DInsertNeighbourCell3D(const unsigned int& cell0DIndex,
                                              const unsigned int& neighbourIndex,
                                              const unsigned int& neigbourCell3DIndex) { throw runtime_error("Not implemented"); }
      inline unsigned int Cell0DNumberNeighbourCell3D(const unsigned int& cell0DIndex) const { return 0; }
      inline unsigned int Cell0DNeighbourCell3D(const unsigned int& cell0DIndex,
                                                const unsigned int& neighbourIndex) const { throw runtime_error("Not implemented"); }
      inline bool Cell0DHasNeighbourCell3D(const unsigned int& cell0DIndex,
                                           const unsigned int& neighbourIndex) const { throw runtime_error("Not implemented"); }
      inline void Cell0DResetNeighbourCell3D(const unsigned int& cell0DIndex,
                                             const unsigned int& neighbourIndex)
      {
        throw runtime_error("Not implemented");
      }

      void Cell0DInitializeDoubleProperties(const unsigned int& numberDoubleProperties);
      unsigned int Cell0DAddDoubleProperty(const string& propertyId);
      void Cell0DInitializeDoublePropertyValues(const unsigned int& cell0DIndex,
                                                const unsigned int& propertyIndex,
                                                const unsigned int& porpertySize);
      inline void Cell0DInsertDoublePropertyValue(const unsigned int& cell0DIndex,
                                                  const unsigned int& propertyIndex,
                                                  const unsigned int& propertyValueIndex,
                                                  const double& propertyValue)
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());

        _mesh.Cell0DDoublePropertyValues[propertyIndex][_mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex] +
            propertyValueIndex] = propertyValue;
      }
      inline unsigned int Cell0DNumberDoubleProperties() const
      {
        return _mesh.Cell0DDoublePropertyIds.size();
      }
      inline string Cell0DDoublePropertyId(const unsigned int& propertyIndex) const
      {
        return _mesh.Cell0DDoublePropertyIds[propertyIndex];
      }
      inline bool Cell0DDoublePropertyExists(const string& propertyId) const
      {
        return _mesh.Cell0DDoublePropertyIndices.find(propertyId) !=
            _mesh.Cell0DDoublePropertyIndices.end();
      }
      inline unsigned int Cell0DDoublePropertyIndex(const string& propertyId) const
      {
        Output::Assert(Cell0DDoublePropertyExists(propertyId));
        return _mesh.Cell0DDoublePropertyIndices.at(propertyId);
      }
      inline unsigned int Cell0DDoublePropertySize(const unsigned int& cell0DIndex,
                                                   const unsigned int& propertyIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());
        return _mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex + 1] -
            _mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex];
      }
      inline double Cell0DDoublePropertyValue(const unsigned int& cell0DIndex,
                                              const unsigned int& propertyIndex,
                                              const unsigned int& propertyValueIndex) const
      {
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());
        Output::Assert(propertyValueIndex < Cell0DDoublePropertySize(cell0DIndex,
                                                                     propertyIndex));

        return _mesh.Cell0DDoublePropertyValues[propertyIndex][_mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex] +
            propertyValueIndex];
      }

      void Cell1DsInitialize(const unsigned int& numberCell1Ds);
      unsigned int Cell1DAppend(const unsigned int& numberCell1Ds);
      void Cell1DRemove(const unsigned int& cell1DIndex);
      inline void Cell1DInsertExtremes(const unsigned int& cell1DIndex,
                                       const unsigned int& originCell0DIndex,
                                       const unsigned int& endCell0DIndex)
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Output::Assert(originCell0DIndex < Cell0DTotalNumber());
        Output::Assert(endCell0DIndex < Cell0DTotalNumber());
        _mesh.Cell1DVertices[2 * cell1DIndex] = originCell0DIndex;
        _mesh.Cell1DVertices[2 * cell1DIndex + 1] = endCell0DIndex;
        _mesh.Cell1DAdjacency.insert(originCell0DIndex,
                                     endCell0DIndex) = cell1DIndex + 1;
        _mesh.Cell1DAdjacency.makeCompressed();
      }
      inline bool Cell1DExists(const unsigned int& originCell0DIndex,
                               const unsigned int& endCell0DIndex) const
      {
        Output::Assert(originCell0DIndex < Cell0DTotalNumber());
        Output::Assert(endCell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell1DAdjacency.coeff(originCell0DIndex, endCell0DIndex) > 0;
      }

      inline unsigned int Cell1DByExtremes(const unsigned int& originCell0DIndex,
                                           const unsigned int& endCell0DIndex) const
      {
        Output::Assert(Cell1DExists(originCell0DIndex, endCell0DIndex));
        return _mesh.Cell1DAdjacency.coeff(originCell0DIndex, endCell0DIndex) - 1;
      }

      inline void Cell1DSetMarker(const unsigned int& cell1DIndex,
                                  const unsigned int& marker)
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        _mesh.Cell1DMarkers[cell1DIndex] = marker;
      }
      inline void Cell1DSetState(const unsigned int& cell1DIndex,
                                 const bool& state)
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        _mesh.ActiveCell1D[cell1DIndex] = state;
      }
      void Cell1DInitializeNeighbourCell2Ds(const unsigned int& cell1DIndex,
                                            const unsigned int& numberNeighbourCell2Ds);
      inline void Cell1DInsertNeighbourCell2D(const unsigned int& cell1DIndex,
                                              const unsigned int& neighbourIndex,
                                              const unsigned int& neigbourCell2DIndex)
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell2D(cell1DIndex));
        Output::Assert(neigbourCell2DIndex < Cell2DTotalNumber());

        _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] +
            neighbourIndex] = neigbourCell2DIndex;
      }
      inline unsigned int Cell1DTotalNumber() const
      { return _mesh.NumberCell1D; }
      inline unsigned int Cell1DVertex(const unsigned int& cell1DIndex,
                                       const unsigned int& vertexIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Output::Assert(vertexIndex < 2);
        return _mesh.Cell1DVertices[2 * cell1DIndex + vertexIndex];
      }
      inline unsigned int Cell1DOrigin(const unsigned int& cell1DIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DVertices[2 * cell1DIndex];
      }
      inline unsigned int Cell1DEnd(const unsigned int& cell1DIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DVertices[2 * cell1DIndex + 1];
      }
      inline Eigen::Vector3d Cell1DOriginCoordinates(const unsigned int& cell1DIndex) const
      {
        return Cell0DCoordinates(Cell1DOrigin(cell1DIndex));
      }
      inline Eigen::Vector3d Cell1DEndCoordinates(const unsigned int& cell1DIndex) const
      {
        return Cell0DCoordinates(Cell1DEnd(cell1DIndex));
      }
      inline unsigned int Cell1DNumberNeighbourCell2D(const unsigned int& cell1DIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.NumberCell1DNeighbourCell2D[cell1DIndex + 1] -
            _mesh.NumberCell1DNeighbourCell2D[cell1DIndex];
      }
      inline unsigned int Cell1DNeighbourCell2D(const unsigned int& cell1DIndex,
                                                const unsigned int& neighbourIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell2D(cell1DIndex));
        return _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] +
            neighbourIndex];
      }
      bool Cell1DHasNeighbourCell2D(const unsigned int& cell1DIndex,
                                    const unsigned int& neighbourIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell2D(cell1DIndex));
        return _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] +
            neighbourIndex] < _mesh.NumberCell2D;
      }
      inline void Cell1DResetNeighbourCell2D(const unsigned int& cell1DIndex,
                                             const unsigned int& neighbourIndex)
      {
        _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] +
            neighbourIndex] = _mesh.NumberCell2D;
      }

      inline unsigned int Cell1DMarker(const unsigned int& cell1DIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DMarkers[cell1DIndex];
      }
      inline bool Cell1DIsActive(const unsigned int& cell1DIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.ActiveCell1D[cell1DIndex];
      }
      inline bool Cell1DHasUpdatedCell1Ds(const unsigned int& cell1DIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.UpdatedCell1Ds.find(cell1DIndex) != _mesh.UpdatedCell1Ds.end();
      }
      inline unsigned int Cell1DNumberUpdatedCell1Ds(const unsigned int& cell1DIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.UpdatedCell1Ds.at(cell1DIndex).size();
      }
      inline bool Cell1DHasUpdatedCell1D(const unsigned int& cell1DIndex,
                                         const unsigned int& updatedCell1DIdex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Output::Assert(updatedCell1DIdex < Cell1DTotalNumber());
        return _mesh.UpdatedCell1Ds.at(cell1DIndex).find(updatedCell1DIdex) != _mesh.UpdatedCell1Ds.at(cell1DIndex).end();
      }
      void Cell1DInsertUpdatedCell1D(const unsigned int& cell1DIndex,
                                     const unsigned int& updatedCell1DIdex);
      bool Cell1DUpdatedCell1Ds(const unsigned int& cell1DIndex,
                                list<unsigned int>& updatedCell1DIds) const;
      void Cell1DInitializeDoubleProperties(const unsigned int& numberDoubleProperties);

      inline void Cell1DSetId(const unsigned int& cell1DIndex,
                              const unsigned int& id) { return; }
      inline unsigned int Cell1DId(const unsigned int& cell1DIndex) const { return cell1DIndex; }
      inline void Cell1DInitializeNeighbourCell3Ds(const unsigned int& cell1DIndex,
                                                   const unsigned int& numberNeighbourCell3Ds) {  return; }
      inline void Cell1DInsertNeighbourCell3D(const unsigned int& cell1DIndex,
                                              const unsigned int& neighbourIndex, const unsigned int& neigbourCell3DIndex) { throw runtime_error("Not implemented"); }
      inline unsigned int Cell1DNumberNeighbourCell3D(const unsigned int& cell1DIndex) const { return 0; }
      inline unsigned int Cell1DNeighbourCell3D(const unsigned int& cell1DIndex,
                                                const unsigned int& neighbourIndex) const { throw runtime_error("Not implemented"); }
      inline bool Cell1DHasNeighbourCell3D(const unsigned int& cell1DIndex,
                                           const unsigned int& neighbourIndex) const { throw runtime_error("Not implemented"); }
      inline void Cell1DResetNeighbourCell3D(const unsigned int& cell1DIndex,
                                             const unsigned int& neighbourIndex)
      {
        throw runtime_error("Not implemented");
      }

      unsigned int Cell1DAddDoubleProperty(const string& propertyId);
      void Cell1DInitializeDoublePropertyValues(const unsigned int& cell1DIndex,
                                                const unsigned int& propertyIndex,
                                                const unsigned int& porpertySize);
      inline void Cell1DInsertDoublePropertyValue(const unsigned int& cell1DIndex,
                                                  const unsigned int& propertyIndex,
                                                  const unsigned int& propertyValueIndex,
                                                  const double& propertyValue)
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());

        _mesh.Cell1DDoublePropertyValues[propertyIndex][_mesh.Cell1DDoublePropertySizes[propertyIndex][cell1DIndex] +
            propertyValueIndex] = propertyValue;
      }
      inline unsigned int Cell1DNumberDoubleProperties() const
      {
        return _mesh.Cell1DDoublePropertyIds.size();
      }
      inline string Cell1DDoublePropertyId(const unsigned int& propertyIndex) const
      {
        return _mesh.Cell1DDoublePropertyIds[propertyIndex];
      }
      inline bool Cell1DDoublePropertyExists(const string& propertyId) const
      {
        return _mesh.Cell1DDoublePropertyIndices.find(propertyId) !=
            _mesh.Cell1DDoublePropertyIndices.end();
      }
      inline unsigned int Cell1DDoublePropertyIndex(const string& propertyId) const
      {
        Output::Assert(Cell1DDoublePropertyExists(propertyId));
        return _mesh.Cell1DDoublePropertyIndices.at(propertyId);
      }
      inline unsigned int Cell1DDoublePropertySize(const unsigned int& cell1DIndex,
                                                   const unsigned int& propertyIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());
        return _mesh.Cell1DDoublePropertySizes[propertyIndex][cell1DIndex + 1] -
            _mesh.Cell1DDoublePropertySizes[propertyIndex][cell1DIndex];
      }
      inline double Cell1DDoublePropertyValue(const unsigned int& cell1DIndex,
                                              const unsigned int& propertyIndex,
                                              const unsigned int& propertyValueIndex) const
      {
        Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());
        Output::Assert(propertyValueIndex < Cell1DDoublePropertySize(cell1DIndex,
                                                                     propertyIndex));

        return _mesh.Cell1DDoublePropertyValues[propertyIndex][_mesh.Cell1DDoublePropertySizes[propertyIndex][cell1DIndex] +
            propertyValueIndex];
      }


      void Cell2DsInitialize(const unsigned int& numberCell2Ds);
      unsigned int Cell2DAppend(const unsigned int& numberCell2Ds);
      void Cell2DRemove(const unsigned int& cell2DIndex);
      inline void Cell2DSetMarker(const unsigned int& cell2DIndex,
                                  const unsigned int& marker)
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        _mesh.Cell2DMarkers[cell2DIndex] = marker;
      }
      inline void Cell2DSetState(const unsigned int& cell2DIndex,
                                 const bool& state)
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        _mesh.ActiveCell2D[cell2DIndex] = state;
      }
      void Cell2DInitializeVertices(const unsigned int& cell2DIndex,
                                    const unsigned int& numberCell2DVertices);
      void Cell2DInitializeEdges(const unsigned int& cell2DIndex,
                                 const unsigned int& numberCell2DEdges);
      inline void Cell2DInsertVertex(const unsigned int& cell2DIndex,
                                     const unsigned int& vertexIndex,
                                     const unsigned int& vertexCell0DIndex)
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Output::Assert(vertexIndex < Cell2DNumberVertices(cell2DIndex));
        Output::Assert(vertexCell0DIndex < Cell0DTotalNumber());
        _mesh.Cell2DVertices[_mesh.NumberCell2DVertices[cell2DIndex] +
            vertexIndex] = vertexCell0DIndex;
      }
      inline void Cell2DInsertEdge(const unsigned int& cell2DIndex,
                                   const unsigned int& edgeIndex,
                                   const unsigned int& edgeCell1DIndex)
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Output::Assert(edgeIndex < Cell2DNumberEdges(cell2DIndex));
        Output::Assert(edgeCell1DIndex < Cell1DTotalNumber());
        _mesh.Cell2DEdges[_mesh.NumberCell2DEdges[cell2DIndex] +
            edgeIndex] = edgeCell1DIndex;
      }
      void Cell2DAddVertices(const unsigned int& cell2DIndex,
                             const vector<unsigned int>& verticesCell0DIndices);
      void Cell2DAddEdges(const unsigned int& cell2DIndex,
                          const vector<unsigned int>& edgesCell0DIndices);

      inline unsigned int Cell2DTotalNumber() const
      { return _mesh.NumberCell2D; }
      inline unsigned int Cell2DNumberVertices(const unsigned int& cell2DIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.NumberCell2DVertices[cell2DIndex + 1] -
            _mesh.NumberCell2DVertices[cell2DIndex];
      }
      inline unsigned int Cell2DNumberEdges(const unsigned int& cell2DIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.NumberCell2DEdges[cell2DIndex + 1] - _mesh.NumberCell2DEdges[cell2DIndex];
      }
      inline unsigned int Cell2DVertex(const unsigned int& cell2DIndex,
                                       const unsigned int& vertexIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Output::Assert(vertexIndex < Cell2DNumberVertices(cell2DIndex));
        return _mesh.Cell2DVertices[_mesh.NumberCell2DVertices[cell2DIndex] + vertexIndex];
      }
      inline Eigen::Vector3d Cell2DVertexCoordinates(const unsigned int& cell2DIndex,
                                                     const unsigned int& vertexIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Output::Assert(vertexIndex < Cell2DNumberVertices(cell2DIndex));
        return Cell0DCoordinates(Cell2DVertex(cell2DIndex, vertexIndex));
      }
      Eigen::MatrixXd Cell2DVerticesCoordinates(const unsigned int& cell2DIndex) const;
      inline unsigned int Cell2DEdge(const unsigned int& cell2DIndex,
                                     const unsigned int& edgeIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Output::Assert(edgeIndex < Cell2DNumberEdges(cell2DIndex));
        return _mesh.Cell2DEdges[_mesh.NumberCell2DEdges[cell2DIndex] + edgeIndex];
      }
      inline unsigned int Cell2DMarker(const unsigned int& cell2DIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.Cell2DMarkers[cell2DIndex];
      }
      inline bool Cell2DIsActive(const unsigned int& cell2DIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.ActiveCell2D[cell2DIndex];
      }

      inline bool Cell2DHasUpdatedCell2Ds(const unsigned int& cell2DIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.UpdatedCell2Ds.find(cell2DIndex) != _mesh.UpdatedCell2Ds.end();
      }
      inline unsigned int Cell2DNumberUpdatedCell2Ds(const unsigned int& cell2DIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.UpdatedCell2Ds.at(cell2DIndex).size();
      }
      inline bool Cell2DHasUpdatedCell2D(const unsigned int& cell2DIndex,
                                         const unsigned int& updatedCell2DIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Output::Assert(updatedCell2DIndex < Cell2DTotalNumber());
        return _mesh.UpdatedCell2Ds.at(cell2DIndex).find(updatedCell2DIndex) != _mesh.UpdatedCell2Ds.at(cell2DIndex).end();
      }
      void Cell2DInsertUpdatedCell2D(const unsigned int& cell2DIndex,
                                     const unsigned int& updatedCell2DIdex);
      inline bool Cell2DHasOriginalCell2D(const unsigned int& updatedCell2DIndex) const
      {
        Output::Assert(updatedCell2DIndex < Cell2DTotalNumber());
        return _mesh.Cell2DOriginalCell2Ds.at(updatedCell2DIndex) < _mesh.NumberCell2D;
      }
      inline unsigned int Cell2DOriginalCell2D(const unsigned int& updatedCell2DIndex) const
      {
        Output::Assert(updatedCell2DIndex < Cell2DTotalNumber());
        return _mesh.Cell2DOriginalCell2Ds.at(updatedCell2DIndex);
      }
      bool Cell2DUpdatedCell2Ds(const unsigned int& cell2DIndex,
                                list<unsigned int>& updatedCell2DIds) const;

      inline void Cell2DSetId(const unsigned int& cell2DIndex,
                              const unsigned int& id) { return; }
      inline unsigned int Cell2DId(const unsigned int& cell2DIndex) const { return cell2DIndex; }
      inline void Cell2DInitializeNeighbourCell3Ds(const unsigned int& cell2DIndex,
                                                   const unsigned int& numberNeighbourCell3Ds) { return; }
      inline void Cell2DInsertNeighbourCell3D(const unsigned int& cell2DIndex,
                                              const unsigned int& neighbourIndex,
                                              const unsigned int& neigbourCell3DIndex) { throw runtime_error("Not implemented"); }
      inline unsigned int Cell2DNumberNeighbourCell3D(const unsigned int& cell2DIndex) const { return 0; }
      inline unsigned int Cell2DNeighbourCell3D(const unsigned int& cell2DIndex,
                                                const unsigned int& neighbourIndex) const { throw runtime_error("Not implemented"); }
      inline bool Cell2DHasNeighbourCell3D(const unsigned int& cell2DIndex,
                                           const unsigned int& neighbourIndex) const { throw runtime_error("Not implemented"); }
      inline void Cell2DResetNeighbourCell3D(const unsigned int& cell2DIndex,
                                             const unsigned int& neighbourIndex)
      {
        throw runtime_error("Not implemented");
      }
      void Cell2DInitializeDoubleProperties(const unsigned int& numberDoubleProperties);
      unsigned int Cell2DAddDoubleProperty(const string& propertyId);
      void Cell2DInitializeDoublePropertyValues(const unsigned int& cell2DIndex,
                                                const unsigned int& propertyIndex,
                                                const unsigned int& porpertySize);
      inline void Cell2DInsertDoublePropertyValue(const unsigned int& cell2DIndex,
                                                  const unsigned int& propertyIndex,
                                                  const unsigned int& propertyValueIndex,
                                                  const double& propertyValue)
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());

        _mesh.Cell2DDoublePropertyValues[propertyIndex][_mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex] +
            propertyValueIndex] = propertyValue;
      }
      inline unsigned int Cell2DNumberDoubleProperties() const
      {
        return _mesh.Cell2DDoublePropertyIds.size();
      }
      inline string Cell2DDoublePropertyId(const unsigned int& propertyIndex) const
      {
        return _mesh.Cell2DDoublePropertyIds[propertyIndex];
      }
      inline bool Cell2DDoublePropertyExists(const string& propertyId) const
      {
        return _mesh.Cell2DDoublePropertyIndices.find(propertyId) !=
            _mesh.Cell2DDoublePropertyIndices.end();
      }
      inline unsigned int Cell2DDoublePropertyIndex(const string& propertyId) const
      {
        Output::Assert(Cell2DDoublePropertyExists(propertyId));
        return _mesh.Cell2DDoublePropertyIndices.at(propertyId);
      }
      inline unsigned int Cell2DDoublePropertySize(const unsigned int& cell2DIndex,
                                                   const unsigned int& propertyIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());
        return _mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex + 1] -
            _mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex];
      }
      inline double Cell2DDoublePropertyValue(const unsigned int& cell2DIndex,
                                              const unsigned int& propertyIndex,
                                              const unsigned int& propertyValueIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());
        Output::Assert(propertyValueIndex < Cell2DDoublePropertySize(cell2DIndex,
                                                                     propertyIndex));

        return _mesh.Cell2DDoublePropertyValues[propertyIndex][_mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex] +
            propertyValueIndex];
      }

      void Cell2DInitializeSubDivision(const unsigned int& cell2DIndex,
                                       const unsigned int& numberSubDivision);
      inline void Cell2DInsertSubDivision(const unsigned int& cell2DIndex,
                                          const unsigned int& subDivisionIndex,
                                          const unsigned int& cell0DIndex)
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Output::Assert(subDivisionIndex < Cell2DNumberSubDivision(cell2DIndex));
        Output::Assert(cell0DIndex < Cell0DTotalNumber());
        _mesh.Cell2DSubdivision[_mesh.NumberCell2DSubdivision[cell2DIndex] +
            subDivisionIndex] = cell0DIndex;
      }
      inline unsigned int Cell2DNumberSubDivision(const unsigned int& cell2DIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.NumberCell2DSubdivision[cell2DIndex + 1] -
            _mesh.NumberCell2DSubdivision[cell2DIndex];
      }
      inline unsigned int Cell2DSubDivisionCell0D(const unsigned int& cell2DIndex,
                                                  const unsigned int& subDivisionIndex) const
      {
        Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Output::Assert(subDivisionIndex < Cell2DNumberSubDivision(cell2DIndex));
        return _mesh.Cell2DSubdivision[_mesh.NumberCell2DSubdivision[cell2DIndex] + subDivisionIndex];
      }


      void Cell3DsInitialize(const unsigned int& numberCell3Ds);
      unsigned int Cell3DAppend(const unsigned int& numberCell3Ds);
      void Cell3DRemove(const unsigned int& cell3DIndex);
      inline void Cell3DSetMarker(const unsigned int& cell3DIndex,
                                  const unsigned int& marker)
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        _mesh.Cell3DMarkers[cell3DIndex] = marker;
      }
      inline void Cell3DSetState(const unsigned int& cell3DIndex,
                                 const bool& state)
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        _mesh.ActiveCell3D[cell3DIndex] = state;
      }
      void Cell3DInitializeVertices(const unsigned int& cell3DIndex,
                                    const unsigned int& numberCell3DVertices);
      void Cell3DInitializeEdges(const unsigned int& cell3DIndex,
                                 const unsigned int& numberCell3DEdges);
      void Cell3DInitializeFaces(const unsigned int& cell3DIndex,
                                 const unsigned int& numberCell3DFaces);
      inline void Cell3DInsertVertex(const unsigned int& cell3DIndex,
                                     const unsigned int& vertexIndex,
                                     const unsigned int& vertexCell0DIndex)
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Output::Assert(vertexIndex < Cell3DNumberVertices(cell3DIndex));
        Output::Assert(vertexCell0DIndex < Cell0DTotalNumber());

        _mesh.Cell3DVertices[_mesh.NumberCell3DVertices[cell3DIndex] +
            vertexIndex] = vertexCell0DIndex;
      }
      inline void Cell3DInsertEdge(const unsigned int& cell3DIndex,
                                   const unsigned int& edgeIndex,
                                   const unsigned int& edgeCell1DIndex)
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Output::Assert(edgeIndex < Cell3DNumberEdges(cell3DIndex));
        Output::Assert(edgeCell1DIndex < Cell1DTotalNumber());

        _mesh.Cell3DEdges[_mesh.NumberCell3DEdges[cell3DIndex] +
            edgeIndex] = edgeCell1DIndex;
      }
      inline void Cell3DInsertFace(const unsigned int& cell3DIndex,
                                   const unsigned int& faceIndex,
                                   const unsigned int& faceCell2DIndex)
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Output::Assert(faceIndex < Cell3DNumberFaces(cell3DIndex));
        Output::Assert(faceCell2DIndex < Cell2DTotalNumber());

        _mesh.Cell3DFaces[_mesh.NumberCell3DFaces[cell3DIndex] +
            faceIndex] = faceCell2DIndex;
      }
      void Cell3DAddVertices(const unsigned int& cell3DIndex,
                             const vector<unsigned int>& verticesCell0DIndices);
      void Cell3DAddEdges(const unsigned int& cell3DIndex,
                          const vector<unsigned int>& edgesCell0DIndices);
      void Cell3DAddFaces(const unsigned int& cell3DIndex,
                          const vector<unsigned int>& facesCell0DIndices);

      inline unsigned int Cell3DTotalNumber() const
      {
        return _mesh.NumberCell3D;
      }
      inline unsigned int Cell3DNumberVertices(const unsigned int& cell3DIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.NumberCell3DVertices[cell3DIndex + 1] -
            _mesh.NumberCell3DVertices[cell3DIndex];
      }
      inline unsigned int Cell3DNumberEdges(const unsigned int& cell3DIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.NumberCell3DEdges[cell3DIndex + 1] -
            _mesh.NumberCell3DEdges[cell3DIndex];
      }
      inline unsigned int Cell3DNumberFaces(const unsigned int& cell3DIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.NumberCell3DFaces[cell3DIndex + 1] -
            _mesh.NumberCell3DFaces[cell3DIndex];
      }
      inline unsigned int Cell3DVertex(const unsigned int& cell3DIndex,
                                       const unsigned int& vertexIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Output::Assert(vertexIndex < Cell3DNumberVertices(cell3DIndex));

        return _mesh.Cell3DVertices[_mesh.NumberCell3DVertices[cell3DIndex] + vertexIndex];
      }
      inline unsigned int Cell3DEdge(const unsigned int& cell3DIndex,
                                     const unsigned int& edgeIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Output::Assert(edgeIndex < Cell3DNumberEdges(cell3DIndex));

        return _mesh.Cell3DEdges[_mesh.NumberCell3DEdges[cell3DIndex] + edgeIndex];
      }
      inline unsigned int Cell3DFace(const unsigned int& cell3DIndex,
                                     const unsigned int& faceIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Output::Assert(faceIndex < Cell3DNumberFaces(cell3DIndex));

        return _mesh.Cell3DFaces[_mesh.NumberCell3DFaces[cell3DIndex] + faceIndex];
      }
      inline unsigned int Cell3DMarker(const unsigned int& cell3DIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.Cell3DMarkers[cell3DIndex];
      }
      inline bool Cell3DIsActive(const unsigned int& cell3DIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.ActiveCell3D[cell3DIndex];
      }

      inline bool Cell3DHasUpdatedCell3Ds(const unsigned int& cell3DIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.UpdatedCell3Ds.find(cell3DIndex) != _mesh.UpdatedCell3Ds.end();
      }
      inline unsigned int Cell3DNumberUpdatedCell3Ds(const unsigned int& cell3DIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.UpdatedCell3Ds.at(cell3DIndex).size();
      }
      inline bool Cell3DHasUpdatedCell3D(const unsigned int& cell3DIndex,
                                         const unsigned int& updatedCell3DIdex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Output::Assert(updatedCell3DIdex < Cell3DTotalNumber());
        return _mesh.UpdatedCell3Ds.at(cell3DIndex).find(updatedCell3DIdex) != _mesh.UpdatedCell3Ds.at(cell3DIndex).end();
      }
      void Cell3DInsertUpdatedCell3D(const unsigned int& cell3DIndex,
                                     const unsigned int& updatedCell3DIdex);
      bool Cell3DUpdatedCell3Ds(const unsigned int& cell3DIndex,
                                list<unsigned int>& updatedCell3DIds) const;

      inline void Cell3DSetId(const unsigned int& cell3DIndex,
                              const unsigned int& id) { return; }
      inline unsigned int Cell3DId(const unsigned int& cell3DIndex) const { return cell3DIndex; }

      void Cell3DInitializeDoubleProperties(const unsigned int& numberDoubleProperties);
      unsigned int Cell3DAddDoubleProperty(const string& propertyId);
      void Cell3DInitializeDoublePropertyValues(const unsigned int& cell3DIndex,
                                                const unsigned int& propertyIndex,
                                                const unsigned int& porpertySize);
      inline void Cell3DInsertDoublePropertyValue(const unsigned int& cell3DIndex,
                                                  const unsigned int& propertyIndex,
                                                  const unsigned int& propertyValueIndex,
                                                  const double& propertyValue)
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());

        _mesh.Cell3DDoublePropertyValues[propertyIndex][_mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex] +
            propertyValueIndex] = propertyValue;
      }
      inline unsigned int Cell3DNumberDoubleProperties() const
      {
        return _mesh.Cell3DDoublePropertyIds.size();
      }
      inline string Cell3DDoublePropertyId(const unsigned int& propertyIndex) const
      {
        return _mesh.Cell3DDoublePropertyIds[propertyIndex];
      }
      inline bool Cell3DDoublePropertyExists(const string& propertyId) const
      {
        return _mesh.Cell3DDoublePropertyIndices.find(propertyId) !=
            _mesh.Cell3DDoublePropertyIndices.end();
      }
      inline unsigned int Cell3DDoublePropertyIndex(const string& propertyId) const
      {
        Output::Assert(Cell3DDoublePropertyExists(propertyId));
        return _mesh.Cell3DDoublePropertyIndices.at(propertyId);
      }
      inline unsigned int Cell3DDoublePropertySize(const unsigned int& cell3DIndex,
                                                   const unsigned int& propertyIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());
        return _mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex + 1] -
            _mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex];
      }
      inline double Cell3DDoublePropertyValue(const unsigned int& cell3DIndex,
                                              const unsigned int& propertyIndex,
                                              const unsigned int& propertyValueIndex) const
      {
        Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());
        Output::Assert(propertyValueIndex < Cell3DDoublePropertySize(cell3DIndex,
                                                                     propertyIndex));

        return _mesh.Cell3DDoublePropertyValues[propertyIndex][_mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex] +
            propertyValueIndex];
      }

      void Compress();

      string ToString();

  };
}

#endif
