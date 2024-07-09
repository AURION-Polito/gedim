#ifndef __MeshMatricesWrapper_H
#define __MeshMatricesWrapper_H

#include "IOUtilities.hpp"
#include "Eigen/Eigen"
#include "MeshMatrices.hpp"
#include "IMeshDAO.hpp"

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
      void AlignMapContainerHigherElements(std::unordered_map<unsigned int, Container>& elements,
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
      void InitializeNuberVectorWithConstantElements(std::vector<unsigned int>& numberElementVector,
                                                     std::vector<T>& elementVector,
                                                     const unsigned int numberElementSize,
                                                     const unsigned int numberElements,
                                                     const T& elementInitialization = T());
      template<typename T>
      void InitializeNumberVector(std::vector<unsigned int>& numberElementVector,
                                  std::vector<T>& elementVector,
                                  const std::vector<unsigned int>& numberElements,
                                  const T& elementInitialization = T());

      template<typename T>
      void ResizeNumberVectorWithNewNumberElements(std::vector<unsigned int>& numberElementVector,
                                                   std::vector<T>& elementVector,
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
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        _mesh.Cell0DMarkers[cell0DIndex] = marker;
      }
      inline void Cell0DSetState(const unsigned int& cell0DIndex,
                                 const bool& state)
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        _mesh.ActiveCell0D[cell0DIndex] = state;
      }

      inline unsigned int Cell0DTotalNumber() const
      { return _mesh.NumberCell0D; }
      inline double Cell0DCoordinateX(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DCoordinates[3 * cell0DIndex];
      }
      inline double Cell0DCoordinateY(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DCoordinates[3 * cell0DIndex + 1];
      }
      inline double Cell0DCoordinateZ(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DCoordinates[3 * cell0DIndex + 2];
      }
      inline Eigen::Vector3d Cell0DCoordinates(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return Eigen::Vector3d(Cell0DCoordinateX(cell0DIndex),
                               Cell0DCoordinateY(cell0DIndex),
                               Cell0DCoordinateZ(cell0DIndex));
      }
      Eigen::MatrixXd Cell0DsCoordinates() const;
      inline unsigned int Cell0DMarker(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DMarkers[cell0DIndex];
      }
      inline std::vector<unsigned int> Cell0DsMarker() const
      { return _mesh.Cell0DMarkers; }
      inline bool Cell0DIsActive(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.ActiveCell0D[cell0DIndex];
      }

      inline bool Cell0DHasUpdatedCell0Ds(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.UpdatedCell0Ds.find(cell0DIndex) != _mesh.UpdatedCell0Ds.end();
      }
      inline unsigned int Cell0DNumberUpdatedCell0Ds(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.UpdatedCell0Ds.at(cell0DIndex).size();
      }
      inline bool Cell0DHasUpdatedCell0D(const unsigned int& cell0DIndex,
                                         const unsigned int& updatedCell0DIdex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(updatedCell0DIdex < Cell0DTotalNumber());
        return _mesh.UpdatedCell0Ds.at(cell0DIndex).find(updatedCell0DIdex) != _mesh.UpdatedCell0Ds.at(cell0DIndex).end();
      }
      void Cell0DInsertUpdatedCell0D(const unsigned int& cell0DIndex,
                                     const unsigned int& updatedCell0DIdex);
      bool Cell0DUpdatedCell0Ds(const unsigned int& cell0DIndex,
                                std::list<unsigned int>& updatedCell0DIds) const;

      std::vector<std::vector<unsigned int>> Cell0DsNeighbourCell1Ds() const;
      inline void Cell0DsInitializeNeighbourCell1Ds(const std::vector<unsigned int>& numberNeighbourCell1Ds);
      inline void Cell0DInitializeNeighbourCell1Ds(const unsigned int& cell0DIndex,
                                                   const unsigned int& numberNeighbourCell1Ds);
      inline void Cell0DInsertNeighbourCell1D(const unsigned int& cell0DIndex,
                                              const unsigned int& neighbourIndex,
                                              const unsigned int& neigbourCell1DIndex)
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell1D(cell0DIndex));
        Gedim::Output::Assert(neigbourCell1DIndex < Cell1DTotalNumber());

        _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] +
            neighbourIndex] = neigbourCell1DIndex;
      }
      inline unsigned int Cell0DNumberNeighbourCell1D(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.NumberCell0DNeighbourCell1D[cell0DIndex + 1] -
            _mesh.NumberCell0DNeighbourCell1D[cell0DIndex];
      }
      inline unsigned int Cell0DNeighbourCell1D(const unsigned int& cell0DIndex,
                                                const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell1D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] +
            neighbourIndex];
      }
      inline std::vector<unsigned int> Cell0DNeighbourCell1Ds(const unsigned int& cell0DIndex) const
      {
        const unsigned int numNeighs = Cell0DNumberNeighbourCell1D(cell0DIndex);
        std::vector<unsigned int> neighbours(numNeighs);
        for (unsigned int n = 0; n < numNeighs; n++)
          neighbours[n] = Cell0DNeighbourCell1D(cell0DIndex, n);

        return neighbours;
      }
      inline bool Cell0DHasNeighbourCell1D(const unsigned int& cell0DIndex,
                                           const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell1D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] +
            neighbourIndex] < _mesh.NumberCell1D;
      }
      inline void Cell0DResetNeighbourCell1D(const unsigned int& cell0DIndex,
                                             const unsigned int& neighbourIndex)
      {
        _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] +
            neighbourIndex] = std::numeric_limits<unsigned int>::max();
      }

      std::vector<std::vector<unsigned int>> Cell0DsNeighbourCell2Ds() const;
      inline void Cell0DsInitializeNeighbourCell2Ds(const std::vector<unsigned int>& numberNeighbourCell2Ds);
      inline void Cell0DInitializeNeighbourCell2Ds(const unsigned int& cell0DIndex,
                                                   const unsigned int& numberNeighbourCell2Ds);
      inline void Cell0DInsertNeighbourCell2D(const unsigned int& cell0DIndex,
                                              const unsigned int& neighbourIndex,
                                              const unsigned int& neigbourCell2DIndex)
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell2D(cell0DIndex));
        Gedim::Output::Assert(neigbourCell2DIndex < Cell2DTotalNumber());

        _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] +
            neighbourIndex] = neigbourCell2DIndex;
      }
      inline unsigned int Cell0DNumberNeighbourCell2D(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.NumberCell0DNeighbourCell2D[cell0DIndex + 1] -
            _mesh.NumberCell0DNeighbourCell2D[cell0DIndex];
      }
      inline unsigned int Cell0DNeighbourCell2D(const unsigned int& cell0DIndex,
                                                const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell2D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] +
            neighbourIndex];
      }
      inline std::vector<unsigned int> Cell0DNeighbourCell2Ds(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        const unsigned int numNeighs = Cell0DNumberNeighbourCell2D(cell0DIndex);
        std::vector<unsigned int> neighbours(numNeighs);
        for (unsigned int n = 0; n < numNeighs; n++)
          neighbours[n] = Cell0DNeighbourCell2D(cell0DIndex, n);
        return neighbours;
      }
      inline bool Cell0DHasNeighbourCell2D(const unsigned int& cell0DIndex,
                                           const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell2D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] +
            neighbourIndex] < _mesh.NumberCell2D;
      }
      inline void Cell0DResetNeighbourCell2D(const unsigned int& cell0DIndex,
                                             const unsigned int& neighbourIndex)
      {
        _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] +
            neighbourIndex] = std::numeric_limits<unsigned int>::max();
      }
      std::vector<std::vector<unsigned int>> Cell0DsNeighbourCell3Ds() const;
      void Cell0DsInitializeNeighbourCell3Ds(const std::vector<unsigned int>& numberNeighbourCell3Ds);
      void Cell0DInitializeNeighbourCell3Ds(const unsigned int& cell0DIndex,
                                            const unsigned int& numberNeighbourCell3Ds);
      void Cell0DInitializeNeighbourCell3Ds(const unsigned int& cell0DIndex,
                                            const std::vector<unsigned int>& neighbourCell3Ds);
      inline void Cell0DInsertNeighbourCell3D(const unsigned int& cell0DIndex,
                                              const unsigned int& neighbourIndex,
                                              const unsigned int& neigbourCell3DIndex)
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell3D(cell0DIndex));
        Gedim::Output::Assert(neigbourCell3DIndex < Cell3DTotalNumber());

        _mesh.Cell0DNeighbourCell3Ds[_mesh.NumberCell0DNeighbourCell3D[cell0DIndex] +
            neighbourIndex] = neigbourCell3DIndex;
      }
      inline unsigned int Cell0DNumberNeighbourCell3D(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.NumberCell0DNeighbourCell3D[cell0DIndex + 1] -
            _mesh.NumberCell0DNeighbourCell3D[cell0DIndex];
      }
      inline unsigned int Cell0DNeighbourCell3D(const unsigned int& cell0DIndex,
                                                const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell3D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell3Ds[_mesh.NumberCell0DNeighbourCell3D[cell0DIndex] +
            neighbourIndex];
      }
      inline std::vector<unsigned int> Cell0DNeighbourCell3Ds(const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        const unsigned int numNeighs = Cell0DNumberNeighbourCell3D(cell0DIndex);
        std::vector<unsigned int> neighbours(numNeighs);
        for (unsigned int n = 0; n < numNeighs; n++)
          neighbours[n] = Cell0DNeighbourCell3D(cell0DIndex, n);
        return neighbours;
      }
      inline bool Cell0DHasNeighbourCell3D(const unsigned int& cell0DIndex,
                                           const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell3D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell3Ds[_mesh.NumberCell0DNeighbourCell3D[cell0DIndex] +
            neighbourIndex] < _mesh.NumberCell3D;
      }
      inline void Cell0DResetNeighbourCell3D(const unsigned int& cell0DIndex,
                                             const unsigned int& neighbourIndex)
      {
        _mesh.Cell0DNeighbourCell3Ds[_mesh.NumberCell0DNeighbourCell3D[cell0DIndex] +
            neighbourIndex] = std::numeric_limits<unsigned int>::max();
      }

      void Cell0DInitializeDoubleProperties(const unsigned int& numberDoubleProperties);
      unsigned int Cell0DAddDoubleProperty(const std::string& propertyId);
      inline void Cell0DsInitializeDoublePropertyValues(const unsigned int& propertyIndex,
                                                        const std::vector<unsigned int>& propertySizes);
      void Cell0DInitializeDoublePropertyValues(const unsigned int& cell0DIndex,
                                                const unsigned int& propertyIndex,
                                                const unsigned int& propertySize);
      inline void Cell0DInsertDoublePropertyValue(const unsigned int& cell0DIndex,
                                                  const unsigned int& propertyIndex,
                                                  const unsigned int& propertyValueIndex,
                                                  const double& propertyValue)
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());

        _mesh.Cell0DDoublePropertyValues[propertyIndex][_mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex] +
            propertyValueIndex] = propertyValue;
      }
      inline unsigned int Cell0DNumberDoubleProperties() const
      {
        return _mesh.Cell0DDoublePropertyIds.size();
      }
      inline std::string Cell0DDoublePropertyId(const unsigned int& propertyIndex) const
      {
        return _mesh.Cell0DDoublePropertyIds[propertyIndex];
      }
      inline bool Cell0DDoublePropertyExists(const std::string& propertyId) const
      {
        return _mesh.Cell0DDoublePropertyIndices.find(propertyId) !=
            _mesh.Cell0DDoublePropertyIndices.end();
      }
      inline unsigned int Cell0DDoublePropertyIndex(const std::string& propertyId) const
      {
        Gedim::Output::Assert(Cell0DDoublePropertyExists(propertyId));
        return _mesh.Cell0DDoublePropertyIndices.at(propertyId);
      }
      inline unsigned int Cell0DDoublePropertySize(const unsigned int& cell0DIndex,
                                                   const unsigned int& propertyIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());
        return _mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex + 1] -
            _mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex];
      }
      inline double Cell0DDoublePropertyValue(const unsigned int& cell0DIndex,
                                              const unsigned int& propertyIndex,
                                              const unsigned int& propertyValueIndex) const
      {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());
        Gedim::Output::Assert(propertyValueIndex < Cell0DDoublePropertySize(cell0DIndex,
                                                                            propertyIndex));

        return _mesh.Cell0DDoublePropertyValues[propertyIndex][_mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex] +
            propertyValueIndex];
      }

      void Cell1DsInitialize(const unsigned int& numberCell1Ds);
      unsigned int Cell1DAppend(const unsigned int& numberCell1Ds);
      void Cell1DRemove(const unsigned int& cell1DIndex);
      void Cell1DInsertExtremes(const unsigned int& cell1DIndex,
                                const unsigned int& originCell0DIndex,
                                const unsigned int& endCell0DIndex);

      void Cell1DsInsertExtremes(const Eigen::MatrixXi& cell1DExtremes);

      Eigen::MatrixXi Cell1DsExtremes() const;
      /// \return the extrems as Eigen MatrixXi of cell1D, size 2
      /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
      inline Eigen::VectorXi Cell1DExtremes(const unsigned int& cell1DIndex) const
      { return (Eigen::VectorXi(2)<< Cell1DOrigin(cell1DIndex), Cell1DEnd(cell1DIndex)).finished(); }

      unsigned int Cell1DByExtremes(const unsigned int& originCell0DIndex,
                                    const unsigned int& endCell0DIndex) const;

      inline void Cell1DSetMarker(const unsigned int& cell1DIndex,
                                  const unsigned int& marker)
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        _mesh.Cell1DMarkers[cell1DIndex] = marker;
      }
      inline void Cell1DSetState(const unsigned int& cell1DIndex,
                                 const bool& state)
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        _mesh.ActiveCell1D[cell1DIndex] = state;
      }
      std::vector<std::vector<unsigned int>> Cell1DsNeighbourCell2Ds() const;
      inline void Cell1DsInitializeNeighbourCell2Ds(const std::vector<unsigned int>& numberNeighbourCell2Ds);
      inline void Cell1DsInitializeNeighbourCell2Ds(const unsigned int& numberNeighbourCell2Ds);
      void Cell1DInitializeNeighbourCell2Ds(const unsigned int& cell1DIndex,
                                            const unsigned int& numberNeighbourCell2Ds);
      inline void Cell1DInsertNeighbourCell2D(const unsigned int& cell1DIndex,
                                              const unsigned int& neighbourIndex,
                                              const unsigned int& neigbourCell2DIndex)
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell2D(cell1DIndex));
        Gedim::Output::Assert(neigbourCell2DIndex < Cell2DTotalNumber());

        _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] +
            neighbourIndex] = neigbourCell2DIndex;
      }
      inline unsigned int Cell1DTotalNumber() const
      { return _mesh.NumberCell1D; }
      inline unsigned int Cell1DVertex(const unsigned int& cell1DIndex,
                                       const unsigned int& vertexIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(vertexIndex < 2);
        return _mesh.Cell1DVertices[2 * cell1DIndex + vertexIndex];
      }
      inline unsigned int Cell1DOrigin(const unsigned int& cell1DIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DVertices[2 * cell1DIndex];
      }
      inline unsigned int Cell1DEnd(const unsigned int& cell1DIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DVertices[2 * cell1DIndex + 1];
      }
      inline unsigned int Cell1DFindExtreme(const unsigned int& cell1DIndex,
                                            const unsigned int& cell0DIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        if (_mesh.Cell1DVertices[2 * cell1DIndex] == cell0DIndex)
          return 0;
        else if (_mesh.Cell1DVertices[2 * cell1DIndex + 1] == cell0DIndex)
          return 1;
        else
          return 2;
      }
      inline Eigen::MatrixXd Cell1DCoordinates(const unsigned int& cell1DIndex) const
      {
        return (Eigen::MatrixXd(3, 2)<< Cell1DOriginCoordinates(cell1DIndex), Cell1DEndCoordinates(cell1DIndex)).finished();
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
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.NumberCell1DNeighbourCell2D[cell1DIndex + 1] -
            _mesh.NumberCell1DNeighbourCell2D[cell1DIndex];
      }
      inline unsigned int Cell1DNeighbourCell2D(const unsigned int& cell1DIndex,
                                                const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell2D(cell1DIndex));
        return _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] +
            neighbourIndex];
      }
      inline std::vector<unsigned int> Cell1DNeighbourCell2Ds(const unsigned int& cell1DIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        const unsigned int numNeighs = Cell1DNumberNeighbourCell2D(cell1DIndex);
        std::vector<unsigned int> neighbours(numNeighs);
        for (unsigned int n = 0; n < numNeighs; n++)
          neighbours[n] = Cell1DNeighbourCell2D(cell1DIndex, n);
        return neighbours;
      }
      inline bool Cell1DHasNeighbourCell2D(const unsigned int& cell1DIndex,
                                           const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell2D(cell1DIndex));
        return _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] +
            neighbourIndex] < _mesh.NumberCell2D;
      }
      inline void Cell1DResetNeighbourCell2D(const unsigned int& cell1DIndex,
                                             const unsigned int& neighbourIndex)
      {
        _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] +
            neighbourIndex] = std::numeric_limits<unsigned int>::max();
      }

      inline unsigned int Cell1DMarker(const unsigned int& cell1DIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DMarkers[cell1DIndex];
      }
      inline std::vector<unsigned int> Cell1DsMarker() const
      { return _mesh.Cell1DMarkers; }
      inline bool Cell1DIsActive(const unsigned int& cell1DIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.ActiveCell1D[cell1DIndex];
      }
      inline bool Cell1DHasOriginalCell1D(const unsigned int& updatedCell1DIndex) const
      {
        Gedim::Output::Assert(updatedCell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DOriginalCell1Ds.at(updatedCell1DIndex) < _mesh.NumberCell1D;
      }
      inline unsigned int Cell1DOriginalCell1D(const unsigned int& updatedCell1DIndex) const
      {
        Gedim::Output::Assert(updatedCell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DOriginalCell1Ds.at(updatedCell1DIndex);
      }
      inline bool Cell1DHasUpdatedCell1Ds(const unsigned int& cell1DIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.UpdatedCell1Ds.find(cell1DIndex) != _mesh.UpdatedCell1Ds.end();
      }
      inline unsigned int Cell1DNumberUpdatedCell1Ds(const unsigned int& cell1DIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.UpdatedCell1Ds.at(cell1DIndex).size();
      }
      inline bool Cell1DHasUpdatedCell1D(const unsigned int& cell1DIndex,
                                         const unsigned int& updatedCell1DIdex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(updatedCell1DIdex < Cell1DTotalNumber());
        return _mesh.UpdatedCell1Ds.at(cell1DIndex).find(updatedCell1DIdex) != _mesh.UpdatedCell1Ds.at(cell1DIndex).end();
      }
      void Cell1DInsertUpdatedCell1D(const unsigned int& cell1DIndex,
                                     const unsigned int& updatedCell1DIdex);
      bool Cell1DUpdatedCell1Ds(const unsigned int& cell1DIndex,
                                std::list<unsigned int>& updatedCell1DIds) const;
      void Cell1DInitializeDoubleProperties(const unsigned int& numberDoubleProperties);

      std::vector<std::vector<unsigned int>> Cell1DsNeighbourCell3Ds() const;
      void Cell1DsInitializeNeighbourCell3Ds(const std::vector<unsigned int>& numberNeighbourCell3Ds);
      void Cell1DInitializeNeighbourCell3Ds(const unsigned int& cell1DIndex,
                                            const unsigned int& numberNeighbourCell3Ds);
      inline void Cell1DInsertNeighbourCell3D(const unsigned int& cell1DIndex,
                                              const unsigned int& neighbourIndex,
                                              const unsigned int& neigbourCell3DIndex)
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell3D(cell1DIndex));
        Gedim::Output::Assert(neigbourCell3DIndex < Cell3DTotalNumber());

        _mesh.Cell1DNeighbourCell3Ds[_mesh.NumberCell1DNeighbourCell3D[cell1DIndex] +
            neighbourIndex] = neigbourCell3DIndex;
      }
      inline unsigned int Cell1DNumberNeighbourCell3D(const unsigned int& cell1DIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.NumberCell1DNeighbourCell3D[cell1DIndex + 1] -
            _mesh.NumberCell1DNeighbourCell3D[cell1DIndex];
      }
      inline unsigned int Cell1DNeighbourCell3D(const unsigned int& cell1DIndex,
                                                const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell3D(cell1DIndex));
        return _mesh.Cell1DNeighbourCell3Ds[_mesh.NumberCell1DNeighbourCell3D[cell1DIndex] +
            neighbourIndex];
      }
      inline bool Cell1DHasNeighbourCell3D(const unsigned int& cell1DIndex,
                                           const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell3D(cell1DIndex));
        return _mesh.Cell1DNeighbourCell3Ds[_mesh.NumberCell1DNeighbourCell3D[cell1DIndex] +
            neighbourIndex] < _mesh.NumberCell3D;
      }
      inline void Cell1DResetNeighbourCell3D(const unsigned int& cell1DIndex,
                                             const unsigned int& neighbourIndex)
      {
        _mesh.Cell1DNeighbourCell3Ds[_mesh.NumberCell1DNeighbourCell3D[cell1DIndex] +
            neighbourIndex] = std::numeric_limits<unsigned int>::max();
      }

      unsigned int Cell1DAddDoubleProperty(const std::string& propertyId);
      inline void Cell1DsInitializeDoublePropertyValues(const unsigned int& propertyIndex,
                                                        const std::vector<unsigned int>& propertySizes);
      void Cell1DInitializeDoublePropertyValues(const unsigned int& cell1DIndex,
                                                const unsigned int& propertyIndex,
                                                const unsigned int& propertySize);
      inline void Cell1DInsertDoublePropertyValue(const unsigned int& cell1DIndex,
                                                  const unsigned int& propertyIndex,
                                                  const unsigned int& propertyValueIndex,
                                                  const double& propertyValue)
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());

        _mesh.Cell1DDoublePropertyValues[propertyIndex][_mesh.Cell1DDoublePropertySizes[propertyIndex][cell1DIndex] +
            propertyValueIndex] = propertyValue;
      }
      inline unsigned int Cell1DNumberDoubleProperties() const
      {
        return _mesh.Cell1DDoublePropertyIds.size();
      }
      inline std::string Cell1DDoublePropertyId(const unsigned int& propertyIndex) const
      {
        return _mesh.Cell1DDoublePropertyIds[propertyIndex];
      }
      inline bool Cell1DDoublePropertyExists(const std::string& propertyId) const
      {
        return _mesh.Cell1DDoublePropertyIndices.find(propertyId) !=
            _mesh.Cell1DDoublePropertyIndices.end();
      }
      inline unsigned int Cell1DDoublePropertyIndex(const std::string& propertyId) const
      {
        Gedim::Output::Assert(Cell1DDoublePropertyExists(propertyId));
        return _mesh.Cell1DDoublePropertyIndices.at(propertyId);
      }
      inline unsigned int Cell1DDoublePropertySize(const unsigned int& cell1DIndex,
                                                   const unsigned int& propertyIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());
        return _mesh.Cell1DDoublePropertySizes[propertyIndex][cell1DIndex + 1] -
            _mesh.Cell1DDoublePropertySizes[propertyIndex][cell1DIndex];
      }
      inline double Cell1DDoublePropertyValue(const unsigned int& cell1DIndex,
                                              const unsigned int& propertyIndex,
                                              const unsigned int& propertyValueIndex) const
      {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());
        Gedim::Output::Assert(propertyValueIndex < Cell1DDoublePropertySize(cell1DIndex,
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
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        _mesh.Cell2DMarkers[cell2DIndex] = marker;
      }
      inline void Cell2DSetState(const unsigned int& cell2DIndex,
                                 const bool& state)
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        _mesh.ActiveCell2D[cell2DIndex] = state;
      }
      inline void Cell2DsInitializeVertices(const unsigned int& numberCell2DVertices);
      inline void Cell2DsInitializeVertices(const std::vector<unsigned int>& numberCell2DsVertices);
      void Cell2DInitializeVertices(const unsigned int& cell2DIndex,
                                    const unsigned int& numberCell2DVertices);
      inline void Cell2DsInitializeEdges(const unsigned int& numberCell2DEdges);
      inline void Cell2DsInitializeEdges(const std::vector<unsigned int>& numberCell2DsEdges);
      void Cell2DInitializeEdges(const unsigned int& cell2DIndex,
                                 const unsigned int& numberCell2DEdges);
      inline void Cell2DInsertVertices(const unsigned int& cell2DIndex,
                                       const std::vector<unsigned int>& verticesCell0DIndices)
      {
        for (unsigned int v = 0; v < verticesCell0DIndices.size(); v++)
          Cell2DInsertVertex(cell2DIndex, v, verticesCell0DIndices[v]);
      }
      inline void Cell2DInsertVertex(const unsigned int& cell2DIndex,
                                     const unsigned int& vertexIndex,
                                     const unsigned int& vertexCell0DIndex)
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell2DNumberVertices(cell2DIndex));
        Gedim::Output::Assert(vertexCell0DIndex < Cell0DTotalNumber());
        _mesh.Cell2DVertices[_mesh.NumberCell2DVertices[cell2DIndex] +
            vertexIndex] = vertexCell0DIndex;
      }
      inline void Cell2DInsertEdges(const unsigned int& cell2DIndex,
                                    const std::vector<unsigned int>& edgesCell1DIndices)
      {
        for (unsigned int e = 0; e < edgesCell1DIndices.size(); e++)
          Cell2DInsertEdge(cell2DIndex, e, edgesCell1DIndices[e]);
      }
      inline void Cell2DInsertEdge(const unsigned int& cell2DIndex,
                                   const unsigned int& edgeIndex,
                                   const unsigned int& edgeCell1DIndex)
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(edgeIndex < Cell2DNumberEdges(cell2DIndex));
        Gedim::Output::Assert(edgeCell1DIndex < Cell1DTotalNumber());
        _mesh.Cell2DEdges[_mesh.NumberCell2DEdges[cell2DIndex] +
            edgeIndex] = edgeCell1DIndex;
      }
      void Cell2DAddVertices(const unsigned int& cell2DIndex,
                             const std::vector<unsigned int>& verticesCell0DIndices);
      void Cell2DAddEdges(const unsigned int& cell2DIndex,
                          const std::vector<unsigned int>& edgesCell1DIndices);
      void Cell2DAddVerticesAndEdges(const unsigned int& cell2DIndex,
                                     const Eigen::MatrixXi& verticesAndEdgesIndices);

      inline unsigned int Cell2DTotalNumber() const
      { return _mesh.NumberCell2D; }
      inline unsigned int Cell2DNumberVertices(const unsigned int& cell2DIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.NumberCell2DVertices[cell2DIndex + 1] -
            _mesh.NumberCell2DVertices[cell2DIndex];
      }
      inline unsigned int Cell2DNumberEdges(const unsigned int& cell2DIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.NumberCell2DEdges[cell2DIndex + 1] - _mesh.NumberCell2DEdges[cell2DIndex];
      }

      inline std::vector<unsigned int> Cell2DVertices(const unsigned int& cell2DIndex) const
      {
        return std::vector<unsigned int>(_mesh.Cell2DVertices.begin() + _mesh.NumberCell2DVertices[cell2DIndex],
                                         _mesh.Cell2DVertices.begin() + _mesh.NumberCell2DVertices[cell2DIndex] + Cell2DNumberVertices(cell2DIndex));
      }
      std::vector<std::vector<unsigned int>> Cell2DsVertices() const;
      std::vector<Eigen::MatrixXi> Cell2DsExtremes() const;

      inline unsigned int Cell2DVertex(const unsigned int& cell2DIndex,
                                       const unsigned int& vertexIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell2DNumberVertices(cell2DIndex));
        return _mesh.Cell2DVertices[_mesh.NumberCell2DVertices[cell2DIndex] + vertexIndex];
      }
      inline Eigen::Vector3d Cell2DVertexCoordinates(const unsigned int& cell2DIndex,
                                                     const unsigned int& vertexIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell2DNumberVertices(cell2DIndex));
        return Cell0DCoordinates(Cell2DVertex(cell2DIndex, vertexIndex));
      }
      Eigen::MatrixXd Cell2DVerticesCoordinates(const unsigned int& cell2DIndex) const;
      unsigned int Cell2DFindVertex(const unsigned int& cell2DIndex,
                                    const unsigned int& cell0DIndex) const;

      inline std::vector<unsigned int> Cell2DEdges(const unsigned int& cell2DIndex) const
      {
        return std::vector<unsigned int>(_mesh.Cell2DEdges.begin() + _mesh.NumberCell2DEdges[cell2DIndex],
                                         _mesh.Cell2DEdges.begin() + _mesh.NumberCell2DEdges[cell2DIndex] + Cell2DNumberEdges(cell2DIndex));
      }

      inline unsigned int Cell2DEdge(const unsigned int& cell2DIndex,
                                     const unsigned int& edgeIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(edgeIndex < Cell2DNumberEdges(cell2DIndex));
        return _mesh.Cell2DEdges[_mesh.NumberCell2DEdges[cell2DIndex] + edgeIndex];
      }
      unsigned int Cell2DFindEdge(const unsigned int& cell2DIndex,
                                  const unsigned int& cell1DIndex) const;
      unsigned int Cell2DFindEdgeByExtremes(const unsigned int& cell2DIndex,
                                            const unsigned int& originCell0DIndex,
                                            const unsigned int& endCell0DIndex) const;
      inline unsigned int Cell2DMarker(const unsigned int& cell2DIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.Cell2DMarkers[cell2DIndex];
      }
      inline std::vector<unsigned int> Cell2DsMarker() const
      { return _mesh.Cell2DMarkers; }
      inline bool Cell2DIsActive(const unsigned int& cell2DIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.ActiveCell2D[cell2DIndex];
      }

      inline bool Cell2DHasUpdatedCell2Ds(const unsigned int& cell2DIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.UpdatedCell2Ds.find(cell2DIndex) != _mesh.UpdatedCell2Ds.end();
      }
      inline unsigned int Cell2DNumberUpdatedCell2Ds(const unsigned int& cell2DIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.UpdatedCell2Ds.at(cell2DIndex).size();
      }
      inline bool Cell2DHasUpdatedCell2D(const unsigned int& cell2DIndex,
                                         const unsigned int& updatedCell2DIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(updatedCell2DIndex < Cell2DTotalNumber());
        return _mesh.UpdatedCell2Ds.at(cell2DIndex).find(updatedCell2DIndex) != _mesh.UpdatedCell2Ds.at(cell2DIndex).end();
      }
      void Cell2DInsertUpdatedCell2D(const unsigned int& cell2DIndex,
                                     const unsigned int& updatedCell2DIdex);
      inline bool Cell2DHasOriginalCell2D(const unsigned int& updatedCell2DIndex) const
      {
        Gedim::Output::Assert(updatedCell2DIndex < Cell2DTotalNumber());
        return _mesh.Cell2DOriginalCell2Ds.at(updatedCell2DIndex) < _mesh.NumberCell2D;
      }
      inline unsigned int Cell2DOriginalCell2D(const unsigned int& updatedCell2DIndex) const
      {
        Gedim::Output::Assert(updatedCell2DIndex < Cell2DTotalNumber());
        return _mesh.Cell2DOriginalCell2Ds.at(updatedCell2DIndex);
      }
      bool Cell2DUpdatedCell2Ds(const unsigned int& cell2DIndex,
                                std::list<unsigned int>& updatedCell2DIds) const;

      std::vector<std::vector<unsigned int>> Cell2DsNeighbourCell3Ds() const;
      inline void Cell2DsInitializeNeighbourCell3Ds(const std::vector<unsigned int>& numberNeighbourCell3Ds);
      void Cell2DInitializeNeighbourCell3Ds(const unsigned int& cell2DIndex,
                                            const unsigned int& numberNeighbourCell3Ds);
      inline void Cell2DInsertNeighbourCell3D(const unsigned int& cell2DIndex,
                                              const unsigned int& neighbourIndex,
                                              const unsigned int& neigbourCell3DIndex)
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell2DNumberNeighbourCell3D(cell2DIndex));
        Gedim::Output::Assert(neigbourCell3DIndex < Cell3DTotalNumber());

        _mesh.Cell2DNeighbourCell3Ds[_mesh.NumberCell2DNeighbourCell3D[cell2DIndex] +
            neighbourIndex] = neigbourCell3DIndex;
      }
      inline unsigned int Cell2DNumberNeighbourCell3D(const unsigned int& cell2DIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.NumberCell2DNeighbourCell3D[cell2DIndex + 1] -
            _mesh.NumberCell2DNeighbourCell3D[cell2DIndex];
      }
      inline unsigned int Cell2DNeighbourCell3D(const unsigned int& cell2DIndex,
                                                const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell2DNumberNeighbourCell3D(cell2DIndex));
        return _mesh.Cell2DNeighbourCell3Ds[_mesh.NumberCell2DNeighbourCell3D[cell2DIndex] +
            neighbourIndex];
      }
      inline bool Cell2DHasNeighbourCell3D(const unsigned int& cell2DIndex,
                                           const unsigned int& neighbourIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell2DNumberNeighbourCell3D(cell2DIndex));
        return _mesh.Cell2DNeighbourCell3Ds[_mesh.NumberCell2DNeighbourCell3D[cell2DIndex] +
            neighbourIndex] < _mesh.NumberCell3D;
      }
      inline void Cell2DResetNeighbourCell3D(const unsigned int& cell2DIndex,
                                             const unsigned int& neighbourIndex)
      {
        _mesh.Cell2DNeighbourCell3Ds[_mesh.NumberCell2DNeighbourCell3D[cell2DIndex] +
            neighbourIndex] = std::numeric_limits<unsigned int>::max();
      }
      void Cell2DInitializeDoubleProperties(const unsigned int& numberDoubleProperties);
      unsigned int Cell2DAddDoubleProperty(const std::string& propertyId);
      inline void Cell2DsInitializeDoublePropertyValues(const unsigned int& propertyIndex,
                                                        const std::vector<unsigned int>& propertySizes);
      void Cell2DInitializeDoublePropertyValues(const unsigned int& cell2DIndex,
                                                const unsigned int& propertyIndex,
                                                const unsigned int& propertySize);
      inline void Cell2DInsertDoublePropertyValue(const unsigned int& cell2DIndex,
                                                  const unsigned int& propertyIndex,
                                                  const unsigned int& propertyValueIndex,
                                                  const double& propertyValue)
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());

        _mesh.Cell2DDoublePropertyValues[propertyIndex][_mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex] +
            propertyValueIndex] = propertyValue;
      }
      inline unsigned int Cell2DNumberDoubleProperties() const
      {
        return _mesh.Cell2DDoublePropertyIds.size();
      }
      inline std::string Cell2DDoublePropertyId(const unsigned int& propertyIndex) const
      {
        return _mesh.Cell2DDoublePropertyIds[propertyIndex];
      }
      inline bool Cell2DDoublePropertyExists(const std::string& propertyId) const
      {
        return _mesh.Cell2DDoublePropertyIndices.find(propertyId) !=
            _mesh.Cell2DDoublePropertyIndices.end();
      }
      inline unsigned int Cell2DDoublePropertyIndex(const std::string& propertyId) const
      {
        Gedim::Output::Assert(Cell2DDoublePropertyExists(propertyId));
        return _mesh.Cell2DDoublePropertyIndices.at(propertyId);
      }
      inline unsigned int Cell2DDoublePropertySize(const unsigned int& cell2DIndex,
                                                   const unsigned int& propertyIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());
        return _mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex + 1] -
            _mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex];
      }
      inline double Cell2DDoublePropertyValue(const unsigned int& cell2DIndex,
                                              const unsigned int& propertyIndex,
                                              const unsigned int& propertyValueIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());
        Gedim::Output::Assert(propertyValueIndex < Cell2DDoublePropertySize(cell2DIndex,
                                                                            propertyIndex));

        return _mesh.Cell2DDoublePropertyValues[propertyIndex][_mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex] +
            propertyValueIndex];
      }

      inline void Cell2DsInitializeSubDivision(const std::vector<unsigned int>& numberSubDivisions);

      void Cell2DInitializeSubDivision(const unsigned int& cell2DIndex,
                                       const unsigned int& numberSubDivision);
      inline void Cell2DInsertSubDivision(const unsigned int& cell2DIndex,
                                          const unsigned int& subDivisionIndex,
                                          const unsigned int& cell0DIndex)
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(subDivisionIndex < Cell2DNumberSubDivision(cell2DIndex));
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        _mesh.Cell2DSubdivision[_mesh.NumberCell2DSubdivision[cell2DIndex] +
            subDivisionIndex] = cell0DIndex;
      }
      inline unsigned int Cell2DNumberSubDivision(const unsigned int& cell2DIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.NumberCell2DSubdivision[cell2DIndex + 1] -
            _mesh.NumberCell2DSubdivision[cell2DIndex];
      }
      inline unsigned int Cell2DSubDivisionCell0D(const unsigned int& cell2DIndex,
                                                  const unsigned int& subDivisionIndex) const
      {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(subDivisionIndex < Cell2DNumberSubDivision(cell2DIndex));
        return _mesh.Cell2DSubdivision[_mesh.NumberCell2DSubdivision[cell2DIndex] + subDivisionIndex];
      }


      void Cell3DsInitialize(const unsigned int& numberCell3Ds);
      unsigned int Cell3DAppend(const unsigned int& numberCell3Ds);
      void Cell3DRemove(const unsigned int& cell3DIndex);
      inline void Cell3DSetMarker(const unsigned int& cell3DIndex,
                                  const unsigned int& marker)
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        _mesh.Cell3DMarkers[cell3DIndex] = marker;
      }
      inline void Cell3DSetState(const unsigned int& cell3DIndex,
                                 const bool& state)
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        _mesh.ActiveCell3D[cell3DIndex] = state;
      }
      inline void Cell3DsInitializeVertices(const std::vector<unsigned int>& numberCell3DsVertices);
      inline void Cell3DsInitializeEdges(const std::vector<unsigned int>& numberCell3DsEdges);
      inline void Cell3DsInitializeFaces(const std::vector<unsigned int>& numberCell3DsFaces);
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
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell3DNumberVertices(cell3DIndex));
        Gedim::Output::Assert(vertexCell0DIndex < Cell0DTotalNumber());

        _mesh.Cell3DVertices[_mesh.NumberCell3DVertices[cell3DIndex] +
            vertexIndex] = vertexCell0DIndex;
      }
      inline void Cell3DInsertEdge(const unsigned int& cell3DIndex,
                                   const unsigned int& edgeIndex,
                                   const unsigned int& edgeCell1DIndex)
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(edgeIndex < Cell3DNumberEdges(cell3DIndex));
        Gedim::Output::Assert(edgeCell1DIndex < Cell1DTotalNumber());

        _mesh.Cell3DEdges[_mesh.NumberCell3DEdges[cell3DIndex] +
            edgeIndex] = edgeCell1DIndex;
      }
      inline void Cell3DInsertFace(const unsigned int& cell3DIndex,
                                   const unsigned int& faceIndex,
                                   const unsigned int& faceCell2DIndex)
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(faceIndex < Cell3DNumberFaces(cell3DIndex));
        Gedim::Output::Assert(faceCell2DIndex < Cell2DTotalNumber());

        _mesh.Cell3DFaces[_mesh.NumberCell3DFaces[cell3DIndex] +
            faceIndex] = faceCell2DIndex;
      }
      void Cell3DAddVertices(const unsigned int& cell3DIndex,
                             const std::vector<unsigned int>& verticesCell0DIndices);
      void Cell3DAddEdges(const unsigned int& cell3DIndex,
                          const std::vector<unsigned int>& edgesCell1DIndices);
      void Cell3DAddFaces(const unsigned int& cell3DIndex,
                          const std::vector<unsigned int>& facesCell2DIndices);

      unsigned int Cell3DFindVertex(const unsigned int& cell3DIndex,
                                    const unsigned int& cell0DIndex) const;
      unsigned int Cell3DFindEdge(const unsigned int& cell3DIndex,
                                  const unsigned int& cell1DIndex) const;
      unsigned int Cell3DFindFace(const unsigned int& cell3DIndex,
                                  const unsigned int& cell2DIndex) const;

      unsigned int Cell3DFindEdgeByExtremes(const unsigned int& cell3DIndex,
                                            const unsigned int& originCell0DIndex,
                                            const unsigned int& endCell0DIndex) const;

      inline unsigned int Cell3DTotalNumber() const
      {
        return _mesh.NumberCell3D;
      }
      inline unsigned int Cell3DNumberVertices(const unsigned int& cell3DIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.NumberCell3DVertices[cell3DIndex + 1] -
            _mesh.NumberCell3DVertices[cell3DIndex];
      }
      inline unsigned int Cell3DNumberEdges(const unsigned int& cell3DIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.NumberCell3DEdges[cell3DIndex + 1] -
            _mesh.NumberCell3DEdges[cell3DIndex];
      }
      inline unsigned int Cell3DNumberFaces(const unsigned int& cell3DIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.NumberCell3DFaces[cell3DIndex + 1] -
            _mesh.NumberCell3DFaces[cell3DIndex];
      }
      std::vector<unsigned int> Cell3DVertices(const unsigned int& cell3DIndex) const;
      inline unsigned int Cell3DVertex(const unsigned int& cell3DIndex,
                                       const unsigned int& vertexIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell3DNumberVertices(cell3DIndex));

        return _mesh.Cell3DVertices[_mesh.NumberCell3DVertices[cell3DIndex] + vertexIndex];
      }
      inline Eigen::Vector3d Cell3DVertexCoordinates(const unsigned int& cell3DIndex,
                                                     const unsigned int& vertexIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell3DNumberVertices(cell3DIndex));
        return Cell0DCoordinates(Cell3DVertex(cell3DIndex, vertexIndex));
      }
      Eigen::MatrixXd Cell3DVerticesCoordinates(const unsigned int& cell3DIndex) const;
      std::vector<unsigned int> Cell3DEdges(const unsigned int& cell3DIndex) const;
      inline unsigned int Cell3DEdge(const unsigned int& cell3DIndex,
                                     const unsigned int& edgeIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(edgeIndex < Cell3DNumberEdges(cell3DIndex));

        return _mesh.Cell3DEdges[_mesh.NumberCell3DEdges[cell3DIndex] + edgeIndex];
      }
      std::vector<unsigned int> Cell3DFaces(const unsigned int& cell3DIndex) const;
      inline unsigned int Cell3DFace(const unsigned int& cell3DIndex,
                                     const unsigned int& faceIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(faceIndex < Cell3DNumberFaces(cell3DIndex));

        return _mesh.Cell3DFaces[_mesh.NumberCell3DFaces[cell3DIndex] + faceIndex];
      }
      std::vector<std::vector<std::vector<unsigned int>>> Cell3DsFacesVertices() const;
      std::vector<std::vector<unsigned int>> Cell3DsVertices() const;
      std::vector<std::vector<unsigned int>> Cell3DsEdges() const;
      std::vector<std::vector<unsigned int>> Cell3DsFaces() const;
      inline unsigned int Cell3DMarker(const unsigned int& cell3DIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.Cell3DMarkers[cell3DIndex];
      }
      inline std::vector<unsigned int> Cell3DsMarker() const
      { return _mesh.Cell3DMarkers; }
      inline bool Cell3DIsActive(const unsigned int& cell3DIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.ActiveCell3D[cell3DIndex];
      }

      inline bool Cell3DHasUpdatedCell3Ds(const unsigned int& cell3DIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.UpdatedCell3Ds.find(cell3DIndex) != _mesh.UpdatedCell3Ds.end();
      }
      inline unsigned int Cell3DNumberUpdatedCell3Ds(const unsigned int& cell3DIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.UpdatedCell3Ds.at(cell3DIndex).size();
      }
      inline bool Cell3DHasUpdatedCell3D(const unsigned int& cell3DIndex,
                                         const unsigned int& updatedCell3DIdex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(updatedCell3DIdex < Cell3DTotalNumber());
        return _mesh.UpdatedCell3Ds.at(cell3DIndex).find(updatedCell3DIdex) != _mesh.UpdatedCell3Ds.at(cell3DIndex).end();
      }
      void Cell3DInsertUpdatedCell3D(const unsigned int& cell3DIndex,
                                     const unsigned int& updatedCell3DIdex);
      bool Cell3DUpdatedCell3Ds(const unsigned int& cell3DIndex,
                                std::list<unsigned int>& updatedCell3DIds) const;
      inline bool Cell3DHasOriginalCell3D(const unsigned int& updatedCell3DIndex) const
      {
        Gedim::Output::Assert(updatedCell3DIndex < Cell3DTotalNumber());
        return _mesh.Cell3DOriginalCell3Ds.at(updatedCell3DIndex) < _mesh.NumberCell3D;
      }
      inline unsigned int Cell3DOriginalCell3D(const unsigned int& updatedCell3DIndex) const
      {
        Gedim::Output::Assert(updatedCell3DIndex < Cell3DTotalNumber());
        return _mesh.Cell3DOriginalCell3Ds.at(updatedCell3DIndex);
      }

      void Cell3DInitializeDoubleProperties(const unsigned int& numberDoubleProperties);
      unsigned int Cell3DAddDoubleProperty(const std::string& propertyId);
      inline void Cell3DsInitializeDoublePropertyValues(const unsigned int& propertyIndex,
                                                        const std::vector<unsigned int>& propertySizes);
      void Cell3DInitializeDoublePropertyValues(const unsigned int& cell3DIndex,
                                                const unsigned int& propertyIndex,
                                                const unsigned int& propertySize);
      inline void Cell3DInsertDoublePropertyValue(const unsigned int& cell3DIndex,
                                                  const unsigned int& propertyIndex,
                                                  const unsigned int& propertyValueIndex,
                                                  const double& propertyValue)
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());

        _mesh.Cell3DDoublePropertyValues[propertyIndex][_mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex] +
            propertyValueIndex] = propertyValue;
      }
      inline unsigned int Cell3DNumberDoubleProperties() const
      {
        return _mesh.Cell3DDoublePropertyIds.size();
      }
      inline std::string Cell3DDoublePropertyId(const unsigned int& propertyIndex) const
      {
        return _mesh.Cell3DDoublePropertyIds[propertyIndex];
      }
      inline bool Cell3DDoublePropertyExists(const std::string& propertyId) const
      {
        return _mesh.Cell3DDoublePropertyIndices.find(propertyId) !=
            _mesh.Cell3DDoublePropertyIndices.end();
      }
      inline unsigned int Cell3DDoublePropertyIndex(const std::string& propertyId) const
      {
        Gedim::Output::Assert(Cell3DDoublePropertyExists(propertyId));
        return _mesh.Cell3DDoublePropertyIndices.at(propertyId);
      }
      inline unsigned int Cell3DDoublePropertySize(const unsigned int& cell3DIndex,
                                                   const unsigned int& propertyIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());
        return _mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex + 1] -
            _mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex];
      }
      inline double Cell3DDoublePropertyValue(const unsigned int& cell3DIndex,
                                              const unsigned int& propertyIndex,
                                              const unsigned int& propertyValueIndex) const
      {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());
        Gedim::Output::Assert(propertyValueIndex < Cell3DDoublePropertySize(cell3DIndex,
                                                                            propertyIndex));

        return _mesh.Cell3DDoublePropertyValues[propertyIndex][_mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex] +
            propertyValueIndex];
      }

      void Compress();

      std::string ToString();

  };
}

#endif
