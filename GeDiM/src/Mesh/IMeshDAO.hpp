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

#ifndef __IMeshWrapper_H
#define __IMeshWrapper_H

#include "Eigen/Eigen"
#include "Gedim_Macro.hpp"
#include "IOUtilities.hpp"
#include "MeshMatrices.hpp"

namespace Gedim
{
/// \brief The IMeshDAO (mesh data access object) class to read and write mesh data
class IMeshDAO
{
  public:
    virtual ~IMeshDAO()
    {
    }

    /// \brief Initialize the mesh dimension
    virtual void InitializeDimension(const unsigned int &dimension) = 0;
    /// \return the geometric dimension of the mesh
    virtual unsigned int Dimension() const = 0;

    /// \brief Initialize the Cell0Ds container
    /// \param numberCell0Ds the total number of Cell0Ds
    /// \note No reset of Cell0Ds is performed
    virtual void Cell0DsInitialize(const unsigned int &numberCell0Ds) = 0;
    /// \brief Append Cell0Ds to the Cell0Ds container
    /// \param numberCell0Ds the number of Cell0Ds to append
    /// \return the previous number of Cell0Ds before the append operation
    virtual unsigned int Cell0DAppend(const unsigned int &numberCell0Ds) = 0;

    /// \brief Remove the Cell0D from the mesh
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \note the cell0D is removed and no integrity check in the mesh are performed
    virtual void Cell0DRemove(const unsigned int &cell0DIndex) = 0;
    /// \brief Add the Cell0D Coordinates
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param coordinates the coordinates of the Cell0D
    virtual void Cell0DInsertCoordinates(const unsigned int &cell0DIndex, const Eigen::Vector3d &coordinates) = 0;
    /// \brief Add the Cell0Ds Coordinates
    /// \param coordinates the coordinates of the Cell0Ds, size 3 x Cell0DTotalNumber()
    virtual void Cell0DsInsertCoordinates(const Eigen::MatrixXd &coordinates) = 0;
    /// \brief Set the Cell0D Marker
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param marker the marker of the Cell0D
    virtual void Cell0DSetMarker(const unsigned int &cell0DIndex, const unsigned int &marker) = 0;
    /// \brief Set the Cell0D state
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param state true if Cell0D is active, false otherwise
    virtual void Cell0DSetState(const unsigned int &cell0DIndex, const bool &state) = 0;
    /// \return the total number of Cell0Ds
    virtual unsigned int Cell0DTotalNumber() const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the X coordinate of cell0D
    virtual double Cell0DCoordinateX(const unsigned int &cell0DIndex) const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the Y coordinate of cell0D
    virtual double Cell0DCoordinateY(const unsigned int &cell0DIndex) const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the Z coordinate of cell0D
    virtual double Cell0DCoordinateZ(const unsigned int &cell0DIndex) const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the coordinates as Eigen Vector3d of cell0D, size 3x1
    virtual Eigen::Vector3d Cell0DCoordinates(const unsigned int &cell0DIndex) const = 0;
    /// \return the coordinates as Eigen MatrixXd of cell0D, size 3xCell0DTotalNumber()
    virtual Eigen::MatrixXd Cell0DsCoordinates() const = 0;
    virtual Eigen::MatrixXd Cell0DsCoordinates(const std::vector<unsigned int> &cell0Ds) const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return if the cell0D is active
    virtual bool Cell0DIsActive(const unsigned int &cell0DIndex) const = 0;
    /// \return the activation state of all cell0Ds
    virtual std::vector<bool> Cell0DsState() const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the cell0D marker
    virtual unsigned int Cell0DMarker(const unsigned int &cell0DIndex) const = 0;
    virtual std::vector<unsigned int> Cell0DsMarker() const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return if the cell0D has new cell0Ds associated
    virtual bool Cell0DHasUpdatedCell0Ds(const unsigned int &cell0DIndex) const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the number of new cell0Ds associated to cell0DIndex
    virtual unsigned int Cell0DNumberUpdatedCell0Ds(const unsigned int &cell0DIndex) const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \param updatedCell0DIdex the index of the new cell0D from 0 to Cell0DTotalNumber()
    /// \return if the cell0D has the updatedCell0DIdex associated
    virtual bool Cell0DHasUpdatedCell0D(const unsigned int &cell0DIndex, const unsigned int &updatedCell0DIdex) const = 0;
    /// \brief Add the new Cell0D to an existing Cell0D
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param updatedCell0DIdex the index of the new cell0D from 0 to Cell0DTotalNumber()
    virtual void Cell0DInsertUpdatedCell0D(const unsigned int &cell0DIndex, const unsigned int &updatedCell0DIdex) = 0;
    /// \brief return the updated Cell0D Ids for cell0DIndex
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param updatedCell0DIds the list of the new Cell0D Ids associated to cell0DIndex
    /// \return true if the cell0DIndex is contained in the updatedCell0DIds list, false otherwise
    virtual bool Cell0DUpdatedCell0Ds(const unsigned int &cell0DIndex, std::list<unsigned int> &updatedCell0DIds) const = 0;

    virtual std::vector<std::vector<unsigned int>> Cell0DsNeighbourCell1Ds() const = 0;
    /// \brief Initialize the Cell0Ds Cell1D neighbours number
    /// \param numbersNeighbourCell1Ds the number of Cell1D neighbours of each Cell0D, size 1 x Cell0DTotalNumber()
    virtual void Cell0DsInitializeNeighbourCell1Ds(const std::vector<unsigned int> &numbersNeighbourCell1Ds) = 0;
    /// \brief Initialize the Cell0D Cell1D neighbours number
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param numberNeighbourCell1Ds the number of Cell1D neighbours of the Cell0D
    virtual void Cell0DInitializeNeighbourCell1Ds(const unsigned int &cell0DIndex, const unsigned int &numberNeighbourCell1Ds) = 0;
    virtual void Cell0DInitializeNeighbourCell1Ds(const unsigned int &cell0DIndex,
                                                  const std::vector<unsigned int> &neighbourCell1Ds) = 0;
    /// \brief Insert the Cell0D Cell1D neighbour
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of Cell1D neighbour of the Cell0D from 0 to
    /// Cell0DNumberNeighbourCell1D(cell0DIndex) \param neigbourCell1DIndex the Cell1D neighbour index from 0 to
    /// Cell1DTotalNumber() \note Cell0DInitializeNeighbourCell1Ds() shall be called before
    virtual void Cell0DInsertNeighbourCell1D(const unsigned int &cell0DIndex,
                                             const unsigned int &neighbourIndex,
                                             const unsigned int &neigbourCell1DIndex) = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the number of Neighbour Cell1Ds of Cell0D
    virtual unsigned int Cell0DNumberNeighbourCell1D(const unsigned int &cell0DIndex) const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell1D from 0 to Cell0DNumberNeighbourCell1D(cell0DIndex)
    /// \return the Cell1D index of Neighbour Cell1Ds of Cell0D from 0 to Cell1DTotalNumber()
    virtual unsigned int Cell0DNeighbourCell1D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const = 0;

    virtual std::vector<unsigned int> Cell0DNeighbourCell1Ds(const unsigned int &cell0DIndex) const = 0;

    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell1D from 0 to Cell0DNumberNeighbourCell1D(cell0DIndex)
    /// \return true if Neighbour Cell1Ds of Cell0D at position neighbourIndex exists
    virtual bool Cell0DHasNeighbourCell1D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const = 0;

    /// \brief Reset the Cell0D Cell1D neighbour to empty value (Cell0DHasNeighbourCell1D is false)
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of Cell1D neighbour of the Cell0D from 0 to
    /// Cell0DNumberNeighbourCell1D(cell0DIndex)
    virtual void Cell0DResetNeighbourCell1D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) = 0;

    virtual std::vector<std::vector<unsigned int>> Cell0DsNeighbourCell2Ds() const = 0;
    /// \brief Initialize the Cell0Ds Cell2D neighbours number
    /// \param numbersNeighbourCell2Ds the number of Cell2D neighbours of each Cell0D, size 1 x Cell0DTotalNumber()
    virtual void Cell0DsInitializeNeighbourCell2Ds(const std::vector<unsigned int> &numbersNeighbourCell2Ds) = 0;
    /// \brief Initialize the Cell0D Cell2D neighbours number
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param numberNeighbourCell2Ds the number of Cell2D neighbours of the Cell0D
    virtual void Cell0DInitializeNeighbourCell2Ds(const unsigned int &cell0DIndex, const unsigned int &numberNeighbourCell2Ds) = 0;
    virtual void Cell0DInitializeNeighbourCell2Ds(const unsigned int &cell0DIndex,
                                                  const std::vector<unsigned int> &neighbourCell2Ds) = 0;
    /// \brief Insert the Cell0D Cell2D neighbour
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of Cell2D neighbour of the Cell0D from 0 to
    /// Cell0DNumberNeighbourCell2D(cell0DIndex) \param neigbourCell2DIndex the Cell2D neighbour index from 0 to
    /// Cell2DTotalNumber() \note Cell0DInitializeNeighbourCell2Ds() shall be called before
    virtual void Cell0DInsertNeighbourCell2D(const unsigned int &cell0DIndex,
                                             const unsigned int &neighbourIndex,
                                             const unsigned int &neigbourCell2DIndex) = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the number of Neighbour Cell2Ds of Cell0D
    virtual unsigned int Cell0DNumberNeighbourCell2D(const unsigned int &cell0DIndex) const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell2D from 0 to Cell0DNumberNeighbourCell2D(cell0DIndex)
    /// \return the Cell2D index of Neighbour Cell2Ds of Cell0D from 0 to Cell2DTotalNumber()
    virtual unsigned int Cell0DNeighbourCell2D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const = 0;

    virtual std::vector<unsigned int> Cell0DNeighbourCell2Ds(const unsigned int &cell0DIndex) const = 0;

    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell2D from 0 to Cell0DNumberNeighbourCell2D(cell0DIndex)
    /// \return true if Neighbour Cell2Ds of Cell0D at position neighbourIndex exists
    virtual bool Cell0DHasNeighbourCell2D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const = 0;
    /// \brief Reset the Cell0D Cell2D neighbour to empty value (Cell0DHasNeighbourCell2D is false)
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of Cell2D neighbour of the Cell0D from 0 to
    /// Cell0DNumberNeighbourCell2D(cell0DIndex)
    virtual void Cell0DResetNeighbourCell2D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) = 0;

    virtual std::vector<std::vector<unsigned int>> Cell0DsNeighbourCell3Ds() const = 0;
    /// \brief Initialize the Cell0Ds Cell3D neighbours number
    /// \param numbersNeighbourCell3Ds the number of Cell3D neighbours of each Cell0D, size 1 x Cell0DTotalNumber()
    virtual void Cell0DsInitializeNeighbourCell3Ds(const std::vector<unsigned int> &numbersNeighbourCell3Ds) = 0;
    /// \brief Initialize the Cell0D Cell3D neighbours number
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param numberNeighbourCell3Ds the number of Cell3D neighbours of the Cell0D
    virtual void Cell0DInitializeNeighbourCell3Ds(const unsigned int &cell0DIndex, const unsigned int &numberNeighbourCell3Ds) = 0;
    virtual void Cell0DInitializeNeighbourCell3Ds(const unsigned int &cell0DIndex,
                                                  const std::vector<unsigned int> &neighbourCell3Ds) = 0;
    /// \brief Insert the Cell0D Cell3D neighbour
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of Cell3D neighbour of the Cell0D from 0 to
    /// Cell0DNumberNeighbourCell3D(cell0DIndex) \param neigbourCell3DIndex the Cell3D neighbour index from 0 to
    /// Cell3DTotalNumber() \note Cell0DInitializeNeighbourCell3Ds() shall be called before
    virtual void Cell0DInsertNeighbourCell3D(const unsigned int &cell0DIndex,
                                             const unsigned int &neighbourIndex,
                                             const unsigned int &neigbourCell3DIndex) = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the number of Neighbour Cell3Ds of Cell0D
    virtual unsigned int Cell0DNumberNeighbourCell3D(const unsigned int &cell0DIndex) const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell3D from 0 to Cell0DNumberNeighbourCell3D(cell0DIndex)
    /// \return the Cell3D index of Neighbour Cell3Ds of Cell0D from 0 to Cell3DTotalNumber()
    virtual unsigned int Cell0DNeighbourCell3D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const = 0;
    virtual std::vector<unsigned int> Cell0DNeighbourCell3Ds(const unsigned int &cell0DIndex) const = 0;
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell3D from 0 to Cell0DNumberNeighbourCell3D(cell0DIndex)
    /// \return true if Neighbour Cell3Ds of Cell0D at position neighbourIndex exists
    virtual bool Cell0DHasNeighbourCell3D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const = 0;
    /// \brief Reset the Cell0D Cell3D neighbour to empty value (Cell0DHasNeighbourCell3D is false)
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param neighbourIndex the number of Cell3D neighbour of the Cell0D from 0 to
    /// Cell0DNumberNeighbourCell3D(cell0DIndex)
    virtual void Cell0DResetNeighbourCell3D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) = 0;

    /// \brief Initialize the Cell0Ds double properties
    /// \param numberDoubleProperties the total number of Cell0Ds properties
    /// \note No reset of Cell0Ds is performed
    virtual void Cell0DInitializeDoubleProperties(const unsigned int &numberDoubleProperties) = 0;
    /// \brief Add the Cell0Ds double property identified by id
    /// \param propertyId the id of Cell0Ds property
    /// \return the double property position
    virtual unsigned int Cell0DAddDoubleProperty(const std::string &propertyId) = 0;
    /// \brief Initialize the Cell0Ds double property sizes
    /// \param propertyIndex the index of Cell0D double property from 0 to Cell0DNumberProperties()
    /// \param porpertySize the double property size of each Cell0D, size 1 x Cell0DTotalNumber()
    virtual void Cell0DsInitializeDoublePropertyValues(const unsigned int &propertyIndex,
                                                       const std::vector<unsigned int> &porpertySizes) = 0;
    /// \brief Initialize the Cell0Ds double property size
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param propertyIndex the index of Cell0D double property from 0 to Cell0DNumberProperties()
    /// \param porpertySize the double property size of Cell0D
    virtual void Cell0DInitializeDoublePropertyValues(const unsigned int &cell0DIndex,
                                                      const unsigned int &propertyIndex,
                                                      const unsigned int &porpertySize) = 0;
    /// \brief Insert the Cell0Ds double property value at position
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param propertyIndex the index of Cell0D double property from 0 to Cell0DNumberProperties()
    /// \param propertyIndex the index of Cell0D double property from 0 to Cell0DNumberProperties()
    virtual void Cell0DInsertDoublePropertyValue(const unsigned int &cell0DIndex,
                                                 const unsigned int &propertyIndex,
                                                 const unsigned int &propertyValueIndex,
                                                 const double &propertyValue) = 0;

    /// \return the total number of double properties of Cell0Ds
    virtual unsigned int Cell0DNumberDoubleProperties() const = 0;
    /// \return the id of the double property of Cell0Ds
    /// \param propertyIndex the index of Cell0D double property from 0 to Cell0DNumberProperties()
    virtual std::string Cell0DDoublePropertyId(const unsigned int &propertyIndex) const = 0;
    /// \return true if the double propertyId of Cell0Ds exists
    /// \param propertyId the id of Cell0D double property
    virtual bool Cell0DDoublePropertyExists(const std::string &propertyId) const = 0;
    /// \return the propertyIndex of the double property of Cell0Ds from 0 to Cell0DNumberProperties()
    /// \param propertyId the id of Cell0D double property
    virtual unsigned int Cell0DDoublePropertyIndex(const std::string &propertyId) const = 0;
    /// \return the size of the double property of Cell0D
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param propertyIndex the index of Cell0D double property from 0 to Cell0DNumberProperties()
    virtual unsigned int Cell0DDoublePropertySize(const unsigned int &cell0DIndex, const unsigned int &propertyIndex) const = 0;
    /// \return the value of the double property at valueIndex of Cell0D
    /// \param cell0DIndex the index of Cell0D from 0 to Cell0DTotalNumber()
    /// \param propertyIndex the index of Cell0D double property from 0 to Cell0DNumberProperties()
    /// \param propertyValueIndex the index of Cell0D double property value from 0 to Cell0DDoublePropertySize()
    virtual double Cell0DDoublePropertyValue(const unsigned int &cell0DIndex,
                                             const unsigned int &propertyIndex,
                                             const unsigned int &propertyValueIndex) const = 0;

    /// \brief Initialize the Cell1Ds container
    /// \param numberCell1Ds the total number of Cell1Ds
    /// \note No reset of Cell1Ds is performed
    virtual void Cell1DsInitialize(const unsigned int &numberCell1Ds) = 0;
    /// \brief Append Cell1Ds to the Cell1Ds container
    /// \param numberCell1Ds the number of Cell1Ds to append
    /// \return the previous number of Cell1Ds before the append operation
    virtual unsigned int Cell1DAppend(const unsigned int &numberCell1Ds) = 0;
    /// \brief Remove the Cell1D from the mesh
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \note the cell1D is removed and no integrity check in the mesh are performed
    virtual void Cell1DRemove(const unsigned int &cell1DIndex) = 0;
    /// \brief Set the Cell1D Origin and End
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param originCell0DIndex the Cell0D index of Cell1D origin from 0 to Cell0DTotalNumber()
    /// \param endCell0DIndex the Cell0D index of Cell1D end from 0 to Cell0DTotalNumber()
    virtual void Cell1DInsertExtremes(const unsigned int &cell1DIndex,
                                      const unsigned int &originCell0DIndex,
                                      const unsigned int &endCell0DIndex) = 0;

    /// \brief Set the Cell1D Extremes for the whole mesh edges
    /// \param cell1DExtremes the origin and end indices of all the edges, size 2 x Cell1DTotalNumber()
    virtual void Cell1DsInsertExtremes(const Eigen::MatrixXi &cell1DExtremes) = 0;

    /// \return the extrems as Eigen MatrixXi of cell1Ds, size 2xCell1DTotalNumber()
    virtual Eigen::MatrixXi Cell1DsExtremes() const = 0;
    /// \return the extrems as Eigen MatrixXi of cell1D, size 2
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    virtual Eigen::VectorXi Cell1DExtremes(const unsigned int &cell1DIndex) const = 0;
    virtual Eigen::MatrixXi Cell1DsExtremes(const std::vector<unsigned int> &cell1Ds) const = 0;
    /// \return the Cell1D Index if Cell1D (origin->end) exists, Cell1DTotalNumber() otherwise
    /// \param originCell0DIndex the Cell0D Id of origin from 0 to Cell0DTotalNumber()
    /// \param endCell0DIndex the Cell0D Id of origin from 0 to Cell0DTotalNumber()
    virtual unsigned int Cell1DByExtremes(const unsigned int &originCell0DIndex, const unsigned int &endCell0DIndex) const = 0;
    /// \brief Set the Cell1D Marker
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param marker the marker of the Cell1D
    virtual void Cell1DSetMarker(const unsigned int &cell1DIndex, const unsigned int &marker) = 0;
    /// \brief Set the Cell1D state
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param state true if Cell1D is active, false otherwise
    virtual void Cell1DSetState(const unsigned int &cell1DIndex, const bool &state) = 0;
    /// \return the total number of Cell1Ds
    virtual unsigned int Cell1DTotalNumber() const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \param vertexIndex the index of the vertex from 0 to 2
    /// \return Cell0D index of the vertex of Cell1D
    virtual unsigned int Cell1DVertex(const unsigned int &cell1DIndex, const unsigned int &vertexIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return the coordinates of Cell1D, size 3x2
    virtual Eigen::MatrixXd Cell1DCoordinates(const unsigned int &cell1DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return the origin coordinates of Cell1D
    virtual Eigen::Vector3d Cell1DOriginCoordinates(const unsigned int &cell1DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return the end coordinates of Cell1D
    virtual Eigen::Vector3d Cell1DEndCoordinates(const unsigned int &cell1DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return the origin Cell0D index of Cell1D from 0 to Cell0DTotalNumber()
    virtual unsigned int Cell1DOrigin(const unsigned int &cell1DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return the end Cell0D index of Cell1D from 0 to Cell0DTotalNumber()
    virtual unsigned int Cell1DEnd(const unsigned int &cell1DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the index of the cell0DIndex on the cell1D from 0 to 1, 2 if not found
    virtual unsigned int Cell1DFindExtreme(const unsigned int &cell1DIndex, const unsigned int &cell0DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return the cell1D marker
    virtual unsigned int Cell1DMarker(const unsigned int &cell1DIndex) const = 0;
    virtual std::vector<unsigned int> Cell1DsMarker() const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return if the cell1D is active
    virtual bool Cell1DIsActive(const unsigned int &cell1DIndex) const = 0;
    /// \return the activation state of all cell1Ds
    virtual std::vector<bool> Cell1DsState() const = 0;

    /// \param updatedCell1DIndex the updated cell1D index, from 0 to Cell1DTotalNumber()
    /// \return true if has an original cell, false otherwise (the original cell is itself)
    virtual bool Cell1DHasOriginalCell1D(const unsigned int &updatedCell1DIndex) const = 0;
    /// \param updatedCell1DIndex the updated cell1D index, from 0 to Cell1DTotalNumber()
    /// \return the original cell1D index, from 0 to Cell1DTotalNumber()
    virtual unsigned int Cell1DOriginalCell1D(const unsigned int &updatedCell1DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return if the cell1D has new cell1Ds associated
    virtual bool Cell1DHasUpdatedCell1Ds(const unsigned int &cell1DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return the number of new cell1Ds associated to cell1DIndex
    virtual unsigned int Cell1DNumberUpdatedCell1Ds(const unsigned int &cell1DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \param updatedCell1DIdex the index of the new Cell1D from 0 to Cell1DTotalNumber()
    /// \return if the Cell1D has the updatedCell1DIdex associated
    virtual bool Cell1DHasUpdatedCell1D(const unsigned int &cell1DIndex, const unsigned int &updatedCell1DIdex) const = 0;
    /// \brief Add the new Cell1D to an existing Cell1D
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param updatedCell1DIdex the index of the new Cell1D from 0 to Cell1DTotalNumber()
    virtual void Cell1DInsertUpdatedCell1D(const unsigned int &cell1DIndex, const unsigned int &updatedCell1DIdex) = 0;
    /// \brief return the updated Cell1D Ids for cell1DIndex
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param updatedCell1DIds the list of the new Cell1D Ids associated to cell1DIndex
    /// \return true if the cell1DIndex is contained in the updatedCell1DIds list, false otherwise
    virtual bool Cell1DUpdatedCell1Ds(const unsigned int &cell1DIndex, std::list<unsigned int> &updatedCell1DIds) const = 0;

    virtual std::vector<std::vector<unsigned int>> Cell1DsNeighbourCell2Ds() const = 0;
    /// \brief Initialize the Cell1Ds Cell2D neighbours number
    /// \param numbersNeighbourCell2Ds the number of Cell2D neighbours of each Cell1D, size 1 x Cell1DTotalNumber()
    virtual void Cell1DsInitializeNeighbourCell2Ds(const std::vector<unsigned int> &numbersNeighbourCell2Ds) = 0;
    /// \brief Initialize the Cell1Ds Cell2D neighbours number
    /// \param numberNeighbourCell2Ds the number of Cell2D neighbours of the Cell1D
    virtual void Cell1DsInitializeNeighbourCell2Ds(const unsigned int &numberNeighbourCell2Ds) = 0;

    /// \brief Initialize the Cell1D Cell2D neighbours number
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param numberNeighbourCell2Ds the number of Cell2D neighbours of the Cell1D
    virtual void Cell1DInitializeNeighbourCell2Ds(const unsigned int &cell1DIndex, const unsigned int &numberNeighbourCell2Ds) = 0;
    virtual void Cell1DInitializeNeighbourCell2Ds(const unsigned int &cell1DIndex,
                                                  const std::vector<unsigned int> &neighbourCell2Ds) = 0;
    /// \brief Insert the Cell1D Cell2D neighbour
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param neighbourIndex the number of Cell2D neighbour of the Cell1D from 0 to
    /// Cell1DNumberNeighbourCell2D(cell1DIndex) \param neigbourCell2DIndex the Cell2D neighbour index from 0 to
    /// Cell2DTotalNumber() \note Cell1DInitializeNeighbourCell2Ds() shall be called before
    virtual void Cell1DInsertNeighbourCell2D(const unsigned int &cell1DIndex,
                                             const unsigned int &neighbourIndex,
                                             const unsigned int &neigbourCell2DIndex) = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return the number of Neighbour Cell2Ds of Cell1D
    virtual unsigned int Cell1DNumberNeighbourCell2D(const unsigned int &cell1DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell2D from 0 to Cell1DNumberNeighbourCell2D(cell1DIndex)
    /// \return the Cell2D index of Neighbour Cell2Ds of Cell1D from 0 to Cell2DTotalNumber()
    virtual unsigned int Cell1DNeighbourCell2D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex) const = 0;

    virtual inline std::vector<unsigned int> Cell1DNeighbourCell2Ds(const unsigned int &cell1DIndex) const = 0;

    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell2D from 0 to Cell1DNumberNeighbourCell2D(cell1DIndex)
    /// \return true if Neighbour Cell2Ds of Cell1D at position neighbourIndex exists
    virtual bool Cell1DHasNeighbourCell2D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex) const = 0;
    /// \brief Reset the Cell1D Cell2D neighbour to empty value (Cell1DHasNeighbourCell2D is false)
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param neighbourIndex the number of Cell2D neighbour of the Cell1D from 0 to
    /// Cell1DNumberNeighbourCell2D(cell1DIndex)
    virtual void Cell1DResetNeighbourCell2D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex) = 0;

    virtual std::vector<std::vector<unsigned int>> Cell1DsNeighbourCell3Ds() const = 0;
    /// \brief Initialize the Cell1Ds Cell3D neighbours number
    /// \param numbersNeighbourCell3Ds the number of Cell2D neighbours of each Cell1D, size 1 x Cell1DTotalNumber()
    virtual void Cell1DsInitializeNeighbourCell3Ds(const std::vector<unsigned int> &numbersNeighbourCell3Ds) = 0;
    /// \brief Initialize the Cell1D Cell3D neighbours number
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param numberNeighbourCell3Ds the number of Cell3D neighbours of the Cell1D
    virtual void Cell1DInitializeNeighbourCell3Ds(const unsigned int &cell1DIndex, const unsigned int &numberNeighbourCell3Ds) = 0;
    virtual void Cell1DInitializeNeighbourCell3Ds(const unsigned int &cell1DIndex,
                                                  const std::vector<unsigned int> &neighbourCell3Ds) = 0;
    /// \brief Insert the Cell1D Cell3D neighbour
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param neighbourIndex the number of Cell3D neighbour of the Cell1D from 0 to
    /// Cell1DNumberNeighbourCell3D(cell1DIndex) \param neigbourCell3DIndex the Cell3D neighbour index from 0 to
    /// Cell3DTotalNumber() \note Cell1DInitializeNeighbourCell3Ds() shall be called before
    virtual void Cell1DInsertNeighbourCell3D(const unsigned int &cell1DIndex,
                                             const unsigned int &neighbourIndex,
                                             const unsigned int &neigbourCell3DIndex) = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return the number of Neighbour Cell3Ds of Cell1D
    virtual unsigned int Cell1DNumberNeighbourCell3D(const unsigned int &cell1DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell3D from 0 to Cell1DNumberNeighbourCell3D(cell1DIndex)
    /// \return the Cell3D index of Neighbour Cell3Ds of Cell1D from 0 to Cell3DTotalNumber()
    virtual unsigned int Cell1DNeighbourCell3D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex) const = 0;
    virtual std::vector<unsigned int> Cell1DNeighbourCell3Ds(const unsigned int &cell1DIndex) const = 0;
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell3D from 0 to Cell1DNumberNeighbourCell3D(cell1DIndex)
    /// \return true if Neighbour Cell3Ds of Cell1D at position neighbourIndex exists
    virtual bool Cell1DHasNeighbourCell3D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex) const = 0;
    /// \brief Reset the Cell1D Cell3D neighbour to empty value (Cell1DHasNeighbourCell3D is false)
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param neighbourIndex the number of Cell3D neighbour of the Cell1D from 0 to
    /// Cell1DNumberNeighbourCell3D(cell1DIndex)
    virtual void Cell1DResetNeighbourCell3D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex) = 0;

    /// \brief Initialize the Cell1Ds double properties
    /// \param numberDoubleProperties the total number of Cell1Ds properties
    /// \note No reset of Cell1Ds is performed
    virtual void Cell1DInitializeDoubleProperties(const unsigned int &numberDoubleProperties) = 0;
    /// \brief Add the Cell1Ds double property identified by id
    /// \param propertyId the id of Cell1Ds property
    /// \return the double property position
    virtual unsigned int Cell1DAddDoubleProperty(const std::string &propertyId) = 0;
    /// \brief Initialize the Cell1Ds double property sizes
    /// \param propertyIndex the index of Cell1D double property from 0 to Cell1DNumberProperties()
    /// \param porpertySize the double property size of each Cell1D, size 1 x Cell1DTotalNumber()
    virtual void Cell1DsInitializeDoublePropertyValues(const unsigned int &propertyIndex,
                                                       const std::vector<unsigned int> &porpertySizes) = 0;
    /// \brief Initialize the Cell1Ds double property size
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param propertyIndex the index of Cell1D double property from 0 to Cell1DNumberProperties()
    /// \param porpertySize the double property size of Cell1D
    virtual void Cell1DInitializeDoublePropertyValues(const unsigned int &cell1DIndex,
                                                      const unsigned int &propertyIndex,
                                                      const unsigned int &porpertySize) = 0;
    /// \brief Insert the Cell1Ds double property value at position
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param propertyIndex the index of Cell1D double property from 0 to Cell1DNumberProperties()
    /// \param propertyIndex the index of Cell1D double property from 0 to Cell1DNumberProperties()
    virtual void Cell1DInsertDoublePropertyValue(const unsigned int &cell1DIndex,
                                                 const unsigned int &propertyIndex,
                                                 const unsigned int &propertyValueIndex,
                                                 const double &propertyValue) = 0;

    /// \return the total number of double properties of Cell1Ds
    virtual unsigned int Cell1DNumberDoubleProperties() const = 0;
    /// \return the id of the double property of Cell1Ds
    /// \param propertyIndex the index of Cell1D double property from 0 to Cell1DNumberProperties()
    virtual std::string Cell1DDoublePropertyId(const unsigned int &propertyIndex) const = 0;
    /// \return true if the double propertyId of Cell1Ds exists
    /// \param propertyId the id of Cell1D double property
    virtual bool Cell1DDoublePropertyExists(const std::string &propertyId) const = 0;
    /// \return the propertyIndex of the double property of Cell1Ds from 0 to Cell1DNumberProperties()
    /// \param propertyId the id of Cell1D double property
    virtual unsigned int Cell1DDoublePropertyIndex(const std::string &propertyId) const = 0;
    /// \return the size of the double property of Cell1D
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param propertyIndex the index of Cell1D double property from 0 to Cell1DNumberProperties()
    virtual unsigned int Cell1DDoublePropertySize(const unsigned int &cell1DIndex, const unsigned int &propertyIndex) const = 0;
    /// \return the value of the double property at valueIndex of Cell1D
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param propertyIndex the index of Cell1D double property from 0 to Cell1DNumberProperties()
    /// \param propertyValueIndex the index of Cell1D double property value from 0 to Cell1DDoublePropertySize()
    virtual double Cell1DDoublePropertyValue(const unsigned int &cell1DIndex,
                                             const unsigned int &propertyIndex,
                                             const unsigned int &propertyValueIndex) const = 0;

    /// \brief Initialize the Cell2Ds container
    /// \param numberCell2Ds the total number of Cell2Ds
    /// \note No reset of Cell2Ds is performed
    virtual void Cell2DsInitialize(const unsigned int &numberCell2Ds) = 0;
    /// \brief Append Cell2Ds to the Cell2Ds container
    /// \param numberCell2Ds the number of Cell2Ds to append
    /// \return the previous number of Cell2Ds before the append operation
    virtual unsigned int Cell2DAppend(const unsigned int &numberCell2Ds) = 0;
    /// \brief Remove the Cell2D from the mesh
    /// \param cell2DIndex the index of Cell0D from 0 to Cell2DTotalNumber()
    /// \note the cell2D is removed and no integrity check in the mesh are performed
    virtual void Cell2DRemove(const unsigned int &cell2DIndex) = 0;
    /// \brief Initialize the Cell2Ds vertices number
    /// \param numberCell2DVertices the number of vertices of all Cell2Ds
    virtual void Cell2DsInitializeVertices(const unsigned int &numberCell2DVertices) = 0;
    /// \brief Initialize the Cell2Ds vertices number
    /// \param numberCell2DsVertices the number of vertices of each Cell2D
    virtual void Cell2DsInitializeVertices(const std::vector<unsigned int> &numberCell2DsVertices) = 0;
    /// \brief Initialize the Cell2D vertices  number
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param numberCell2DVertices the number of vertices of Cell2D
    virtual void Cell2DInitializeVertices(const unsigned int &cell2DIndex, const unsigned int &numberCell2DVertices) = 0;
    /// \brief Initialize the Cell2Ds edges number
    /// \param numberCell2DEdges the number of edges of all Cell2Ds
    virtual void Cell2DsInitializeEdges(const unsigned int &numberCell2DEdges) = 0;
    /// \brief Initialize the Cell2Ds edges number
    /// \param numberCell2DsEdges the number of edges of each Cell2D
    virtual void Cell2DsInitializeEdges(const std::vector<unsigned int> &numberCell2DsEdges) = 0;
    /// \brief Initialize the Cell2D edges number
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param numberCell2DEdges the number of edges of Cell2D
    virtual void Cell2DInitializeEdges(const unsigned int &cell2DIndex, const unsigned int &numberCell2DEdges) = 0;
    /// \brief Insert the Cell2D vertex
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param verticesCell0DIndices the Cell0D vertices index from 0 to Cell0DTotalNumber()
    /// \note Cell2DInitializeVertices() should be called before using this method
    virtual void Cell2DInsertVertices(const unsigned int &cell2DIndex, const std::vector<unsigned int> &verticesCell0DIndices) = 0;
    /// \brief Insert the Cell2D vertex
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param vertexIndex the number of vertex of the Cell2D from 0 to Cell2DNumberVertices(cell2DIndex)
    /// \param vertexCell0DIndex the Cell0D vertex index from 0 to Cell0DTotalNumber()
    /// \note Cell2DInitializeVertices() should be called before using this method
    virtual void Cell2DInsertVertex(const unsigned int &cell2DIndex,
                                    const unsigned int &vertexIndex,
                                    const unsigned int &vertexCell0DIndex) = 0;
    /// \brief Add the Cell2D vertices
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param vertexCell0DIndices the Cell0D vertices indices from 0 to Cell0DTotalNumber()
    /// \note No itialization is necessary
    virtual void Cell2DAddVertices(const unsigned int &cell2DIndex, const std::vector<unsigned int> &verticesCell0DIndices) = 0;
    /// \brief Insert the Cell2D edge
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param edgesCell1DIndices the Cell1D edges indices from 0 to Cell1DTotalNumber()
    /// \note Cell2DInitializeEdges() should be called before using this method
    virtual void Cell2DInsertEdges(const unsigned int &cell2DIndex, const std::vector<unsigned int> &edgesCell1DIndices) = 0;
    /// \brief Insert the Cell2D edge
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param edgeIndex the number of edge of the Cell2D from 0 to Cell2DNumberEdges(cell2DIndex)
    /// \param edgeCell0DIndex the Cell1D edge index from 0 to Cell1DTotalNumber()
    /// \note Cell2DInitializeEdges() should be called before using this method
    virtual void Cell2DInsertEdge(const unsigned int &cell2DIndex, const unsigned int &edgeIndex, const unsigned int &edgeCell1DIndex) = 0;
    /// \brief Add the Cell2D edges
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param edgesCell1DIndices the Cell1D edges indices from 0 to Cell1DTotalNumber()
    /// \note No itialization is necessary
    virtual void Cell2DAddEdges(const unsigned int &cell2DIndex, const std::vector<unsigned int> &edgesCell1DIndices) = 0;

    /// \brief Cell2D Add Vertices And Edges
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param verticesAndEdgesIndices the matrix of Cell0Ds and Cell1Ds indices
    /// \note No itialization is necessary
    virtual void Cell2DAddVerticesAndEdges(const unsigned int &cell2DIndex, const Eigen::MatrixXi &verticesAndEdgesIndices) = 0;

    /// \brief Set the Cell2D Marker
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param marker the marker of the Cell2D
    virtual void Cell2DSetMarker(const unsigned int &cell2DIndex, const unsigned int &marker) = 0;
    /// \brief Set the Cell1D state
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param state true if Cell1D is active, false otherwise
    virtual void Cell2DSetState(const unsigned int &cell2DIndex, const bool &state) = 0;
    /// \return the total number of Cell2Ds
    virtual unsigned int Cell2DTotalNumber() const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \return the number of vertices of Cell2D
    virtual unsigned int Cell2DNumberVertices(const unsigned int &cell2DIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \return the number of edges of Cell2D
    virtual unsigned int Cell2DNumberEdges(const unsigned int &cell2DIndex) const = 0;
    /// \return the Cell0D index collections of all Cell2Ds, size Cell2DTotalNumber() x
    /// Cell2DNumberVertices(cell2DIndex)
    virtual std::vector<std::vector<unsigned int>> Cell2DsVertices() const = 0;
    /// \return the Cell0Ds and Cell1Ds index collections of all Cell2Ds, size Cell2DTotalNumber() x (2 x
    /// Cell2DNumberVertices(cell2DIndex))
    virtual std::vector<Eigen::MatrixXi> Cell2DsExtremes() const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \return the Cell0D index collections of Cell2D from 0 to Cell0DTotalNumber(), size
    /// Cell2DNumberVertices(cell2DIndex)
    virtual std::vector<unsigned int> Cell2DVertices(const unsigned int &cell2DIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \param vertexIndex the index of cell0D vertex from 0 to NumberCell2DVertices(cell2DIndex)
    /// \return the Cell0D index of vertex of Cell2D from 0 to Cell0DTotalNumber()
    virtual unsigned int Cell2DVertex(const unsigned int &cell2DIndex, const unsigned int &vertexIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \param vertexIndex the index of cell0D vertex from 0 to NumberCell2DVertices(cell2DIndex)
    /// \return the Cell0D coordinates of vertex of Cell2D, size 3 x 1
    virtual Eigen::Vector3d Cell2DVertexCoordinates(const unsigned int &cell2DIndex, const unsigned int &vertexIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \return the Cell0D coordinates of all the vertices of Cell2D, size 3 x NumberCell2DVertices(cell2DIndex)
    virtual Eigen::MatrixXd Cell2DVerticesCoordinates(const unsigned int &cell2DIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the index of the cell0DIndex on the cell2D from 0 to NumberCell2DVertices(cell2DIndex),
    /// NumberCell2DVertices(cell2DIndex) if not found
    virtual unsigned int Cell2DFindVertex(const unsigned int &cell2DIndex, const unsigned int &cell0DIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \return the Cell1D index collections of Cell2D from 0 to Cell1DTotalNumber(), size
    /// Cell2DNumberEdges(cell2DIndex)
    virtual std::vector<unsigned int> Cell2DEdges(const unsigned int &cell2DIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \param edgeIndex the index of cell1D edge from 0 to NumberCell2DEdges(cell2DIndex)
    /// \return the Cell1D index of edge of Cell2D from 0 to Cell1DTotalNumber()
    virtual unsigned int Cell2DEdge(const unsigned int &cell2DIndex, const unsigned int &edgeIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return the index of the cell1DIndex on the cell2D from 0 to NumberCell2DEdges(cell2DIndex),
    /// NumberCell2DEdges(cell2DIndex) if not found
    virtual unsigned int Cell2DFindEdge(const unsigned int &cell2DIndex, const unsigned int &cell1DIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \param originCell0DIndex the Cell0D Id of origin from 0 to Cell0DTotalNumber()
    /// \param endCell0DIndex the Cell0D Id of origin from 0 to Cell0DTotalNumber()
    /// \return the index of the cell1DIndex on the cell2D from 0 to NumberCell2DEdges(cell2DIndex),
    /// NumberCell2DEdges(cell2DIndex) otherwise
    virtual unsigned int Cell2DFindEdgeByExtremes(const unsigned int &cell2DIndex,
                                                  const unsigned int &originCell0DIndex,
                                                  const unsigned int &endCell0DIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \return the cell2D marker
    virtual unsigned int Cell2DMarker(const unsigned int &cell2DIndex) const = 0;
    virtual std::vector<unsigned int> Cell2DsMarker() const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \return if the cell2D is active
    virtual bool Cell2DIsActive(const unsigned int &cell2DIndex) const = 0;
    /// \return the activation state of all cell2Ds
    virtual std::vector<bool> Cell2DsState() const = 0;

    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \return if the cell2D has new cell2Ds associated
    virtual bool Cell2DHasUpdatedCell2Ds(const unsigned int &cell2DIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \return the number of new cell2Ds associated to cell2DIndex
    virtual unsigned int Cell2DNumberUpdatedCell2Ds(const unsigned int &cell2DIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \param updatedCell2DIndex the index of the new Cell2D from 0 to Cell2DTotalNumber()
    /// \return if the Cell2D has the updatedCell2DIdex associated
    virtual bool Cell2DHasUpdatedCell2D(const unsigned int &cell2DIndex, const unsigned int &updatedCell2DIndex) const = 0;
    /// \brief Add the new Cell2D to an existing Cell2D
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param updatedCell2DIdex the index of the new Cell2D from 0 to Cell2DTotalNumber()
    virtual void Cell2DInsertUpdatedCell2D(const unsigned int &cell2DIndex, const unsigned int &updatedCell2DIdex) = 0;

    /// \param updatedCell2DIndex the updated cell2D index, from 0 to Cell2DTotalNumber()
    /// \return true if has an original cell, false otherwise (the original cell is itself)
    virtual bool Cell2DHasOriginalCell2D(const unsigned int &updatedCell2DIndex) const = 0;
    /// \param updatedCell2DIndex the updated cell2D index, from 0 to Cell2DTotalNumber()
    /// \return the original cell2D index, from 0 to Cell2DTotalNumber()
    virtual unsigned int Cell2DOriginalCell2D(const unsigned int &updatedCell2DIndex) const = 0;
    /// \brief return the updated Cell2D Ids for cell2DIndex
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param updatedCell2DIds the list of the new Cell2D Ids associated to cell2DIndex
    /// \return true if the cell2DIndex is contained in the updatedCell2DIds list, false otherwise
    virtual bool Cell2DUpdatedCell2Ds(const unsigned int &cell2DIndex, std::list<unsigned int> &updatedCell2DIds) const = 0;

    virtual std::vector<std::vector<unsigned int>> Cell2DsNeighbourCell3Ds() const = 0;
    /// \brief Initialize the Cell2Ds Cell3D neighbours number
    /// \param numbersNeighbourCell3Ds the number of Cell3D neighbours of each Cell2D, size 1 x Cell2DTotalNumber()
    virtual void Cell2DsInitializeNeighbourCell3Ds(const std::vector<unsigned int> &numbersNeighbourCell3Ds) = 0;
    /// \brief Initialize the Cell2D Cell3D neighbours number
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param numberNeighbourCell3Ds the number of Cell3D neighbours of the Cell2D
    virtual void Cell2DInitializeNeighbourCell3Ds(const unsigned int &cell2DIndex, const unsigned int &numberNeighbourCell3Ds) = 0;
    virtual void Cell2DInitializeNeighbourCell3Ds(const unsigned int &cell2DIndex,
                                                  const std::vector<unsigned int> &neighbourCell3Ds) = 0;
    /// \brief Insert the Cell2D Cell3D neighbour
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param neighbourIndex the number of Cell3D neighbour of the Cell2D from 0 to
    /// Cell2DNumberNeighbourCell3D(cell2DIndex) \param neigbourCell3DIndex the Cell3D neighbour index from 0 to
    /// Cell3DTotalNumber() \note Cell2DInitializeNeighbourCell3Ds() shall be called before
    virtual void Cell2DInsertNeighbourCell3D(const unsigned int &cell2DIndex,
                                             const unsigned int &neighbourIndex,
                                             const unsigned int &neigbourCell3DIndex) = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \return the number of Neighbour Cell3Ds of Cell2D
    virtual unsigned int Cell2DNumberNeighbourCell3D(const unsigned int &cell2DIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell3D from 0 to Cell2DNumberNeighbourCell3D(cell2DIndex)
    /// \return the Cell3D index of Neighbour Cell3Ds of Cell2D from 0 to Cell3DTotalNumber()
    virtual unsigned int Cell2DNeighbourCell3D(const unsigned int &cell2DIndex, const unsigned int &neighbourIndex) const = 0;
    virtual std::vector<unsigned int> Cell2DNeighbourCell3Ds(const unsigned int &cell2DIndex) const = 0;
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \param neighbourIndex the number of neigbourh Cell3D from 0 to Cell2DNumberNeighbourCell3D(cell2DIndex)
    /// \return true if Neighbour Cell3Ds of Cell2D at position neighbourIndex exists
    virtual bool Cell2DHasNeighbourCell3D(const unsigned int &cell2DIndex, const unsigned int &neighbourIndex) const = 0;
    /// \brief Reset the Cell2D Cell3D neighbour to empty value (Cell2DHasNeighbourCell3D is false)
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param neighbourIndex the number of Cell3D neighbour of the Cell2D from 0 to
    /// Cell2DNumberNeighbourCell3D(cell2DIndex)
    virtual void Cell2DResetNeighbourCell3D(const unsigned int &cell2DIndex, const unsigned int &neighbourIndex) = 0;

    /// \brief Initialize the Cell2Ds double properties
    /// \param numberDoubleProperties the total number of Cell2Ds properties
    /// \note No reset of Cell2Ds is performed
    virtual void Cell2DInitializeDoubleProperties(const unsigned int &numberDoubleProperties) = 0;
    /// \brief Add the Cell2Ds double property identified by id
    /// \param propertyId the id of Cell2Ds property
    /// \return the double property position
    virtual unsigned int Cell2DAddDoubleProperty(const std::string &propertyId) = 0;
    /// \brief Initialize the Cell2Ds double property sizes
    /// \param propertyIndex the index of Cell2D double property from 0 to Cell2DNumberProperties()
    /// \param porpertySize the double property size of each Cell2D, size 1 x Cell2DTotalNumber()
    virtual void Cell2DsInitializeDoublePropertyValues(const unsigned int &propertyIndex,
                                                       const std::vector<unsigned int> &porpertySizes) = 0;
    /// \brief Initialize the Cell2Ds double property size
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param propertyIndex the index of Cell2D double property from 0 to Cell2DNumberProperties()
    /// \param porpertySize the double property size of Cell2D
    virtual void Cell2DInitializeDoublePropertyValues(const unsigned int &cell2DIndex,
                                                      const unsigned int &propertyIndex,
                                                      const unsigned int &porpertySize) = 0;
    /// \brief Insert the Cell2Ds double property value at position
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param propertyIndex the index of Cell2D double property from 0 to Cell2DNumberProperties()
    /// \param propertyIndex the index of Cell2D double property from 0 to Cell2DNumberProperties()
    virtual void Cell2DInsertDoublePropertyValue(const unsigned int &cell2DIndex,
                                                 const unsigned int &propertyIndex,
                                                 const unsigned int &propertyValueIndex,
                                                 const double &propertyValue) = 0;

    /// \return the total number of double properties of Cell2Ds
    virtual unsigned int Cell2DNumberDoubleProperties() const = 0;
    /// \return the id of the double property of Cell2Ds
    /// \param propertyIndex the index of Cell2D double property from 0 to Cell2DNumberProperties()
    virtual std::string Cell2DDoublePropertyId(const unsigned int &propertyIndex) const = 0;
    /// \return true if the double propertyId of Cell2Ds exists
    /// \param propertyId the id of Cell2D double property
    virtual bool Cell2DDoublePropertyExists(const std::string &propertyId) const = 0;
    /// \return the propertyIndex of the double property of Cell2Ds from 0 to Cell2DNumberProperties()
    /// \param propertyId the id of Cell2D double property
    virtual unsigned int Cell2DDoublePropertyIndex(const std::string &propertyId) const = 0;
    /// \return the size of the double property of Cell2D
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param propertyIndex the index of Cell2D double property from 0 to Cell2DNumberProperties()
    virtual unsigned int Cell2DDoublePropertySize(const unsigned int &cell2DIndex, const unsigned int &propertyIndex) const = 0;
    /// \return the value of the double property at valueIndex of Cell2D
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param propertyIndex the index of Cell2D double property from 0 to Cell2DNumberProperties()
    /// \param propertyValueIndex the index of Cell2D double property value from 0 to Cell2DDoublePropertySize()
    virtual double Cell2DDoublePropertyValue(const unsigned int &cell2DIndex,
                                             const unsigned int &propertyIndex,
                                             const unsigned int &propertyValueIndex) const = 0;

    /// \brief Initialize the Cell2D subdivision number for each Cell2D
    /// \param numberSubDivisions the number of sub-polygons for each Cell2D, size 1 x Cell2DTotalNumber()
    /// \note each subdivision is a triangle, thus numberSubDivision shall be a multiple of 3
    virtual void Cell2DsInitializeSubDivision(const std::vector<unsigned int> &numberSubDivisions) = 0;
    /// \brief Initialize the Cell2D subdivision number
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param numberSubDivision the number of sub-polygons of the Cell2D
    /// \note each subdivision is a triangle, thus numberSubDivision shall be a multiple of 3
    virtual void Cell2DInitializeSubDivision(const unsigned int &cell2DIndex, const unsigned int &numberSubDivision) = 0;
    /// \brief Insert the subDivision vertex index
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param subDivisionIndex the subDivision index, from 0 to Cell2DNumberSubDivision(cell2DIndex)
    /// \param cell2DVertexIndex the Cell2D vertex index of the subDivision, from 0 to Cell0DTotalNumber()
    /// \note each subdivision is a triangle
    virtual void Cell2DInsertSubDivision(const unsigned int &cell2DIndex,
                                         const unsigned int &subDivisionIndex,
                                         const unsigned int &cell2DVertexIndex) = 0;

    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \return the total number of vertices of sub-polygons contained in the subdivision, a multiple of 3
    /// \note each subdivision is a triangle
    virtual unsigned int Cell2DNumberSubDivision(const unsigned int &cell2DIndex) const = 0;

    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param subDivisionIndex the subDivision index, from 0 to Cell2DNumberSubDivision(cell2DIndex)
    /// \return the Cell0D index of sub-polygons contained in the subdivision, from 0 to Cell0DTotalNumber()
    /// \note each sub-division shall be a triangle
    virtual unsigned int Cell2DSubDivisionCell0D(const unsigned int &cell2DIndex, const unsigned int &subDivisionIndex) const = 0;

    /// \brief Initialize the Cell3Ds container
    /// \param numberCell3Ds the total number of Cell3Ds
    /// \note No reset of Cell3Ds is performed
    virtual void Cell3DsInitialize(const unsigned int &numberCell3Ds) = 0;
    /// \brief Append Cell3Ds to the Cell3Ds container
    /// \param numberCell3Ds the number of Cell3Ds to append
    /// \return the previous number of Cell3Ds before the append operation
    virtual unsigned int Cell3DAppend(const unsigned int &numberCell3Ds) = 0;
    /// \brief Remove the Cell3D from the mesh
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \note the cell3D is removed and no integrity check in the mesh are performed
    virtual void Cell3DRemove(const unsigned int &cell3DIndex) = 0;
    /// \brief Initialize the Cell3Ds vertices number
    /// \param numberCell3DsVertices the number of vertices of each Cell3D
    virtual void Cell3DsInitializeVertices(const std::vector<unsigned int> &numberCell3DsVertices) = 0;
    /// \brief Initialize the Cell3D vertices  number
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param numberCell3DVertices the number of vertices of Cell3D
    virtual void Cell3DInitializeVertices(const unsigned int &cell3DIndex, const unsigned int &numberCell3DVertices) = 0;
    /// \brief Initialize the Cell3Ds edges number
    /// \param numberCell3DsEdges the number of edges of each Cell3D
    virtual void Cell3DsInitializeEdges(const std::vector<unsigned int> &numberCell3DsEdges) = 0;
    /// \brief Initialize the Cell3D edges number
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param numberCell3DEdges the number of edges of Cell3D
    virtual void Cell3DInitializeEdges(const unsigned int &cell3DIndex, const unsigned int &numberCell3DEdges) = 0;
    /// \brief Initialize the Cell3Ds faces number
    /// \param numberCell3DsFaces the number of faces of each Cell3D
    virtual void Cell3DsInitializeFaces(const std::vector<unsigned int> &numberCell3DsFaces) = 0;
    /// \brief Initialize the Cell3D faces number
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param numberCell3DFaces the number of faces of Cell3D
    virtual void Cell3DInitializeFaces(const unsigned int &cell3DIndex, const unsigned int &numberCell3DFaces) = 0;
    /// \brief Insert the Cell3D vertex
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param vertexIndex the number of vertex of the Cell3D from 0 to Cell3DNumberVertices(cell3DIndex)
    /// \param vertexCell0DIndex the Cell0D vertex index from 0 to Cell0DTotalNumber()
    /// \note Cell3DInitializeVertices() should be called before using this method
    virtual void Cell3DInsertVertex(const unsigned int &cell3DIndex,
                                    const unsigned int &vertexIndex,
                                    const unsigned int &vertexCell0DIndex) = 0;
    /// \brief Insert the Cell3D vertex
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param verticesCell0DIndices the Cell0D vertices index from 0 to Cell0DTotalNumber()
    /// \note Cell3DInitializeVertices() should be called before using this method
    virtual void Cell3DInsertVertices(const unsigned int &cell3DIndex, const std::vector<unsigned int> &verticesCell0DIndices) = 0;
    /// \brief Add the Cell3D vertices
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param vertexCell0DIndices the Cell0D vertices indices from 0 to Cell0DTotalNumber()
    /// \note No itialization is necessary
    virtual void Cell3DAddVertices(const unsigned int &cell3DIndex, const std::vector<unsigned int> &verticesCell0DIndices) = 0;
    /// \brief Insert the Cell3D edge
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param edgeIndex the number of edge of the Cell3D from 0 to Cell3DNumberEdges(cell3DIndex)
    /// \param edgeCell0DIndex the Cell1D edge index from 0 to Cell1DTotalNumber()
    /// \note Cell3DInitializeEdges() should be called before using this method
    virtual void Cell3DInsertEdge(const unsigned int &cell3DIndex, const unsigned int &edgeIndex, const unsigned int &edgeCell1DIndex) = 0;
    /// \brief Insert the Cell3D edges
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param edgesCell1DIndices the Cell1Ds index from 0 to Cell1DTotalNumber()
    /// \note Cell3DInitializeEdges() should be called before using this method
    virtual void Cell3DInsertEdges(const unsigned int &cell3DIndex, const std::vector<unsigned int> &edgesCell1DIndices) = 0;
    /// \brief Add the Cell3D edges
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param edgesCell0DIndices the Cell1D edges indices from 0 to Cell1DTotalNumber()
    /// \note No itialization is necessary
    virtual void Cell3DAddEdges(const unsigned int &cell3DIndex, const std::vector<unsigned int> &edgesCell0DIndices) = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \param cell0DIndex the index of cell0D from 0 to Cell0DTotalNumber()
    /// \return the index of the cell0DIndex on the cell3D from 0 to NumberCell3DVertices(cell3DIndex),
    /// NumberCell3DVertices(cell3DIndex) if not found
    virtual unsigned int Cell3DFindVertex(const unsigned int &cell3DIndex, const unsigned int &cell0DIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \param cell1DIndex the index of cell1D from 0 to Cell1DTotalNumber()
    /// \return the index of the cell1DIndex on the cell3D from 0 to NumberCell3DEdges(cell3DIndex),
    /// NumberCell3DEdges(cell3DIndex) if not found
    virtual unsigned int Cell3DFindEdge(const unsigned int &cell3DIndex, const unsigned int &cell1DIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \param cell2DIndex the index of cell2D from 0 to Cell2DTotalNumber()
    /// \return the index of the cell2DIndex on the cell3D from 0 to NumberCell3DFaces(cell3DIndex),
    /// NumberCell3DFaces(cell3DIndex) if not found
    virtual unsigned int Cell3DFindFace(const unsigned int &cell3DIndex, const unsigned int &cell2DIndex) const = 0;

    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \param originCell0DIndex the Cell0D Id of origin from 0 to Cell0DTotalNumber()
    /// \param endCell0DIndex the Cell0D Id of origin from 0 to Cell0DTotalNumber()
    /// \return the index of the cell1DIndex on the cell2D from 0 to NumberCell2DEdges(cell3DIndex),
    /// NumberCell2DEdges(cell2DIndex) otherwise
    virtual unsigned int Cell3DFindEdgeByExtremes(const unsigned int &cell3DIndex,
                                                  const unsigned int &originCell0DIndex,
                                                  const unsigned int &endCell0DIndex) const = 0;
    /// \brief Insert the Cell3D face
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param faceIndex the number of face of the Cell3D from 0 to Cell3DNumberFaces(cell3DIndex)
    /// \param faceCell0DIndex the Cell2D face index from 0 to Cell2DTotalNumber()
    /// \note Cell3DInitializeFaces() should be called before using this method
    virtual void Cell3DInsertFace(const unsigned int &cell3DIndex, const unsigned int &faceIndex, const unsigned int &faceCell2DIndex) = 0;
    /// \brief Insert the Cell3D faces
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param facesCell2DIndices the Cell2D index from 0 to Cell2DTotalNumber()
    /// \note Cell3DInitializeFaces() should be called before using this method
    virtual void Cell3DInsertFaces(const unsigned int &cell3DIndex, const std::vector<unsigned int> &facesCell2DIndices) = 0;
    /// \brief Add the Cell3D faces
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param facesCell0DIndices the Cell2D faces indices from 0 to Cell2DTotalNumber()
    /// \note No itialization is necessary
    virtual void Cell3DAddFaces(const unsigned int &cell3DIndex, const std::vector<unsigned int> &facesCell0DIndices) = 0;
    /// \brief Set the Cell1D Marker
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param marker the marker of the Cell3D
    virtual void Cell3DSetMarker(const unsigned int &cell3DIndex, const unsigned int &marker) = 0;
    /// \brief Set the Cell3D state
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param state true if Cell3D is active, false otherwise
    virtual void Cell3DSetState(const unsigned int &cell3DIndex, const bool &state) = 0;
    /// \return the total number of Cell3Ds
    virtual unsigned int Cell3DTotalNumber() const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \return the number of vertices of Cell3D
    virtual unsigned int Cell3DNumberVertices(const unsigned int &cell3DIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \return the number of edges of Cell3D
    virtual unsigned int Cell3DNumberEdges(const unsigned int &cell3DIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \return the number of faces of Cell3D
    virtual unsigned int Cell3DNumberFaces(const unsigned int &cell3DIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \return the Cell0D index collections of Cell3D from 0 to Cell0DTotalNumber(), size
    /// Cell3DNumberVertices(cell3DIndex)
    virtual std::vector<unsigned int> Cell3DVertices(const unsigned int &cell3DIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \param vertexIndex the index of cell0D vertex from 0 to NumberCell3DVertices(cell3DIndex)
    /// \return the Cell0D index of vertex of Cell3D from 0 to Cell0DTotalNumber()
    virtual unsigned int Cell3DVertex(const unsigned int &cell3DIndex, const unsigned int &vertexIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \param vertexIndex the index of cell0D vertex from 0 to NumberCell3DVertices(cell3DIndex)
    /// \return the Cell0D coordinates of vertex of Cell3D, size 3 x 1
    virtual Eigen::Vector3d Cell3DVertexCoordinates(const unsigned int &cell3DIndex, const unsigned int &vertexIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \return the Cell0D coordinates of all the vertices of Cell3D, size 3 x NumberCell3DVertices(cell3DIndex)
    virtual Eigen::MatrixXd Cell3DVerticesCoordinates(const unsigned int &cell3DIndex) const = 0;

    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \return the Cell1D index collections of Cell3D from 0 to Cell1DTotalNumber(), size
    /// Cell3DNumberEdges(cell3DIndex)
    virtual std::vector<unsigned int> Cell3DEdges(const unsigned int &cell3DIndex) const = 0;

    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \param edgeIndex the index of cell1D edge from 0 to NumberCell3DEdges(edgeIndex)
    /// \return the Cell1D index of edge of Cell3D from 0 to Cell1DTotalNumber()
    virtual unsigned int Cell3DEdge(const unsigned int &cell3DIndex, const unsigned int &edgeIndex) const = 0;

    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \return the Cell2D index collections of Cell3D from 0 to Cell2DTotalNumber(), size
    /// Cell3DNumberFaces(cell3DIndex)
    virtual std::vector<unsigned int> Cell3DFaces(const unsigned int &cell3DIndex) const = 0;

    /// \return the Cell0D index collections of all the faces of all Cell3Ds, size Cell3DTotalNumber() x
    /// Cell3DNumberFaces(cell3DIndex) x Cell2DNumberVertices(cell2DIndex)
    virtual std::vector<std::vector<std::vector<unsigned int>>> Cell3DsFacesVertices() const = 0;
    /// \return the Cell0D index collections of all Cell3Ds, size Cell3DTotalNumber() x
    /// Cell3DNumberVertices(cell3DIndex)
    virtual std::vector<std::vector<unsigned int>> Cell3DsVertices() const = 0;
    /// \return the Cell1D index collections of all Cell3Ds, size Cell3DTotalNumber() x Cell3DNumberEdges(cell3DIndex)
    virtual std::vector<std::vector<unsigned int>> Cell3DsEdges() const = 0;
    /// \return the Cell2D index collections of all Cell3Ds, size Cell3DTotalNumber() x Cell3DNumberFaces(cell3DIndex)
    virtual std::vector<std::vector<unsigned int>> Cell3DsFaces() const = 0;

    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \param faceIndex the index of cell2D face from 0 to NumberCell3DFaces(cell3DIndex)
    /// \return the Cell2D index of face of Cell3D from 0 to Cell2DTotalNumber()
    virtual unsigned int Cell3DFace(const unsigned int &cell3DIndex, const unsigned int &faceIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \return the cell3D marker
    virtual unsigned int Cell3DMarker(const unsigned int &cell3DIndex) const = 0;
    virtual std::vector<unsigned int> Cell3DsMarker() const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \return if the cell3D is active
    virtual bool Cell3DIsActive(const unsigned int &cell3DIndex) const = 0;
    /// \return the activation state of all cell3Ds
    virtual std::vector<bool> Cell3DsState() const = 0;

    /// \param updatedCell3DIndex the updated cell3D index, from 0 to Cell3DTotalNumber()
    /// \return true if has an original cell, false otherwise (the original cell is itself)
    virtual bool Cell3DHasOriginalCell3D(const unsigned int &updatedCell3DIndex) const = 0;
    /// \param updatedCell3DIndex the updated cell3D index, from 0 to Cell3DTotalNumber()
    /// \return the original cell3D index, from 0 to Cell3DTotalNumber()
    virtual unsigned int Cell3DOriginalCell3D(const unsigned int &updatedCell3DIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \return if the cell3D has new cell3Ds associated
    virtual bool Cell3DHasUpdatedCell3Ds(const unsigned int &cell3DIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \return the number of new cell3Ds associated to cell3DIndex
    virtual unsigned int Cell3DNumberUpdatedCell3Ds(const unsigned int &cell3DIndex) const = 0;
    /// \param cell3DIndex the index of cell3D from 0 to Cell3DTotalNumber()
    /// \param updatedCell3DIdex the index of the new Cell3D from 0 to Cell3DTotalNumber()
    /// \return if the Cell3D has the updatedCell3DIdex associated
    virtual bool Cell3DHasUpdatedCell3D(const unsigned int &cell3DIndex, const unsigned int &updatedCell3DIdex) const = 0;
    /// \brief Add the new Cell3D to an existing Cell3D
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param updatedCell3DIdex the index of the new Cell3D from 0 to Cell3DTotalNumber()
    virtual void Cell3DInsertUpdatedCell3D(const unsigned int &cell3DIndex, const unsigned int &updatedCell3DIdex) = 0;
    /// \brief return the updated Cell3D Ids for cell3DIndex
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param updatedCell3DIds the list of the new Cell3D Ids associated to cell3DIndex
    /// \return true if the cell3DIndex is contained in the updatedCell3DIds list, false otherwise
    virtual bool Cell3DUpdatedCell3Ds(const unsigned int &cell3DIndex, std::list<unsigned int> &updatedCell3DIds) const = 0;

    /// \brief Initialize the Cell3Ds double properties
    /// \param numberDoubleProperties the total number of Cell3Ds properties
    /// \note No reset of Cell3Ds is performed
    virtual void Cell3DInitializeDoubleProperties(const unsigned int &numberDoubleProperties) = 0;
    /// \brief Add the Cell3Ds double property identified by id
    /// \param propertyId the id of Cell3Ds property
    /// \return the double property position
    virtual unsigned int Cell3DAddDoubleProperty(const std::string &propertyId) = 0;
    /// \brief Initialize the Cell3Ds double property sizes
    /// \param propertyIndex the index of Cell3D double property from 0 to Cell3DNumberProperties()
    /// \param porpertySize the double property size of each Cell3D, size 1 x Cell3DTotalNumber()
    virtual void Cell3DsInitializeDoublePropertyValues(const unsigned int &propertyIndex,
                                                       const std::vector<unsigned int> &porpertySizes) = 0;
    /// \brief Initialize the Cell3Ds double property size
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param propertyIndex the index of Cell3D double property from 0 to Cell3DNumberProperties()
    /// \param porpertySize the double property size of Cell3D
    virtual void Cell3DInitializeDoublePropertyValues(const unsigned int &cell3DIndex,
                                                      const unsigned int &propertyIndex,
                                                      const unsigned int &porpertySize) = 0;
    /// \brief Insert the Cell3Ds double property value at position
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param propertyIndex the index of Cell3D double property from 0 to Cell3DNumberProperties()
    /// \param propertyIndex the index of Cell3D double property from 0 to Cell3DNumberProperties()
    virtual void Cell3DInsertDoublePropertyValue(const unsigned int &cell3DIndex,
                                                 const unsigned int &propertyIndex,
                                                 const unsigned int &propertyValueIndex,
                                                 const double &propertyValue) = 0;

    /// \return the total number of double properties of Cell3Ds
    virtual unsigned int Cell3DNumberDoubleProperties() const = 0;
    /// \return the id of the double property of Cell3Ds
    /// \param propertyIndex the index of Cell3D double property from 0 to Cell3DNumberProperties()
    virtual std::string Cell3DDoublePropertyId(const unsigned int &propertyIndex) const = 0;
    /// \return true if the double propertyId of Cell3Ds exists
    /// \param propertyId the id of Cell3D double property
    virtual bool Cell3DDoublePropertyExists(const std::string &propertyId) const = 0;
    /// \return the propertyIndex of the double property of Cell3Ds from 0 to Cell3DNumberProperties()
    /// \param propertyId the id of Cell3D double property
    virtual unsigned int Cell3DDoublePropertyIndex(const std::string &propertyId) const = 0;
    /// \return the size of the double property of Cell3D
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param propertyIndex the index of Cell3D double property from 0 to Cell3DNumberProperties()
    virtual unsigned int Cell3DDoublePropertySize(const unsigned int &cell3DIndex, const unsigned int &propertyIndex) const = 0;
    /// \return the value of the double property at valueIndex of Cell3D
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param propertyIndex the index of Cell3D double property from 0 to Cell3DNumberProperties()
    /// \param propertyValueIndex the index of Cell3D double property value from 0 to Cell3DDoublePropertySize()
    virtual double Cell3DDoublePropertyValue(const unsigned int &cell3DIndex,
                                             const unsigned int &propertyIndex,
                                             const unsigned int &propertyValueIndex) const = 0;

    /// \brief Compact the mesh to save memory
    virtual void Compress() = 0;

    /// \return The mesh converted to string
    virtual std::string ToString() = 0;
};
} // namespace Gedim

namespace Gedim
{
class MeshMatricesDAO : public IMeshDAO
{
  private:
    Gedim::MeshMatrices &_mesh;

    /// \brief for each i,j element of sparse matrix A if A[i,j] > minElement then A[i, j]--
    /// if A[i, j] == minElement the A[i, j] = 0
    /// \param matrix the sparse matrix A
    /// \param minElement the minElement
    /// \param newElementInitialization the new element initialization
    template <typename T> void AlignSparseMatrixHigherElements(Eigen::SparseMatrix<T> &matrix, const T &minElement);

    /// \brief for each i element of container on each map key v if v[i] > minElement then v[i]--
    /// if v[i] == minElement the v[i] = newElementInitialization
    /// \param elements the container map v
    /// \param minElement the minElement
    /// \param newElementInitialization the new element initialization
    template <class Container, class T>
    void AlignMapContainerHigherElements(std::unordered_map<unsigned int, Container> &elements,
                                         const T &minElement,
                                         const T &newElementInitialization);

    /// \brief for each i element of container v if v[i] > minElement then v[i]--
    /// if v[i] == minElement the v[i] = newElementInitialization
    /// \param elements the container v
    /// \param minElement the minElement
    /// \param newElementInitialization the new element initialization
    template <class Container, class T>
    void AlignContainerHigherElements(Container &elements, const T &minElement, const T &newElementInitialization);

    /// \brief for each i element of container v if v[i] == element then v[i] = newElementInitialization
    /// \param elements the container v
    /// \param element the element
    /// \param newElementInitialization the new element initialization
    template <class Container, class T>
    void AlignContainerElements(Container &elements, const T &element, const T &newElementInitialization);

    template <typename T>
    void InitializeNuberVectorWithConstantElements(std::vector<unsigned int> &numberElementVector,
                                                   std::vector<T> &elementVector,
                                                   const unsigned int numberElementSize,
                                                   const unsigned int numberElements,
                                                   const T &elementInitialization = T());
    template <typename T>
    void InitializeNumberVector(std::vector<unsigned int> &numberElementVector,
                                std::vector<T> &elementVector,
                                const std::vector<unsigned int> &numberElements,
                                const T &elementInitialization = T());

    template <typename T>
    void ResizeNumberVectorWithNewNumberElements(std::vector<unsigned int> &numberElementVector,
                                                 std::vector<T> &elementVector,
                                                 const unsigned int &numberElements,
                                                 const unsigned int &vectorIndex,
                                                 const unsigned int &newNumberElements,
                                                 const T &newElementInitialization = T());

  public:
    MeshMatricesDAO(Gedim::MeshMatrices &mesh);
    ~MeshMatricesDAO();

    inline Gedim::MeshMatrices &MeshData()
    {
        return _mesh;
    }
    inline const Gedim::MeshMatrices &MeshData() const
    {
        return _mesh;
    }

    inline void InitializeDimension(const unsigned int &dimension)
    {
        _mesh.Dimension = dimension;
    }
    inline unsigned int Dimension() const
    {
        return _mesh.Dimension;
    }

    void Cell0DsInitialize(const unsigned int &numberCell0Ds);
    unsigned int Cell0DAppend(const unsigned int &numberCell0Ds);
    void Cell0DRemove(const unsigned int &cell0DIndex);

    void Cell0DInsertCoordinates(const unsigned int &cell0DIndex, const Eigen::Vector3d &coordinates);
    void Cell0DsInsertCoordinates(const Eigen::MatrixXd &coordinates);
    inline void Cell0DSetMarker(const unsigned int &cell0DIndex, const unsigned int &marker)
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        _mesh.Cell0DMarkers[cell0DIndex] = marker;
    }
    inline void Cell0DSetState(const unsigned int &cell0DIndex, const bool &state)
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        _mesh.ActiveCell0D[cell0DIndex] = state;
    }

    inline unsigned int Cell0DTotalNumber() const
    {
        return _mesh.NumberCell0D;
    }
    inline double Cell0DCoordinateX(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DCoordinates[3 * cell0DIndex];
    }
    inline double Cell0DCoordinateY(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DCoordinates[3 * cell0DIndex + 1];
    }
    inline double Cell0DCoordinateZ(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DCoordinates[3 * cell0DIndex + 2];
    }
    inline Eigen::Vector3d Cell0DCoordinates(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return Eigen::Vector3d(Cell0DCoordinateX(cell0DIndex), Cell0DCoordinateY(cell0DIndex), Cell0DCoordinateZ(cell0DIndex));
    }
    Eigen::MatrixXd Cell0DsCoordinates() const;
    Eigen::MatrixXd Cell0DsCoordinates(const std::vector<unsigned int> &cell0Ds) const;
    inline unsigned int Cell0DMarker(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.Cell0DMarkers[cell0DIndex];
    }
    inline std::vector<unsigned int> Cell0DsMarker() const
    {
        return _mesh.Cell0DMarkers;
    }
    inline bool Cell0DIsActive(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.ActiveCell0D[cell0DIndex];
    }
    inline std::vector<bool> Cell0DsState() const
    {
        return _mesh.ActiveCell0D;
    }

    inline bool Cell0DHasUpdatedCell0Ds(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.UpdatedCell0Ds.find(cell0DIndex) != _mesh.UpdatedCell0Ds.end();
    }
    inline unsigned int Cell0DNumberUpdatedCell0Ds(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.UpdatedCell0Ds.at(cell0DIndex).size();
    }
    inline bool Cell0DHasUpdatedCell0D(const unsigned int &cell0DIndex, const unsigned int &updatedCell0DIdex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(updatedCell0DIdex < Cell0DTotalNumber());
        return _mesh.UpdatedCell0Ds.at(cell0DIndex).find(updatedCell0DIdex) != _mesh.UpdatedCell0Ds.at(cell0DIndex).end();
    }
    void Cell0DInsertUpdatedCell0D(const unsigned int &cell0DIndex, const unsigned int &updatedCell0DIdex);
    bool Cell0DUpdatedCell0Ds(const unsigned int &cell0DIndex, std::list<unsigned int> &updatedCell0DIds) const;

    std::vector<std::vector<unsigned int>> Cell0DsNeighbourCell1Ds() const;
    inline void Cell0DsInitializeNeighbourCell1Ds(const std::vector<unsigned int> &numberNeighbourCell1Ds);
    inline void Cell0DInitializeNeighbourCell1Ds(const unsigned int &cell0DIndex, const unsigned int &numberNeighbourCell1Ds);
    void Cell0DInitializeNeighbourCell1Ds(const unsigned int &cell0DIndex, const std::vector<unsigned int> &neighbourCell1Ds);
    inline void Cell0DInsertNeighbourCell1D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex, const unsigned int &neigbourCell1DIndex)
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell1D(cell0DIndex));
        Gedim::Output::Assert(neigbourCell1DIndex < Cell1DTotalNumber());

        _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] + neighbourIndex] = neigbourCell1DIndex;
    }
    inline unsigned int Cell0DNumberNeighbourCell1D(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.NumberCell0DNeighbourCell1D[cell0DIndex + 1] - _mesh.NumberCell0DNeighbourCell1D[cell0DIndex];
    }
    inline unsigned int Cell0DNeighbourCell1D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell1D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] + neighbourIndex];
    }
    inline std::vector<unsigned int> Cell0DNeighbourCell1Ds(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());

        std::list<unsigned int> neighbours;
        for (unsigned int n = 0; n < Cell0DNumberNeighbourCell1D(cell0DIndex); n++)
        {
            if (!Cell0DHasNeighbourCell1D(cell0DIndex, n))
                continue;

            neighbours.push_back(Cell0DNeighbourCell1D(cell0DIndex, n));
        }

        return std::vector<unsigned int>(neighbours.begin(), neighbours.end());
    }
    inline bool Cell0DHasNeighbourCell1D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell1D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] + neighbourIndex] < _mesh.NumberCell1D;
    }
    inline void Cell0DResetNeighbourCell1D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex)
    {
        _mesh.Cell0DNeighbourCell1Ds[_mesh.NumberCell0DNeighbourCell1D[cell0DIndex] + neighbourIndex] =
            std::numeric_limits<unsigned int>::max();
    }

    std::vector<std::vector<unsigned int>> Cell0DsNeighbourCell2Ds() const;
    inline void Cell0DsInitializeNeighbourCell2Ds(const std::vector<unsigned int> &numberNeighbourCell2Ds);
    inline void Cell0DInitializeNeighbourCell2Ds(const unsigned int &cell0DIndex, const unsigned int &numberNeighbourCell2Ds);
    void Cell0DInitializeNeighbourCell2Ds(const unsigned int &cell0DIndex, const std::vector<unsigned int> &neighbourCell2Ds);
    inline void Cell0DInsertNeighbourCell2D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex, const unsigned int &neigbourCell2DIndex)
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell2D(cell0DIndex));
        Gedim::Output::Assert(neigbourCell2DIndex < Cell2DTotalNumber());

        _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] + neighbourIndex] = neigbourCell2DIndex;
    }
    inline unsigned int Cell0DNumberNeighbourCell2D(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.NumberCell0DNeighbourCell2D[cell0DIndex + 1] - _mesh.NumberCell0DNeighbourCell2D[cell0DIndex];
    }
    inline unsigned int Cell0DNeighbourCell2D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell2D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] + neighbourIndex];
    }
    inline std::vector<unsigned int> Cell0DNeighbourCell2Ds(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());

        std::list<unsigned int> neighbours;
        for (unsigned int n = 0; n < Cell0DNumberNeighbourCell2D(cell0DIndex); n++)
        {
            if (!Cell0DHasNeighbourCell2D(cell0DIndex, n))
                continue;

            neighbours.push_back(Cell0DNeighbourCell2D(cell0DIndex, n));
        }

        return std::vector<unsigned int>(neighbours.begin(), neighbours.end());
    }
    inline bool Cell0DHasNeighbourCell2D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell2D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] + neighbourIndex] < _mesh.NumberCell2D;
    }
    inline void Cell0DResetNeighbourCell2D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex)
    {
        _mesh.Cell0DNeighbourCell2Ds[_mesh.NumberCell0DNeighbourCell2D[cell0DIndex] + neighbourIndex] =
            std::numeric_limits<unsigned int>::max();
    }
    std::vector<std::vector<unsigned int>> Cell0DsNeighbourCell3Ds() const;
    void Cell0DsInitializeNeighbourCell3Ds(const std::vector<unsigned int> &numberNeighbourCell3Ds);
    void Cell0DInitializeNeighbourCell3Ds(const unsigned int &cell0DIndex, const unsigned int &numberNeighbourCell3Ds);
    void Cell0DInitializeNeighbourCell3Ds(const unsigned int &cell0DIndex, const std::vector<unsigned int> &neighbourCell3Ds);
    inline void Cell0DInsertNeighbourCell3D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex, const unsigned int &neigbourCell3DIndex)
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell3D(cell0DIndex));
        Gedim::Output::Assert(neigbourCell3DIndex < Cell3DTotalNumber());

        _mesh.Cell0DNeighbourCell3Ds[_mesh.NumberCell0DNeighbourCell3D[cell0DIndex] + neighbourIndex] = neigbourCell3DIndex;
    }
    inline unsigned int Cell0DNumberNeighbourCell3D(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        return _mesh.NumberCell0DNeighbourCell3D[cell0DIndex + 1] - _mesh.NumberCell0DNeighbourCell3D[cell0DIndex];
    }
    inline unsigned int Cell0DNeighbourCell3D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell3D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell3Ds[_mesh.NumberCell0DNeighbourCell3D[cell0DIndex] + neighbourIndex];
    }
    inline std::vector<unsigned int> Cell0DNeighbourCell3Ds(const unsigned int &cell0DIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());

        std::list<unsigned int> neighbours;
        for (unsigned int n = 0; n < Cell0DNumberNeighbourCell3D(cell0DIndex); n++)
        {
            if (!Cell0DHasNeighbourCell3D(cell0DIndex, n))
                continue;

            neighbours.push_back(Cell0DNeighbourCell3D(cell0DIndex, n));
        }

        return std::vector<unsigned int>(neighbours.begin(), neighbours.end());
    }
    inline bool Cell0DHasNeighbourCell3D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell0DNumberNeighbourCell3D(cell0DIndex));
        return _mesh.Cell0DNeighbourCell3Ds[_mesh.NumberCell0DNeighbourCell3D[cell0DIndex] + neighbourIndex] < _mesh.NumberCell3D;
    }
    inline void Cell0DResetNeighbourCell3D(const unsigned int &cell0DIndex, const unsigned int &neighbourIndex)
    {
        _mesh.Cell0DNeighbourCell3Ds[_mesh.NumberCell0DNeighbourCell3D[cell0DIndex] + neighbourIndex] =
            std::numeric_limits<unsigned int>::max();
    }

    void Cell0DInitializeDoubleProperties(const unsigned int &numberDoubleProperties);
    unsigned int Cell0DAddDoubleProperty(const std::string &propertyId);
    inline void Cell0DsInitializeDoublePropertyValues(const unsigned int &propertyIndex, const std::vector<unsigned int> &propertySizes);
    void Cell0DInitializeDoublePropertyValues(const unsigned int &cell0DIndex,
                                              const unsigned int &propertyIndex,
                                              const unsigned int &propertySize);
    inline void Cell0DInsertDoublePropertyValue(const unsigned int &cell0DIndex,
                                                const unsigned int &propertyIndex,
                                                const unsigned int &propertyValueIndex,
                                                const double &propertyValue)
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());

        _mesh.Cell0DDoublePropertyValues[propertyIndex][_mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex] + propertyValueIndex] =
            propertyValue;
    }
    inline unsigned int Cell0DNumberDoubleProperties() const
    {
        return _mesh.Cell0DDoublePropertyIds.size();
    }
    inline std::string Cell0DDoublePropertyId(const unsigned int &propertyIndex) const
    {
        return _mesh.Cell0DDoublePropertyIds[propertyIndex];
    }
    inline bool Cell0DDoublePropertyExists(const std::string &propertyId) const
    {
        return _mesh.Cell0DDoublePropertyIndices.find(propertyId) != _mesh.Cell0DDoublePropertyIndices.end();
    }
    inline unsigned int Cell0DDoublePropertyIndex(const std::string &propertyId) const
    {
        Gedim::Output::Assert(Cell0DDoublePropertyExists(propertyId));
        return _mesh.Cell0DDoublePropertyIndices.at(propertyId);
    }
    inline unsigned int Cell0DDoublePropertySize(const unsigned int &cell0DIndex, const unsigned int &propertyIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());
        return _mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex + 1] -
               _mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex];
    }
    inline double Cell0DDoublePropertyValue(const unsigned int &cell0DIndex,
                                            const unsigned int &propertyIndex,
                                            const unsigned int &propertyValueIndex) const
    {
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell0DNumberDoubleProperties());
        Gedim::Output::Assert(propertyValueIndex < Cell0DDoublePropertySize(cell0DIndex, propertyIndex));

        return _mesh.Cell0DDoublePropertyValues[propertyIndex][_mesh.Cell0DDoublePropertySizes[propertyIndex][cell0DIndex] + propertyValueIndex];
    }

    void Cell1DsInitialize(const unsigned int &numberCell1Ds);
    unsigned int Cell1DAppend(const unsigned int &numberCell1Ds);
    void Cell1DRemove(const unsigned int &cell1DIndex);
    void Cell1DInsertExtremes(const unsigned int &cell1DIndex, const unsigned int &originCell0DIndex, const unsigned int &endCell0DIndex);

    void Cell1DsInsertExtremes(const Eigen::MatrixXi &cell1DExtremes);

    Eigen::MatrixXi Cell1DsExtremes() const;
    Eigen::MatrixXi Cell1DsExtremes(const std::vector<unsigned int> &cell1Ds) const;
    /// \return the extrems as Eigen MatrixXi of cell1D, size 2
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    inline Eigen::VectorXi Cell1DExtremes(const unsigned int &cell1DIndex) const
    {
        return (Eigen::VectorXi(2) << Cell1DOrigin(cell1DIndex), Cell1DEnd(cell1DIndex)).finished();
    }

    unsigned int Cell1DByExtremes(const unsigned int &originCell0DIndex, const unsigned int &endCell0DIndex) const;

    inline void Cell1DSetMarker(const unsigned int &cell1DIndex, const unsigned int &marker)
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        _mesh.Cell1DMarkers[cell1DIndex] = marker;
    }
    inline void Cell1DSetState(const unsigned int &cell1DIndex, const bool &state)
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        _mesh.ActiveCell1D[cell1DIndex] = state;
    }
    std::vector<std::vector<unsigned int>> Cell1DsNeighbourCell2Ds() const;
    inline void Cell1DsInitializeNeighbourCell2Ds(const std::vector<unsigned int> &numberNeighbourCell2Ds);
    inline void Cell1DsInitializeNeighbourCell2Ds(const unsigned int &numberNeighbourCell2Ds);
    void Cell1DInitializeNeighbourCell2Ds(const unsigned int &cell1DIndex, const unsigned int &numberNeighbourCell2Ds);
    void Cell1DInitializeNeighbourCell2Ds(const unsigned int &cell1DIndex, const std::vector<unsigned int> &neighbourCell2Ds);
    inline void Cell1DInsertNeighbourCell2D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex, const unsigned int &neigbourCell2DIndex)
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell2D(cell1DIndex));
        Gedim::Output::Assert(neigbourCell2DIndex < Cell2DTotalNumber());

        _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] + neighbourIndex] = neigbourCell2DIndex;
    }
    inline unsigned int Cell1DTotalNumber() const
    {
        return _mesh.NumberCell1D;
    }
    inline unsigned int Cell1DVertex(const unsigned int &cell1DIndex, const unsigned int &vertexIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(vertexIndex < 2);
        return _mesh.Cell1DVertices[2 * cell1DIndex + vertexIndex];
    }
    inline unsigned int Cell1DOrigin(const unsigned int &cell1DIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DVertices[2 * cell1DIndex];
    }
    inline unsigned int Cell1DEnd(const unsigned int &cell1DIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DVertices[2 * cell1DIndex + 1];
    }
    inline unsigned int Cell1DFindExtreme(const unsigned int &cell1DIndex, const unsigned int &cell0DIndex) const
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
    inline Eigen::MatrixXd Cell1DCoordinates(const unsigned int &cell1DIndex) const
    {
        return (Eigen::MatrixXd(3, 2) << Cell1DOriginCoordinates(cell1DIndex), Cell1DEndCoordinates(cell1DIndex)).finished();
    }
    inline Eigen::Vector3d Cell1DOriginCoordinates(const unsigned int &cell1DIndex) const
    {
        return Cell0DCoordinates(Cell1DOrigin(cell1DIndex));
    }
    inline Eigen::Vector3d Cell1DEndCoordinates(const unsigned int &cell1DIndex) const
    {
        return Cell0DCoordinates(Cell1DEnd(cell1DIndex));
    }
    inline unsigned int Cell1DNumberNeighbourCell2D(const unsigned int &cell1DIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.NumberCell1DNeighbourCell2D[cell1DIndex + 1] - _mesh.NumberCell1DNeighbourCell2D[cell1DIndex];
    }
    inline unsigned int Cell1DNeighbourCell2D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell2D(cell1DIndex));
        return _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] + neighbourIndex];
    }
    inline std::vector<unsigned int> Cell1DNeighbourCell2Ds(const unsigned int &cell1DIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());

        std::list<unsigned int> neighbours;
        for (unsigned int n = 0; n < Cell1DNumberNeighbourCell2D(cell1DIndex); n++)
        {
            if (!Cell1DHasNeighbourCell2D(cell1DIndex, n))
                continue;

            neighbours.push_back(Cell1DNeighbourCell2D(cell1DIndex, n));
        }

        return std::vector<unsigned int>(neighbours.begin(), neighbours.end());
    }
    inline bool Cell1DHasNeighbourCell2D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell2D(cell1DIndex));
        return _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] + neighbourIndex] < _mesh.NumberCell2D;
    }
    inline void Cell1DResetNeighbourCell2D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex)
    {
        _mesh.Cell1DNeighbourCell2Ds[_mesh.NumberCell1DNeighbourCell2D[cell1DIndex] + neighbourIndex] =
            std::numeric_limits<unsigned int>::max();
    }

    inline unsigned int Cell1DMarker(const unsigned int &cell1DIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DMarkers[cell1DIndex];
    }
    inline std::vector<unsigned int> Cell1DsMarker() const
    {
        return _mesh.Cell1DMarkers;
    }
    inline bool Cell1DIsActive(const unsigned int &cell1DIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.ActiveCell1D[cell1DIndex];
    }
    inline std::vector<bool> Cell1DsState() const
    {
        return _mesh.ActiveCell1D;
    }
    inline bool Cell1DHasOriginalCell1D(const unsigned int &updatedCell1DIndex) const
    {
        Gedim::Output::Assert(updatedCell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DOriginalCell1Ds.at(updatedCell1DIndex) < _mesh.NumberCell1D;
    }
    inline unsigned int Cell1DOriginalCell1D(const unsigned int &updatedCell1DIndex) const
    {
        Gedim::Output::Assert(updatedCell1DIndex < Cell1DTotalNumber());
        return _mesh.Cell1DOriginalCell1Ds.at(updatedCell1DIndex);
    }
    inline bool Cell1DHasUpdatedCell1Ds(const unsigned int &cell1DIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.UpdatedCell1Ds.find(cell1DIndex) != _mesh.UpdatedCell1Ds.end();
    }
    inline unsigned int Cell1DNumberUpdatedCell1Ds(const unsigned int &cell1DIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.UpdatedCell1Ds.at(cell1DIndex).size();
    }
    inline bool Cell1DHasUpdatedCell1D(const unsigned int &cell1DIndex, const unsigned int &updatedCell1DIdex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(updatedCell1DIdex < Cell1DTotalNumber());
        return _mesh.UpdatedCell1Ds.at(cell1DIndex).find(updatedCell1DIdex) != _mesh.UpdatedCell1Ds.at(cell1DIndex).end();
    }
    void Cell1DInsertUpdatedCell1D(const unsigned int &cell1DIndex, const unsigned int &updatedCell1DIdex);
    bool Cell1DUpdatedCell1Ds(const unsigned int &cell1DIndex, std::list<unsigned int> &updatedCell1DIds) const;
    void Cell1DInitializeDoubleProperties(const unsigned int &numberDoubleProperties);

    std::vector<std::vector<unsigned int>> Cell1DsNeighbourCell3Ds() const;
    void Cell1DsInitializeNeighbourCell3Ds(const std::vector<unsigned int> &numberNeighbourCell3Ds);
    void Cell1DInitializeNeighbourCell3Ds(const unsigned int &cell1DIndex, const unsigned int &numberNeighbourCell3Ds);
    void Cell1DInitializeNeighbourCell3Ds(const unsigned int &cell1DIndex, const std::vector<unsigned int> &neighbourCell3Ds);
    inline void Cell1DInsertNeighbourCell3D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex, const unsigned int &neigbourCell3DIndex)
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell3D(cell1DIndex));
        Gedim::Output::Assert(neigbourCell3DIndex < Cell3DTotalNumber());

        _mesh.Cell1DNeighbourCell3Ds[_mesh.NumberCell1DNeighbourCell3D[cell1DIndex] + neighbourIndex] = neigbourCell3DIndex;
    }
    inline unsigned int Cell1DNumberNeighbourCell3D(const unsigned int &cell1DIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        return _mesh.NumberCell1DNeighbourCell3D[cell1DIndex + 1] - _mesh.NumberCell1DNeighbourCell3D[cell1DIndex];
    }
    inline unsigned int Cell1DNeighbourCell3D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell3D(cell1DIndex));
        return _mesh.Cell1DNeighbourCell3Ds[_mesh.NumberCell1DNeighbourCell3D[cell1DIndex] + neighbourIndex];
    }
    inline std::vector<unsigned int> Cell1DNeighbourCell3Ds(const unsigned int &cell1DIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());

        std::list<unsigned int> neighbours;
        for (unsigned int n = 0; n < Cell1DNumberNeighbourCell3D(cell1DIndex); n++)
        {
            if (!Cell1DHasNeighbourCell3D(cell1DIndex, n))
                continue;

            neighbours.push_back(Cell1DNeighbourCell3D(cell1DIndex, n));
        }

        return std::vector<unsigned int>(neighbours.begin(), neighbours.end());
    }
    inline bool Cell1DHasNeighbourCell3D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell1DNumberNeighbourCell3D(cell1DIndex));
        return _mesh.Cell1DNeighbourCell3Ds[_mesh.NumberCell1DNeighbourCell3D[cell1DIndex] + neighbourIndex] < _mesh.NumberCell3D;
    }
    inline void Cell1DResetNeighbourCell3D(const unsigned int &cell1DIndex, const unsigned int &neighbourIndex)
    {
        _mesh.Cell1DNeighbourCell3Ds[_mesh.NumberCell1DNeighbourCell3D[cell1DIndex] + neighbourIndex] =
            std::numeric_limits<unsigned int>::max();
    }

    unsigned int Cell1DAddDoubleProperty(const std::string &propertyId);
    inline void Cell1DsInitializeDoublePropertyValues(const unsigned int &propertyIndex, const std::vector<unsigned int> &propertySizes);
    void Cell1DInitializeDoublePropertyValues(const unsigned int &cell1DIndex,
                                              const unsigned int &propertyIndex,
                                              const unsigned int &propertySize);
    inline void Cell1DInsertDoublePropertyValue(const unsigned int &cell1DIndex,
                                                const unsigned int &propertyIndex,
                                                const unsigned int &propertyValueIndex,
                                                const double &propertyValue)
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());

        _mesh.Cell1DDoublePropertyValues[propertyIndex][_mesh.Cell1DDoublePropertySizes[propertyIndex][cell1DIndex] + propertyValueIndex] =
            propertyValue;
    }
    inline unsigned int Cell1DNumberDoubleProperties() const
    {
        return _mesh.Cell1DDoublePropertyIds.size();
    }
    inline std::string Cell1DDoublePropertyId(const unsigned int &propertyIndex) const
    {
        return _mesh.Cell1DDoublePropertyIds[propertyIndex];
    }
    inline bool Cell1DDoublePropertyExists(const std::string &propertyId) const
    {
        return _mesh.Cell1DDoublePropertyIndices.find(propertyId) != _mesh.Cell1DDoublePropertyIndices.end();
    }
    inline unsigned int Cell1DDoublePropertyIndex(const std::string &propertyId) const
    {
        Gedim::Output::Assert(Cell1DDoublePropertyExists(propertyId));
        return _mesh.Cell1DDoublePropertyIndices.at(propertyId);
    }
    inline unsigned int Cell1DDoublePropertySize(const unsigned int &cell1DIndex, const unsigned int &propertyIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());
        return _mesh.Cell1DDoublePropertySizes[propertyIndex][cell1DIndex + 1] -
               _mesh.Cell1DDoublePropertySizes[propertyIndex][cell1DIndex];
    }
    inline double Cell1DDoublePropertyValue(const unsigned int &cell1DIndex,
                                            const unsigned int &propertyIndex,
                                            const unsigned int &propertyValueIndex) const
    {
        Gedim::Output::Assert(cell1DIndex < Cell1DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell1DNumberDoubleProperties());
        Gedim::Output::Assert(propertyValueIndex < Cell1DDoublePropertySize(cell1DIndex, propertyIndex));

        return _mesh.Cell1DDoublePropertyValues[propertyIndex][_mesh.Cell1DDoublePropertySizes[propertyIndex][cell1DIndex] + propertyValueIndex];
    }

    void Cell2DsInitialize(const unsigned int &numberCell2Ds);
    unsigned int Cell2DAppend(const unsigned int &numberCell2Ds);
    void Cell2DRemove(const unsigned int &cell2DIndex);
    inline void Cell2DSetMarker(const unsigned int &cell2DIndex, const unsigned int &marker)
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        _mesh.Cell2DMarkers[cell2DIndex] = marker;
    }
    inline void Cell2DSetState(const unsigned int &cell2DIndex, const bool &state)
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        _mesh.ActiveCell2D[cell2DIndex] = state;
    }
    inline void Cell2DsInitializeVertices(const unsigned int &numberCell2DVertices);
    inline void Cell2DsInitializeVertices(const std::vector<unsigned int> &numberCell2DsVertices);
    void Cell2DInitializeVertices(const unsigned int &cell2DIndex, const unsigned int &numberCell2DVertices);
    inline void Cell2DsInitializeEdges(const unsigned int &numberCell2DEdges);
    inline void Cell2DsInitializeEdges(const std::vector<unsigned int> &numberCell2DsEdges);
    void Cell2DInitializeEdges(const unsigned int &cell2DIndex, const unsigned int &numberCell2DEdges);
    inline void Cell2DInsertVertices(const unsigned int &cell2DIndex, const std::vector<unsigned int> &verticesCell0DIndices)
    {
        for (unsigned int v = 0; v < verticesCell0DIndices.size(); v++)
            Cell2DInsertVertex(cell2DIndex, v, verticesCell0DIndices[v]);
    }
    inline void Cell2DInsertVertex(const unsigned int &cell2DIndex, const unsigned int &vertexIndex, const unsigned int &vertexCell0DIndex)
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell2DNumberVertices(cell2DIndex));
        Gedim::Output::Assert(vertexCell0DIndex < Cell0DTotalNumber());
        _mesh.Cell2DVertices[_mesh.NumberCell2DVertices[cell2DIndex] + vertexIndex] = vertexCell0DIndex;
    }
    inline void Cell2DInsertEdges(const unsigned int &cell2DIndex, const std::vector<unsigned int> &edgesCell1DIndices)
    {
        for (unsigned int e = 0; e < edgesCell1DIndices.size(); e++)
            Cell2DInsertEdge(cell2DIndex, e, edgesCell1DIndices[e]);
    }
    inline void Cell2DInsertEdge(const unsigned int &cell2DIndex, const unsigned int &edgeIndex, const unsigned int &edgeCell1DIndex)
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(edgeIndex < Cell2DNumberEdges(cell2DIndex));
        Gedim::Output::Assert(edgeCell1DIndex < Cell1DTotalNumber());
        _mesh.Cell2DEdges[_mesh.NumberCell2DEdges[cell2DIndex] + edgeIndex] = edgeCell1DIndex;
    }
    void Cell2DAddVertices(const unsigned int &cell2DIndex, const std::vector<unsigned int> &verticesCell0DIndices);
    void Cell2DAddEdges(const unsigned int &cell2DIndex, const std::vector<unsigned int> &edgesCell1DIndices);
    void Cell2DAddVerticesAndEdges(const unsigned int &cell2DIndex, const Eigen::MatrixXi &verticesAndEdgesIndices);

    inline unsigned int Cell2DTotalNumber() const
    {
        return _mesh.NumberCell2D;
    }
    inline unsigned int Cell2DNumberVertices(const unsigned int &cell2DIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.NumberCell2DVertices[cell2DIndex + 1] - _mesh.NumberCell2DVertices[cell2DIndex];
    }
    inline unsigned int Cell2DNumberEdges(const unsigned int &cell2DIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.NumberCell2DEdges[cell2DIndex + 1] - _mesh.NumberCell2DEdges[cell2DIndex];
    }

    inline std::vector<unsigned int> Cell2DVertices(const unsigned int &cell2DIndex) const
    {
        return std::vector<unsigned int>(_mesh.Cell2DVertices.begin() + _mesh.NumberCell2DVertices[cell2DIndex],
                                         _mesh.Cell2DVertices.begin() + _mesh.NumberCell2DVertices[cell2DIndex] +
                                             Cell2DNumberVertices(cell2DIndex));
    }
    std::vector<std::vector<unsigned int>> Cell2DsVertices() const;
    std::vector<Eigen::MatrixXi> Cell2DsExtremes() const;

    inline unsigned int Cell2DVertex(const unsigned int &cell2DIndex, const unsigned int &vertexIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell2DNumberVertices(cell2DIndex));
        return _mesh.Cell2DVertices[_mesh.NumberCell2DVertices[cell2DIndex] + vertexIndex];
    }
    inline Eigen::Vector3d Cell2DVertexCoordinates(const unsigned int &cell2DIndex, const unsigned int &vertexIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell2DNumberVertices(cell2DIndex));
        return Cell0DCoordinates(Cell2DVertex(cell2DIndex, vertexIndex));
    }
    Eigen::MatrixXd Cell2DVerticesCoordinates(const unsigned int &cell2DIndex) const;
    unsigned int Cell2DFindVertex(const unsigned int &cell2DIndex, const unsigned int &cell0DIndex) const;

    inline std::vector<unsigned int> Cell2DEdges(const unsigned int &cell2DIndex) const
    {
        return std::vector<unsigned int>(_mesh.Cell2DEdges.begin() + _mesh.NumberCell2DEdges[cell2DIndex],
                                         _mesh.Cell2DEdges.begin() + _mesh.NumberCell2DEdges[cell2DIndex] +
                                             Cell2DNumberEdges(cell2DIndex));
    }

    inline unsigned int Cell2DEdge(const unsigned int &cell2DIndex, const unsigned int &edgeIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(edgeIndex < Cell2DNumberEdges(cell2DIndex));
        return _mesh.Cell2DEdges[_mesh.NumberCell2DEdges[cell2DIndex] + edgeIndex];
    }
    unsigned int Cell2DFindEdge(const unsigned int &cell2DIndex, const unsigned int &cell1DIndex) const;
    unsigned int Cell2DFindEdgeByExtremes(const unsigned int &cell2DIndex,
                                          const unsigned int &originCell0DIndex,
                                          const unsigned int &endCell0DIndex) const;
    inline unsigned int Cell2DMarker(const unsigned int &cell2DIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.Cell2DMarkers[cell2DIndex];
    }
    inline std::vector<unsigned int> Cell2DsMarker() const
    {
        return _mesh.Cell2DMarkers;
    }
    inline bool Cell2DIsActive(const unsigned int &cell2DIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.ActiveCell2D[cell2DIndex];
    }
    inline std::vector<bool> Cell2DsState() const
    {
        return _mesh.ActiveCell2D;
    }

    inline bool Cell2DHasUpdatedCell2Ds(const unsigned int &cell2DIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.UpdatedCell2Ds.find(cell2DIndex) != _mesh.UpdatedCell2Ds.end();
    }
    inline unsigned int Cell2DNumberUpdatedCell2Ds(const unsigned int &cell2DIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.UpdatedCell2Ds.at(cell2DIndex).size();
    }
    inline bool Cell2DHasUpdatedCell2D(const unsigned int &cell2DIndex, const unsigned int &updatedCell2DIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(updatedCell2DIndex < Cell2DTotalNumber());
        return _mesh.UpdatedCell2Ds.at(cell2DIndex).find(updatedCell2DIndex) != _mesh.UpdatedCell2Ds.at(cell2DIndex).end();
    }
    void Cell2DInsertUpdatedCell2D(const unsigned int &cell2DIndex, const unsigned int &updatedCell2DIdex);
    inline bool Cell2DHasOriginalCell2D(const unsigned int &updatedCell2DIndex) const
    {
        Gedim::Output::Assert(updatedCell2DIndex < Cell2DTotalNumber());
        return _mesh.Cell2DOriginalCell2Ds.at(updatedCell2DIndex) < _mesh.NumberCell2D;
    }
    inline unsigned int Cell2DOriginalCell2D(const unsigned int &updatedCell2DIndex) const
    {
        Gedim::Output::Assert(updatedCell2DIndex < Cell2DTotalNumber());
        return _mesh.Cell2DOriginalCell2Ds.at(updatedCell2DIndex);
    }
    bool Cell2DUpdatedCell2Ds(const unsigned int &cell2DIndex, std::list<unsigned int> &updatedCell2DIds) const;

    std::vector<std::vector<unsigned int>> Cell2DsNeighbourCell3Ds() const;
    inline void Cell2DsInitializeNeighbourCell3Ds(const std::vector<unsigned int> &numberNeighbourCell3Ds);
    void Cell2DInitializeNeighbourCell3Ds(const unsigned int &cell2DIndex, const unsigned int &numberNeighbourCell3Ds);
    void Cell2DInitializeNeighbourCell3Ds(const unsigned int &cell2DIndex, const std::vector<unsigned int> &neighbourCell3Ds);
    inline void Cell2DInsertNeighbourCell3D(const unsigned int &cell2DIndex, const unsigned int &neighbourIndex, const unsigned int &neigbourCell3DIndex)
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell2DNumberNeighbourCell3D(cell2DIndex));
        Gedim::Output::Assert(neigbourCell3DIndex < Cell3DTotalNumber());

        _mesh.Cell2DNeighbourCell3Ds[_mesh.NumberCell2DNeighbourCell3D[cell2DIndex] + neighbourIndex] = neigbourCell3DIndex;
    }
    inline unsigned int Cell2DNumberNeighbourCell3D(const unsigned int &cell2DIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.NumberCell2DNeighbourCell3D[cell2DIndex + 1] - _mesh.NumberCell2DNeighbourCell3D[cell2DIndex];
    }
    inline unsigned int Cell2DNeighbourCell3D(const unsigned int &cell2DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell2DNumberNeighbourCell3D(cell2DIndex));
        return _mesh.Cell2DNeighbourCell3Ds[_mesh.NumberCell2DNeighbourCell3D[cell2DIndex] + neighbourIndex];
    }
    inline std::vector<unsigned int> Cell2DNeighbourCell3Ds(const unsigned int &cell2DIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());

        std::list<unsigned int> neighbours;
        for (unsigned int n = 0; n < Cell2DNumberNeighbourCell3D(cell2DIndex); n++)
        {
            if (!Cell2DHasNeighbourCell3D(cell2DIndex, n))
                continue;

            neighbours.push_back(Cell2DNeighbourCell3D(cell2DIndex, n));
        }

        return std::vector<unsigned int>(neighbours.begin(), neighbours.end());
    }
    inline bool Cell2DHasNeighbourCell3D(const unsigned int &cell2DIndex, const unsigned int &neighbourIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(neighbourIndex < Cell2DNumberNeighbourCell3D(cell2DIndex));
        return _mesh.Cell2DNeighbourCell3Ds[_mesh.NumberCell2DNeighbourCell3D[cell2DIndex] + neighbourIndex] < _mesh.NumberCell3D;
    }
    inline void Cell2DResetNeighbourCell3D(const unsigned int &cell2DIndex, const unsigned int &neighbourIndex)
    {
        _mesh.Cell2DNeighbourCell3Ds[_mesh.NumberCell2DNeighbourCell3D[cell2DIndex] + neighbourIndex] =
            std::numeric_limits<unsigned int>::max();
    }
    void Cell2DInitializeDoubleProperties(const unsigned int &numberDoubleProperties);
    unsigned int Cell2DAddDoubleProperty(const std::string &propertyId);
    inline void Cell2DsInitializeDoublePropertyValues(const unsigned int &propertyIndex, const std::vector<unsigned int> &propertySizes);
    void Cell2DInitializeDoublePropertyValues(const unsigned int &cell2DIndex,
                                              const unsigned int &propertyIndex,
                                              const unsigned int &propertySize);
    inline void Cell2DInsertDoublePropertyValue(const unsigned int &cell2DIndex,
                                                const unsigned int &propertyIndex,
                                                const unsigned int &propertyValueIndex,
                                                const double &propertyValue)
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());

        _mesh.Cell2DDoublePropertyValues[propertyIndex][_mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex] + propertyValueIndex] =
            propertyValue;
    }
    inline unsigned int Cell2DNumberDoubleProperties() const
    {
        return _mesh.Cell2DDoublePropertyIds.size();
    }
    inline std::string Cell2DDoublePropertyId(const unsigned int &propertyIndex) const
    {
        return _mesh.Cell2DDoublePropertyIds[propertyIndex];
    }
    inline bool Cell2DDoublePropertyExists(const std::string &propertyId) const
    {
        return _mesh.Cell2DDoublePropertyIndices.find(propertyId) != _mesh.Cell2DDoublePropertyIndices.end();
    }
    inline unsigned int Cell2DDoublePropertyIndex(const std::string &propertyId) const
    {
        Gedim::Output::Assert(Cell2DDoublePropertyExists(propertyId));
        return _mesh.Cell2DDoublePropertyIndices.at(propertyId);
    }
    inline unsigned int Cell2DDoublePropertySize(const unsigned int &cell2DIndex, const unsigned int &propertyIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());
        return _mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex + 1] -
               _mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex];
    }
    inline double Cell2DDoublePropertyValue(const unsigned int &cell2DIndex,
                                            const unsigned int &propertyIndex,
                                            const unsigned int &propertyValueIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell2DNumberDoubleProperties());
        Gedim::Output::Assert(propertyValueIndex < Cell2DDoublePropertySize(cell2DIndex, propertyIndex));

        return _mesh.Cell2DDoublePropertyValues[propertyIndex][_mesh.Cell2DDoublePropertySizes[propertyIndex][cell2DIndex] + propertyValueIndex];
    }

    inline void Cell2DsInitializeSubDivision(const std::vector<unsigned int> &numberSubDivisions);

    void Cell2DInitializeSubDivision(const unsigned int &cell2DIndex, const unsigned int &numberSubDivision);
    inline void Cell2DInsertSubDivision(const unsigned int &cell2DIndex, const unsigned int &subDivisionIndex, const unsigned int &cell0DIndex)
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(subDivisionIndex < Cell2DNumberSubDivision(cell2DIndex));
        Gedim::Output::Assert(cell0DIndex < Cell0DTotalNumber());
        _mesh.Cell2DSubdivision[_mesh.NumberCell2DSubdivision[cell2DIndex] + subDivisionIndex] = cell0DIndex;
    }
    inline unsigned int Cell2DNumberSubDivision(const unsigned int &cell2DIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        return _mesh.NumberCell2DSubdivision[cell2DIndex + 1] - _mesh.NumberCell2DSubdivision[cell2DIndex];
    }
    inline unsigned int Cell2DSubDivisionCell0D(const unsigned int &cell2DIndex, const unsigned int &subDivisionIndex) const
    {
        Gedim::Output::Assert(cell2DIndex < Cell2DTotalNumber());
        Gedim::Output::Assert(subDivisionIndex < Cell2DNumberSubDivision(cell2DIndex));
        return _mesh.Cell2DSubdivision[_mesh.NumberCell2DSubdivision[cell2DIndex] + subDivisionIndex];
    }

    void Cell3DsInitialize(const unsigned int &numberCell3Ds);
    unsigned int Cell3DAppend(const unsigned int &numberCell3Ds);
    void Cell3DRemove(const unsigned int &cell3DIndex);
    inline void Cell3DSetMarker(const unsigned int &cell3DIndex, const unsigned int &marker)
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        _mesh.Cell3DMarkers[cell3DIndex] = marker;
    }
    inline void Cell3DSetState(const unsigned int &cell3DIndex, const bool &state)
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        _mesh.ActiveCell3D[cell3DIndex] = state;
    }
    inline void Cell3DsInitializeVertices(const std::vector<unsigned int> &numberCell3DsVertices);
    inline void Cell3DsInitializeEdges(const std::vector<unsigned int> &numberCell3DsEdges);
    inline void Cell3DsInitializeFaces(const std::vector<unsigned int> &numberCell3DsFaces);
    void Cell3DInitializeVertices(const unsigned int &cell3DIndex, const unsigned int &numberCell3DVertices);
    void Cell3DInitializeEdges(const unsigned int &cell3DIndex, const unsigned int &numberCell3DEdges);
    void Cell3DInitializeFaces(const unsigned int &cell3DIndex, const unsigned int &numberCell3DFaces);
    inline void Cell3DInsertVertex(const unsigned int &cell3DIndex, const unsigned int &vertexIndex, const unsigned int &vertexCell0DIndex)
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell3DNumberVertices(cell3DIndex));
        Gedim::Output::Assert(vertexCell0DIndex < Cell0DTotalNumber());

        _mesh.Cell3DVertices[_mesh.NumberCell3DVertices[cell3DIndex] + vertexIndex] = vertexCell0DIndex;
    }
    inline void Cell3DInsertVertices(const unsigned int &cell3DIndex, const std::vector<unsigned int> &verticesCell0DIndices)
    {
        for (unsigned int v = 0; v < verticesCell0DIndices.size(); v++)
            Cell3DInsertVertex(cell3DIndex, v, verticesCell0DIndices[v]);
    }
    inline void Cell3DInsertEdge(const unsigned int &cell3DIndex, const unsigned int &edgeIndex, const unsigned int &edgeCell1DIndex)
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(edgeIndex < Cell3DNumberEdges(cell3DIndex));
        Gedim::Output::Assert(edgeCell1DIndex < Cell1DTotalNumber());

        _mesh.Cell3DEdges[_mesh.NumberCell3DEdges[cell3DIndex] + edgeIndex] = edgeCell1DIndex;
    }
    inline void Cell3DInsertEdges(const unsigned int &cell3DIndex, const std::vector<unsigned int> &edgesCell1DIndices)
    {
        for (unsigned int e = 0; e < edgesCell1DIndices.size(); e++)
            Cell3DInsertEdge(cell3DIndex, e, edgesCell1DIndices[e]);
    }
    inline void Cell3DInsertFace(const unsigned int &cell3DIndex, const unsigned int &faceIndex, const unsigned int &faceCell2DIndex)
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(faceIndex < Cell3DNumberFaces(cell3DIndex));
        Gedim::Output::Assert(faceCell2DIndex < Cell2DTotalNumber());

        _mesh.Cell3DFaces[_mesh.NumberCell3DFaces[cell3DIndex] + faceIndex] = faceCell2DIndex;
    }
    inline void Cell3DInsertFaces(const unsigned int &cell3DIndex, const std::vector<unsigned int> &facesCell2DIndices)
    {
        for (unsigned int f = 0; f < facesCell2DIndices.size(); f++)
            Cell3DInsertFace(cell3DIndex, f, facesCell2DIndices[f]);
    }
    void Cell3DAddVertices(const unsigned int &cell3DIndex, const std::vector<unsigned int> &verticesCell0DIndices);
    void Cell3DAddEdges(const unsigned int &cell3DIndex, const std::vector<unsigned int> &edgesCell1DIndices);
    void Cell3DAddFaces(const unsigned int &cell3DIndex, const std::vector<unsigned int> &facesCell2DIndices);

    unsigned int Cell3DFindVertex(const unsigned int &cell3DIndex, const unsigned int &cell0DIndex) const;
    unsigned int Cell3DFindEdge(const unsigned int &cell3DIndex, const unsigned int &cell1DIndex) const;
    unsigned int Cell3DFindFace(const unsigned int &cell3DIndex, const unsigned int &cell2DIndex) const;

    unsigned int Cell3DFindEdgeByExtremes(const unsigned int &cell3DIndex,
                                          const unsigned int &originCell0DIndex,
                                          const unsigned int &endCell0DIndex) const;

    inline unsigned int Cell3DTotalNumber() const
    {
        return _mesh.NumberCell3D;
    }
    inline unsigned int Cell3DNumberVertices(const unsigned int &cell3DIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.NumberCell3DVertices[cell3DIndex + 1] - _mesh.NumberCell3DVertices[cell3DIndex];
    }
    inline unsigned int Cell3DNumberEdges(const unsigned int &cell3DIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.NumberCell3DEdges[cell3DIndex + 1] - _mesh.NumberCell3DEdges[cell3DIndex];
    }
    inline unsigned int Cell3DNumberFaces(const unsigned int &cell3DIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.NumberCell3DFaces[cell3DIndex + 1] - _mesh.NumberCell3DFaces[cell3DIndex];
    }
    std::vector<unsigned int> Cell3DVertices(const unsigned int &cell3DIndex) const;
    inline unsigned int Cell3DVertex(const unsigned int &cell3DIndex, const unsigned int &vertexIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell3DNumberVertices(cell3DIndex));

        return _mesh.Cell3DVertices[_mesh.NumberCell3DVertices[cell3DIndex] + vertexIndex];
    }
    inline Eigen::Vector3d Cell3DVertexCoordinates(const unsigned int &cell3DIndex, const unsigned int &vertexIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(vertexIndex < Cell3DNumberVertices(cell3DIndex));
        return Cell0DCoordinates(Cell3DVertex(cell3DIndex, vertexIndex));
    }
    Eigen::MatrixXd Cell3DVerticesCoordinates(const unsigned int &cell3DIndex) const;
    std::vector<unsigned int> Cell3DEdges(const unsigned int &cell3DIndex) const;
    inline unsigned int Cell3DEdge(const unsigned int &cell3DIndex, const unsigned int &edgeIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(edgeIndex < Cell3DNumberEdges(cell3DIndex));

        return _mesh.Cell3DEdges[_mesh.NumberCell3DEdges[cell3DIndex] + edgeIndex];
    }
    std::vector<unsigned int> Cell3DFaces(const unsigned int &cell3DIndex) const;
    inline unsigned int Cell3DFace(const unsigned int &cell3DIndex, const unsigned int &faceIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(faceIndex < Cell3DNumberFaces(cell3DIndex));

        return _mesh.Cell3DFaces[_mesh.NumberCell3DFaces[cell3DIndex] + faceIndex];
    }
    std::vector<std::vector<std::vector<unsigned int>>> Cell3DsFacesVertices() const;
    std::vector<std::vector<unsigned int>> Cell3DsVertices() const;
    std::vector<std::vector<unsigned int>> Cell3DsEdges() const;
    std::vector<std::vector<unsigned int>> Cell3DsFaces() const;
    inline unsigned int Cell3DMarker(const unsigned int &cell3DIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.Cell3DMarkers[cell3DIndex];
    }
    inline std::vector<unsigned int> Cell3DsMarker() const
    {
        return _mesh.Cell3DMarkers;
    }
    inline bool Cell3DIsActive(const unsigned int &cell3DIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.ActiveCell3D[cell3DIndex];
    }
    inline std::vector<bool> Cell3DsState() const
    {
        return _mesh.ActiveCell3D;
    }

    inline bool Cell3DHasUpdatedCell3Ds(const unsigned int &cell3DIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.UpdatedCell3Ds.find(cell3DIndex) != _mesh.UpdatedCell3Ds.end();
    }
    inline unsigned int Cell3DNumberUpdatedCell3Ds(const unsigned int &cell3DIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        return _mesh.UpdatedCell3Ds.at(cell3DIndex).size();
    }
    inline bool Cell3DHasUpdatedCell3D(const unsigned int &cell3DIndex, const unsigned int &updatedCell3DIdex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(updatedCell3DIdex < Cell3DTotalNumber());
        return _mesh.UpdatedCell3Ds.at(cell3DIndex).find(updatedCell3DIdex) != _mesh.UpdatedCell3Ds.at(cell3DIndex).end();
    }
    void Cell3DInsertUpdatedCell3D(const unsigned int &cell3DIndex, const unsigned int &updatedCell3DIdex);
    bool Cell3DUpdatedCell3Ds(const unsigned int &cell3DIndex, std::list<unsigned int> &updatedCell3DIds) const;
    inline bool Cell3DHasOriginalCell3D(const unsigned int &updatedCell3DIndex) const
    {
        Gedim::Output::Assert(updatedCell3DIndex < Cell3DTotalNumber());
        return _mesh.Cell3DOriginalCell3Ds.at(updatedCell3DIndex) < _mesh.NumberCell3D;
    }
    inline unsigned int Cell3DOriginalCell3D(const unsigned int &updatedCell3DIndex) const
    {
        Gedim::Output::Assert(updatedCell3DIndex < Cell3DTotalNumber());
        return _mesh.Cell3DOriginalCell3Ds.at(updatedCell3DIndex);
    }

    void Cell3DInitializeDoubleProperties(const unsigned int &numberDoubleProperties);
    unsigned int Cell3DAddDoubleProperty(const std::string &propertyId);
    inline void Cell3DsInitializeDoublePropertyValues(const unsigned int &propertyIndex, const std::vector<unsigned int> &propertySizes);
    void Cell3DInitializeDoublePropertyValues(const unsigned int &cell3DIndex,
                                              const unsigned int &propertyIndex,
                                              const unsigned int &propertySize);
    inline void Cell3DInsertDoublePropertyValue(const unsigned int &cell3DIndex,
                                                const unsigned int &propertyIndex,
                                                const unsigned int &propertyValueIndex,
                                                const double &propertyValue)
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());

        _mesh.Cell3DDoublePropertyValues[propertyIndex][_mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex] + propertyValueIndex] =
            propertyValue;
    }
    inline unsigned int Cell3DNumberDoubleProperties() const
    {
        return _mesh.Cell3DDoublePropertyIds.size();
    }
    inline std::string Cell3DDoublePropertyId(const unsigned int &propertyIndex) const
    {
        return _mesh.Cell3DDoublePropertyIds[propertyIndex];
    }
    inline bool Cell3DDoublePropertyExists(const std::string &propertyId) const
    {
        return _mesh.Cell3DDoublePropertyIndices.find(propertyId) != _mesh.Cell3DDoublePropertyIndices.end();
    }
    inline unsigned int Cell3DDoublePropertyIndex(const std::string &propertyId) const
    {
        Gedim::Output::Assert(Cell3DDoublePropertyExists(propertyId));
        return _mesh.Cell3DDoublePropertyIndices.at(propertyId);
    }
    inline unsigned int Cell3DDoublePropertySize(const unsigned int &cell3DIndex, const unsigned int &propertyIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());
        return _mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex + 1] -
               _mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex];
    }
    inline double Cell3DDoublePropertyValue(const unsigned int &cell3DIndex,
                                            const unsigned int &propertyIndex,
                                            const unsigned int &propertyValueIndex) const
    {
        Gedim::Output::Assert(cell3DIndex < Cell3DTotalNumber());
        Gedim::Output::Assert(propertyIndex < Cell3DNumberDoubleProperties());
        Gedim::Output::Assert(propertyValueIndex < Cell3DDoublePropertySize(cell3DIndex, propertyIndex));

        return _mesh.Cell3DDoublePropertyValues[propertyIndex][_mesh.Cell3DDoublePropertySizes[propertyIndex][cell3DIndex] + propertyValueIndex];
    }

    void Compress();

    std::string ToString();
};
} // namespace Gedim

#endif
