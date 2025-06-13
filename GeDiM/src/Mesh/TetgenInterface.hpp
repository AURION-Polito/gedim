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

#ifndef __TetgenInterface_H
#define __TetgenInterface_H

#include "Gedim_Macro.hpp"
#include "IMeshDAO.hpp"

#if ENABLE_TETGEN == 1
#include "tetgen.h"
#endif

#include "Eigen/Eigen"

namespace Gedim
{
/// \brief The Tetgen Interface
/// \see http://wias-berlin.de/software/tetgen/files/tetcall.cxx
class TetgenInterface final
{
  private:
#if ENABLE_TETGEN == 1
    void CreateTetgenInput(const Eigen::MatrixXd &polyhedronVertices,
                           const Eigen::MatrixXi &polyhedronEdges,
                           const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                           tetgenio &tetgenInput,
                           const Eigen::MatrixXd &constrainedPoints = Eigen::MatrixXd(),
                           const std::vector<Eigen::VectorXi> &constrainedFaces = std::vector<Eigen::VectorXi>()) const;
    void CreateTetgenInput(const Eigen::MatrixXd &points, const std::vector<std::vector<unsigned int>> &facets, tetgenio &tetgenInput) const;
    void CreateDelaunayInput(const Eigen::MatrixXd &points, const std::vector<unsigned int> &points_marker, tetgenio &tetgenInput) const;
    void CreateTetgenOutput(const double &maxTetrahedronArea,
                            tetgenio &tetgenInput,
                            tetgenio &tetgenOutput,
                            const std::string &tetgenOptions = "Qpqfezna") const;
    void CreateTetgenOutput(tetgenio &tetgenInput, tetgenio &tetgenOutput, const std::string &tetgenOptions) const;
    void ConvertTetgenOutputToMeshDAO(const tetgenio &tetgenOutput, IMeshDAO &mesh) const;

    void DeleteTetgenStructure(tetgenio &tetgenInput, tetgenio &tetgenOutput) const;
    void ExportTetgenOutput(const std::string &nameFolder, const std::string &nameFile, tetgenio &tetgenOutput) const;
#endif

  public:
    TetgenInterface();
    ~TetgenInterface();

    void CreateDelaunay(const Eigen::MatrixXd &points, const std::vector<unsigned int> &points_marker, IMeshDAO &mesh) const;

    void CreateMesh(const Eigen::MatrixXd &polyhedronVertices,
                    const Eigen::MatrixXi &polyhedronEdges,
                    const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                    const double &maxTetrahedronVolume,
                    IMeshDAO &mesh,
                    const std::string &tetgenOptions = "Qpqfezna") const;

    void CreateMesh(const Eigen::MatrixXd &points,
                    const std::vector<std::vector<unsigned int>> &facets,
                    const double &maxTetrahedronVolume,
                    IMeshDAO &mesh,
                    const std::string &tetgenOptions = "Qpqfezna") const;
};

} // namespace Gedim

#endif // __TetgenInterface_H
