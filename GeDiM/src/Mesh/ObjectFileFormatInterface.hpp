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

#ifndef __ObjectFileFormatInterface_H
#define __ObjectFileFormatInterface_H

#include "IMeshDAO.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
class ObjectFileFormatInterface final
{
  public:
    struct OFFMesh final
    {
        unsigned int NumCell0Ds;
        unsigned int NumCell1Ds;
        unsigned int NumCell2Ds;

        Eigen::MatrixXd Cell0Ds;
        std::vector<Eigen::VectorXi> Cell2Ds;
    };

    ObjectFileFormatInterface();
    ~ObjectFileFormatInterface();

    OFFMesh StringsToOFFMesh(const std::vector<std::string> &fileLines) const;
    std::vector<std::string> OFFMeshToStrings(const OFFMesh &mesh) const;
    OFFMesh MeshDAOToOFFMesh(const IMeshDAO &originalMesh) const;
    void OFFMeshToMeshDAO(const OFFMesh &originalMesh, const MeshUtilities &meshUtilities, IMeshDAO &convertedMesh) const;

    void ImportMeshFromFile(const std::string &offFilePath, const MeshUtilities &meshUtilities, IMeshDAO &mesh) const;

    void ExportMeshToFile(const IMeshDAO &mesh, const std::string &offFilePath) const;
};
} // namespace Gedim

#endif // __ObjectFileFormatInterface_H
