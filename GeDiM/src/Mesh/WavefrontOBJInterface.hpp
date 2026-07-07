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

#ifndef __WavefrontOBJInterface_H
#define __WavefrontOBJInterface_H

#include "IMeshDAO.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
class WavefrontOBJInterface final
{
  public:
    struct OBJMesh final
    {
        unsigned int NumCell0Ds;
        unsigned int NumCell1Ds;
        unsigned int NumCell2Ds;

        Eigen::MatrixXd Cell0Ds;
        std::vector<Eigen::VectorXi> Cell2Ds;
    };

    WavefrontOBJInterface();
    ~WavefrontOBJInterface();

    OBJMesh StringsToOBJMesh(const std::vector<std::string> &fileLines) const;
    std::vector<std::string> OBJMeshToStrings(const OBJMesh &mesh) const;
    OBJMesh MeshDAOToOBJMesh(const IMeshDAO &originalMesh) const;
    void OBJMeshToMeshDAO(const OBJMesh &originalMesh, const MeshUtilities &meshUtilities, IMeshDAO &convertedMesh) const;

    std::vector<std::string> ImportMeshFromFile(const std::string &OBJFilePath, const MeshUtilities &meshUtilities, IMeshDAO &mesh) const;

    void ExportMeshToFile(const IMeshDAO &mesh, const std::string &OBJFilePath) const;
};
} // namespace Gedim

#endif // __ObjectFileFormatInterface_H
