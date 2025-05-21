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

#ifndef __GEDIM_IFILEREADER_H
#define __GEDIM_IFILEREADER_H

#include "IOUtilities.hpp"

namespace Gedim
{
/// \brief Interface for File Reader
/// \copyright See top level LICENSE file for details.
class IFileReader
{
  public:
    virtual ~IFileReader()
    {
    }

    /// \return the file path
    virtual std::string Path() = 0;
    /// \brief Open the file
    /// \return if the open is successfull
    virtual bool Open() = 0;
    /// \brief Jump line
    virtual void NextLine() = 0;
    /// \brief Get a single line in file
    virtual void GetLine(std::string &line) = 0;
    /// \brief Get all lines
    virtual void GetAllLines(std::vector<std::string> &lines) = 0;
    /// \brief Close the file
    virtual void Close() = 0;
};
} // namespace Gedim

#endif // __GEDIM_IFILEREADER_H
