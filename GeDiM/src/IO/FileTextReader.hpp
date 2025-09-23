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

#ifndef __GEDIM_FILEREADER_H
#define __GEDIM_FILEREADER_H

#include "IFileTextReader.hpp"

#include <fstream>

namespace Gedim
{
/// \brief C++ File Reader
/// \copyright See top level LICENSE file for details.
class FileReader : public IFileReader
{
  private:
    std::ifstream _file;
    std::string _filePath;

  public:
    FileReader(const std::string &filePath);
    virtual ~FileReader()
    {
    }

    std::string Path()
    {
        return _filePath;
    }
    bool Open();
    void NextLine();
    void GetLine(std::string &line)
    {
        getline(_file, line);
    }
    void GetAllLines(std::vector<std::string> &lines);
    void Close()
    {
        _file.close();
    }
};
} // namespace Gedim

#endif // __GEDIM_FILEREADER_H
