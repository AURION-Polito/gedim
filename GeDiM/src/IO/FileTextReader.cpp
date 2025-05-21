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

#include "FileTextReader.hpp"

namespace Gedim
{
// ***************************************************************************
FileReader::FileReader(const std::string &filePath)
{
    _filePath = filePath;
}
// ***************************************************************************
bool FileReader::Open()
{
    _file.open(_filePath);
    return !_file.fail();
}
// ***************************************************************************
void FileReader::NextLine()
{
    std::string line;
    getline(_file, line);
}
// ***************************************************************************
void FileReader::GetAllLines(std::vector<std::string> &lines)
{
    std::list<std::string> listLines;
    std::string line;
    while (getline(_file, line))
        listLines.push_back(line);

    lines.reserve(listLines.size());
    for (const std::string &line : listLines)
        lines.push_back(line);
}
// ***************************************************************************
} // namespace Gedim
