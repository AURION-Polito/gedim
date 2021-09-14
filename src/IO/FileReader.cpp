#include "FileReader.hpp"

namespace Gedim
{
  // ***************************************************************************
  FileReader::FileReader(const string& filePath)
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
    string line;
    getline(_file, line);
  }
  // ***************************************************************************
  void FileReader::GetAllLines(vector<string>& lines)
  {
    list<string> listLines;
    string line;
    while (getline(_file, line))
      listLines.push_back(line);

    lines.reserve(listLines.size());
    for (const string& line : listLines)
      lines.push_back(line);
  }
  // ***************************************************************************
}
