#ifndef __FILEREADER_H
#define __FILEREADER_H

#include "IFileReader.hpp"

using namespace std;

namespace Gedim
{
  /// \brief C++ File Reader
  /// \copyright See top level LICENSE file for details.
  class FileReader : public IFileReader
  {
    private:
      ifstream _file;
      string _filePath;

    public:
      FileReader(const string& filePath);
      virtual ~FileReader() {}

      string Path() { return _filePath; }
      bool Open();
      void NextLine();
      void GetLine(string& line) { getline(_file, line); }
      void GetAllLines(vector<string>& lines);
      void Close() { _file.close(); }
  };
}

#endif // __FILEREADER_H
