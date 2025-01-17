#ifndef __GEDIM_FILEREADER_H
#define __GEDIM_FILEREADER_H

#include "IFileTextReader.hpp"

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
