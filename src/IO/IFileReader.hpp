#ifndef __GEDIM_IFILEREADER_H
#define __GEDIM_IFILEREADER_H

#include "IOUtilities.hpp"

using namespace std;

namespace Gedim
{
  /// \brief Interface for File Reader
  /// \copyright See top level LICENSE file for details.
  class IFileReader
  {
    public:
      virtual ~IFileReader() {}

      /// \return the file path
      virtual string Path() = 0;
      /// \brief Open the file
      /// \return if the open is successfull
      virtual bool Open() = 0;
      /// \brief Jump line
      virtual void NextLine() = 0;
      /// \brief Get a single line in file
      virtual void GetLine(string& line) = 0;
      /// \brief Get all lines
      virtual void GetAllLines(vector<string>& lines) = 0;
      /// \brief Close the file
      virtual void Close() = 0;
  };
}

#endif // __GEDIM_IFILEREADER_H
