// Copyright (C) 2014 Vicini Fabio
//
// This file is part of the dissertation of the author.
//
// This is a free program: you can redistribute it and/or modify
// it under the terms of the author.
//
// This program is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.
//
// Modified by Vicini Fabio 2014
//
// First added:  2014-10-05

#ifndef __GEDIM_IOUtilities_H
#define __GEDIM_IOUtilities_H

#include "Macro.hpp"
#include <unordered_set>

#if USE_MPI == 1
#include <mpi.h>
#endif

#include <cstdarg>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <unordered_map>

namespace Gedim
{
  class Output;
  class Profiler;
  class LogFile;

  /// @brief Class to print variables
  class Output
  {
    public:
      enum FileFilter
      {
        FilesAndDirectories = 0,
        Files = 1,
        Directories = 2
      };

      static std::string BlueColor;
      static std::string RedColor;
      static std::string GreenColor;
      static std::string YellowColor;
      static std::string EndColor;

      static int MaxElementToPrint; ///< Max elements to print in vectors
      static int StartingIndexToPrint; ///< Starting index to print in vectors

      /// Creating Folders if not exists
      static void CreateFolder(const std::string& nameFolder);
      /// Split a std::string in more strings dived by character
      static std::vector<std::string> StringSplit(const std::string& stringToSplit,
                                                  const char character);
      /// Check if a file exists
      static bool FileExists(const std::string& nameFile);
      /// Get the path from a file path
      static void GetFilePath(const std::string& nameFilePath,
                              std::string& filePath,
                              std::string& fileName,
                              std::string& fileExtension);
      /// Returns a list of files/directories in mainDirectory
      static void GetFileList(const std::string& mainDirectory,
                              std::vector<std::string>& out,
                              const FileFilter& filter = Output::FilesAndDirectories,
                              const bool& hiddenElements = false);
      /// Returns a list of paths containing the file/directory found in mainDirectory
      static void FindPaths(const std::string& mainDirectory,
                            const std::string& objectNameToFind,
                            std::vector<std::string>& paths,
                            const FileFilter& filter = Output::FilesAndDirectories,
                            const bool& hiddenElements = false);

      static void GetBinaryFileSize(const std::string& nameFile,
                                    unsigned int& fileSize,
                                    const unsigned int& sizeOfSingleDataToWrite,
                                    const unsigned int& startingPosition = 0);
      static void ReadBinaryFile(const std::string& nameFile,
                                 void* dataToRead,
                                 const unsigned int& sizeOfSingleDataToWrite,
                                 const unsigned int& dataSizeToRead,
                                 const unsigned int& startingPosition = 0);
      static bool ReadBinaryFile(const std::string& nameFile,
                                 std::vector<double>& dataToRead,
                                 const unsigned int& dataSizeToRead = 0,
                                 const unsigned int& startingPosition = 0);
      static void WriteBinaryFile(const std::string& nameFile,
                                  const void* dataToWrite,
                                  const unsigned int& sizeOfSingleDataToWrite,
                                  const unsigned int& dataSizeToWrite,
                                  const bool& append = false);
      static bool WriteBinaryFile(const std::string& nameFile,
                                  const std::vector<double>& dataToWrite,
                                  const unsigned int& dataSizeToWrite = 0,
                                  const unsigned int& dataStartingPositionToWrite = 0, const bool& append = false);

      /// Print a line of symbol
      static void PrintLine(char character = ' ', bool onlyMaster = true);

      /// Print a line of stars * in \p out.
      static std::ostream& PrintStars(std::ostream& out);
      /// Print a line of lines - in \p out.
      static std::ostream& PrintLines(std::ostream& out);

      /// Used to print to file, or on screen the status of the program.
      static void PrintStatusProgram(const std::string& programStep, ...);
      /// Used to print only on screen a generic message in the same line
      static void PrintGenericMessageOnLine(const std::string& message, ...);
      /// Used to print to file, or on screen a generic message
      static void PrintGenericMessage(const std::string& message, const bool& onlyMaster, ...);
      /// Used to print to file, or on screen an error message
      static void PrintErrorMessage(const std::string& message, const bool& onlyMaster, ...);
      /// Used to print to file, or on screen a warning message
      static void PrintWarningMessage(const std::string& message, const bool& onlyMaster, ...);
      /// Used to print to file, or on screen a debug message
      static void PrintDebugMessage(const std::string& message, const bool& onlyMaster, ...);
      /// Used to print to file, or on screen a success message
      static void PrintSuccessMessage(const std::string& message, const bool& onlyMaster, ...);

      /// Assert for all code, generate exception if something goes wrong
      inline static void Assert(const bool& logicResult)
      {
        if (!logicResult)
          abort();
      }

      /// Assert for all code, generate exception if something goes wrong with a message
      static void Assert(const bool& logicResult, const std::string& message, ...);
  };

  class LogFile
  {
    private:
      static std::string GetDateTime();
      static void PrintMessage(const std::string& type, const std::string& message, va_list args);

    public:
      static int LogMaxFileSize; ///< Max file log size (MB)
      static int LogMaxNumFiles; ///< Max number of backup log files to maintain
      static std::string LogFolder; ///< Folder name where log files of the program are
      static std::string LogNameFile; ///< File name of Log

      static void PrintLine(const char& symbol = ' ');
      static void PrintWarningMessage(const std::string& message, va_list args);
      static void PrintErrorMessage(const std::string& message, va_list args);
      static void PrintInfoMessage(const std::string& message, va_list args);
      static void PrintDebugMessage(const std::string& message, va_list args);

      /// Check if file reach the maximum size and create new file
      static void CheckFileSize(const std::string& nameFile);
      /// Get the file size in MB
      static double GetFileSize(const std::string& nameFile);
  };

  class Profiler
  {
    private:
      static std::map<std::string, double> times; ///< list of times
      static std::map<std::string, double> localTimes; ///< list of times
      static unsigned long totalAvailStartingMemory; ///< Avail Memory on program start
      static int ParseLine(char* line); ///< used in the CheckMemory function

    public:
      static bool ActivateProfiler; ///< Compute time for program profiling
      static std::string TimeFile; ///< File where time profiling data are stored
      static std::string MemoryFile; ///< File where memory profiling data are stored

      static double GetTime(const std::string& nameTime, const bool& localTime = false);
      /// Get time in ms
      static double GetTime();
      static double ComputeTime(const double& startTime,
                                const double& stopTime)
      {
#if USE_MPI == 1
        return stopTime - startTime;
#elif USE_MPI == 0
        return (stopTime - startTime) / (double)CLOCKS_PER_SEC;
#endif // USE_MPI
      }

      /// \brief Get Total Physical Memory from file /proc/meminfo (KB)
      /// \param totalMemory will have 3 values: TotMemory, UsedMemory, AvailMemory
      static void GetTotalMemory(std::vector<unsigned long>& totalMemory);

      /// \brief Get Current Virtual/Physical Memory used by process (KB)
      /// \param memoryUsed will have 2 values: virtual, physical
      static void GetProcessMemory(std::vector<unsigned long>& memoryUsed);

      /// \brief Compute the virtual memory used by the process
      /// \details Get from http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
      static void CheckMemory(const std::string& nameMemory, const bool& checkProcessesMemory = false);

      /// Sleep for n milliseconds
      static void Sleep(const unsigned int& milliseconds);

      /// Start time with nameTime if not exist
      static void StartTime(const std::string& nameTime);
      /// Insert time with nameTime if not exist
      static void InsertTime(const std::string& nameTime,
                             const double& globalTime,
                             const double& localTime);
      /// Split time with nameTime if exists and return the global time and local time
      static void SplitTime(const std::string& nameTime,
                            double& globalTime,
                            double& localTime);
      /// Stop time with nameTime if exists and print to file the result
      static double StopTime(const std::string& nameTime,
                             const bool& printTime = true);

      static void WriteTime(const std::string& nameTime,
                            const double& globalTime,
                            const double& localTime);
  };

  template <typename matrixType>
  std::string MatrixToString(const matrixType& mat,
                             const std::string& matrixTypeStr,
                             const std::string& matrixName)
  {
    std::ostringstream str;
    str.precision(16);
    str<< matrixName<< " = "<< matrixTypeStr<< "("<< mat.rows()<< ", "<< mat.cols()<< ");"<< std::endl;
    for (unsigned int c = 0; c < mat.cols(); c++)
    {
      str<< matrixName<< ".col("<< c<< ")<< ";
      for (unsigned int r = 0; r < mat.rows(); r++)
        str<< std::scientific<< (r == 0 ? "" : ",")<< mat(r, c);
      str<< ";"<< std::endl;
    }

    return str.str();
  }

  template <typename matrixType>
  std::string MatrixCollectionToString(const std::vector<matrixType>& matCollection,
                                       const std::string& matrixTypeStr,
                                       const std::string& matrixName)
  {
    std::ostringstream str;
    str.precision(16);

    str<< matrixName<< " = std::vector<"<< matrixTypeStr<< ">("<< matCollection.size()<< ");"<< std::endl;
    for (unsigned int v = 0; v < matCollection.size(); v++)
    {
      str<< MatrixToString<matrixType>(matCollection.at(v),
                                       matrixTypeStr,
                                       matrixName + "[" + std::to_string(v) + "]");
    }

    return str.str();
  }

  template <typename matrixType>
  std::string MatrixCollectionToString(const std::vector<std::vector<matrixType>>& matCollection,
                                       const std::string& matrixTypeStr,
                                       const std::string& matrixName)
  {
    std::ostringstream str;
    str.precision(16);

    str<< matrixName<< " = std::vector<std::vector<"<< matrixTypeStr<< ">>("<< matCollection.size()<< ");"<< std::endl;
    for (unsigned int v = 0; v < matCollection.size(); v++)
    {
      str<< MatrixCollectionToString<matrixType>(matCollection.at(v),
                                                 matrixTypeStr,
                                                 matrixName + "[" + std::to_string(v) + "]");
    }

    return str.str();
  }

  template <typename matrixType>
  std::string MatrixCollectionToString(const std::vector<std::vector<std::vector<matrixType>>>& matCollection,
                                       const std::string& matrixTypeStr,
                                       const std::string& matrixName)
  {
    std::ostringstream str;
    str.precision(16);

    str<< matrixName<< " = std::vector<std::vector<std::vector<"<< matrixTypeStr<< ">>>("<< matCollection.size()<< ");"<< std::endl;
    for (unsigned int v = 0; v < matCollection.size(); v++)
    {
      str<< MatrixCollectionToString<matrixType>(matCollection.at(v),
                                                 matrixTypeStr,
                                                 matrixName + "[" + std::to_string(v) + "]");
    }

    return str.str();
  }

  /// General print of a vector
  template <typename T>
  std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeVector = (maxElementToPrint == 0 || maxElementToPrint > (int)vec.size()) ? vec.size() : maxElementToPrint;
    unsigned int startIndexVector = (startingIndex >= (int)vec.size()) ? 0 : startingIndex;

    out<< "[";
    for (unsigned int i = startIndexVector; i < startIndexVector + sizeVector && i < vec.size(); i++)
      out<< (i != startIndexVector ? "," : "")<< vec.at(i);
    out<< "]";

    return out;
  }
  /// General print of a vector
  template <typename T>
  std::ostream& operator<<(std::ostream& out, const std::vector<T*>& vec)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeVector = (maxElementToPrint == 0 || maxElementToPrint > (int)vec.size()) ? vec.size() : maxElementToPrint;
    unsigned int startIndexVector = (startingIndex >= (int)vec.size()) ? 0 : startingIndex;

    out<< "[";
    for (unsigned int i = startIndexVector; i < startIndexVector + sizeVector && i < vec.size(); i++)
      out<< (i != startIndexVector ? "," : "")<< *vec.at(i);
    out<< "]";

    return out;
  }
  /// General print of an array
  template <typename T, size_t s>
  std::ostream& operator<<(std::ostream& out, const std::array<T, s>& vec)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeVector = (maxElementToPrint == 0 || maxElementToPrint > (int)vec.size()) ? vec.size() : maxElementToPrint;
    unsigned int startIndexVector = (startingIndex >= (int)vec.size()) ? 0 : startingIndex;

    out<< "[";
    for (unsigned int i = startIndexVector; i < startIndexVector + sizeVector && i < vec.size(); i++)
      out<< (i != startIndexVector ? "," : "")<< vec.at(i);
    out<< "]";

    return out;
  }
  /// General print of an array
  template <typename T, size_t s>
  std::ostream& operator<<(std::ostream& out, const std::array<T*, s>& vec)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeVector = (maxElementToPrint == 0 || maxElementToPrint > (int)vec.size()) ? vec.size() : maxElementToPrint;
    unsigned int startIndexVector = (startingIndex >= (int)vec.size()) ? 0 : startingIndex;

    out<< "[";
    for (unsigned int i = startIndexVector; i < startIndexVector + sizeVector && i < vec.size(); i++)
      out<< (i != startIndexVector ? "," : "")<< *vec.at(i);
    out<< "]";

    return out;
  }
  /// General print of a list
  template <typename T>
  std::ostream& operator<<(std::ostream& out,
                           const std::list<T>& listToPrint)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)listToPrint.size()) ? listToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)listToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out<< "{";
    for(typename std::list<T>::const_iterator iterator = listToPrint.begin(); iterator != listToPrint.end(); iterator++)
    {
      if (counter >= startIndexList)
      {
        out<< (counter != startIndexList ? "," : "")<< *iterator;
        elementPrinted++;
      }

      if (elementPrinted >= sizeList)
        break;

      counter++;
    }
    out<< "}";

    return out;
  }
  /// General print of a list
  template <typename T>
  std::ostream& operator<<(std::ostream& out,
                           const std::list<T*>& listToPrint)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)listToPrint.size()) ? listToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)listToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out<< "{";
    for(typename std::list<T*>::const_iterator iterator = listToPrint.begin(); iterator != listToPrint.end(); iterator++)
    {
      if (counter >= startIndexList)
      {
        out<< (counter != startIndexList ? "," : "")<< **iterator;
        elementPrinted++;
      }

      if (elementPrinted >= sizeList)
        break;

      counter++;
    }
    out<< "}";

    return out;
  }
  /// General print of a unordered_set
  template <typename T>
  std::ostream& operator<<(std::ostream& out,
                           const std::unordered_set<T>& setToPrint)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)setToPrint.size()) ? setToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)setToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out<< "{";
    for(typename std::unordered_set<T>::const_iterator iterator = setToPrint.begin(); iterator != setToPrint.end(); iterator++)
    {
      if (counter >= startIndexList)
      {
        out<< (counter != startIndexList ? "," : "")<< *iterator;
        elementPrinted++;
      }

      if (elementPrinted >= sizeList)
        break;

      counter++;
    }
    out<< "}";

    return out;
  }
  /// General print of a unordered_set
  template <typename T>
  std::ostream& operator<<(std::ostream& out,
                           const std::unordered_set<T*>& setToPrint)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)setToPrint.size()) ? setToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)setToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out<< "{";
    for(typename std::unordered_set<T*>::const_iterator iterator = setToPrint.begin(); iterator != setToPrint.end(); iterator++)
    {
      if (counter >= startIndexList)
      {
        out<< (counter != startIndexList ? "," : "")<< **iterator;
        elementPrinted++;
      }

      if (elementPrinted >= sizeList)
        break;

      counter++;
    }
    out<< "}";

    return out;
  }
  /// General print of a set
  template <typename T>
  std::ostream& operator<<(std::ostream& out,
                           const std::set<T>& setToPrint)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)setToPrint.size()) ? setToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)setToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out<< "{";
    for(typename std::set<T>::const_iterator iterator = setToPrint.begin(); iterator != setToPrint.end(); iterator++)
    {
      if (counter >= startIndexList)
      {
        out<< (counter != startIndexList ? "," : "")<< *iterator;
        elementPrinted++;
      }

      if (elementPrinted >= sizeList)
        break;

      counter++;
    }
    out<< "}";

    return out;
  }
  /// General print of a set
  template <typename T>
  std::ostream& operator<<(std::ostream& out,
                           const std::set<T*>& setToPrint)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)setToPrint.size()) ? setToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)setToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out<< "{";
    for(typename std::set<T*>::const_iterator iterator = setToPrint.begin(); iterator != setToPrint.end(); iterator++)
    {
      if (counter >= startIndexList)
      {
        out<< (counter != startIndexList ? "," : "")<< **iterator;
        elementPrinted++;
      }

      if (elementPrinted >= sizeList)
        break;

      counter++;
    }
    out<< "}";

    return out;
  }
  /// General print of a map
  template<typename map_key, typename map_val>
  std::ostream& operator<<(std::ostream& out,
                           const std::map<map_key, map_val>& mapToPrint)
  {
    unsigned int counter = 0;

    out<< "{";
    for (typename std::map<map_key, map_val>::const_iterator it = mapToPrint.begin();
         it != mapToPrint.end();
         ++it)
    {
      out<< (counter != 0 ? "," : "") << "{"<< it->first << ", " << it->second<< "}";
      counter++;
    }
    out<< "}";

    return out;
  }
  /// General print of a map
  template<typename map_key, typename map_val>
  std::ostream& operator<<(std::ostream& out,
                           const std::map<map_key, map_val*>& mapToPrint)
  {
    unsigned int counter = 0;

    out<< "{";
    for (typename std::map<map_key, map_val*>::const_iterator it = mapToPrint.begin();
         it != mapToPrint.end();
         ++it)
    {
      out<< (counter != 0 ? "," : "") << "{"<< it->first << ", " << (it->second == nullptr ? "NULL" : (*it->second))<< "}";
      counter++;
    }
    out<< "}";

    return out;
  }
  /// General print of a unordered map
  template<typename map_key, typename map_val>
  std::ostream& operator<<(std::ostream& out,
                           const std::unordered_map<map_key, map_val>& mapToPrint)
  {
    unsigned int counter = 0;

    out<< "{";
    for (typename std::unordered_map<map_key, map_val>::const_iterator it = mapToPrint.begin();
         it != mapToPrint.end();
         ++it)
    {
      out<< (counter != 0 ? "," : "") << "{"<< it->first << ", " << it->second<< "}";
      counter++;
    }
    out<< "}";

    return out;
  }

}

#endif
