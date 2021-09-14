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

#ifndef __OUTPUT_H
#define __OUTPUT_H

#include "MacroDefinitions.hpp"

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

using namespace std;

namespace GeDiM
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

      /// \brief Exit codes for methods throughout GeDiM
      enum ExitCodes
      {
        Success = 1, ///< Success flag.
        Abort = -1, ///< Abort flag.
        MpiError = -2, ///< MPI error flag.
        PartitionError = -3, ///< Partitioning error flag.
        FileError = -4, ///< File error flag.
        GenericError = -5, ///< Generic error flag.
        UnimplementedMethod = -6 ///< Flag for an unimplemented feature.
      };

      static string BlueColor;
      static string RedColor;
      static string GreenColor;
      static string YellowColor;
      static string EndColor;

      static int MaxElementToPrint; ///< Max elements to print in vectors
      static int StartingIndexToPrint; ///< Starting index to print in vectors

      /// Creating Folders if not exists
      static void CreateFolder(const string& nameFolder);
      /// Split a string in more strings dived by character
      static vector<string> StringSplit(const string& stringToSplit, const char character);
      /// Check if a file exists
      static bool FileExists(const string& nameFile);
      /// Get the path from a file path
      static void GetFilePath(const string& nameFilePath, string& filePath, string& fileName, string& fileExtension);
      /// Returns a list of files/directories in mainDirectory
      static void GetFileList(const string& mainDirectory, vector<string>& out, const FileFilter& filter = Output::FilesAndDirectories, const bool& hiddenElements = false);
      /// Returns a list of paths containing the file/directory found in mainDirectory
      static void FindPaths(const string& mainDirectory, const string& objectNameToFind, vector<string>& paths, const FileFilter& filter = Output::FilesAndDirectories, const bool& hiddenElements = false);

      static Output::ExitCodes GetBinaryFileSize(const string& nameFile, unsigned int& fileSize, const unsigned int& sizeOfSingleDataToWrite, const unsigned int& startingPosition = 0);
      static Output::ExitCodes ReadBinaryFile(const string& nameFile, void* dataToRead, const unsigned int& sizeOfSingleDataToWrite, const unsigned int& dataSizeToRead, const unsigned int& startingPosition = 0);
      static bool ReadBinaryFile(const string& nameFile, vector<double>& dataToRead, const unsigned int& dataSizeToRead = 0, const unsigned int& startingPosition = 0);
      static Output::ExitCodes WriteBinaryFile(const string& nameFile, const void* dataToWrite, const unsigned int& sizeOfSingleDataToWrite, const unsigned int& dataSizeToWrite, const bool& append = false);
      static bool WriteBinaryFile(const string& nameFile, const vector<double>& dataToWrite, const unsigned int& dataSizeToWrite = 0, const unsigned int& dataStartingPositionToWrite = 0, const bool& append = false);

      /// Print a line of symbol
      static void PrintLine(char character = ' ', bool onlyMaster = true);

      /// Print a line of stars * in \p out.
      static ostream& PrintStars(ostream& out);
      /// Print a line of lines - in \p out.
      static ostream& PrintLines(ostream& out);

      /// Used to print to file, or on screen the status of the program.
      static void PrintStatusProgram(const string& programStep, ...);
      /// Used to print only on screen a generic message in the same line
      static void PrintGenericMessageOnLine(const string& message, ...);
      /// Used to print to file, or on screen a generic message
      static void PrintGenericMessage(const string& message, const bool& onlyMaster, ...);
      /// Used to print to file, or on screen an error message
      static void PrintErrorMessage(const string& message, const bool& onlyMaster, ...);
      /// Used to print to file, or on screen a warning message
      static void PrintWarningMessage(const string& message, const bool& onlyMaster, ...);
      /// Used to print to file, or on screen a debug message
      static void PrintDebugMessage(const string& message, const bool& onlyMaster, ...);
      /// Used to print to file, or on screen a success message
      static void PrintSuccessMessage(const string& message, const bool& onlyMaster, ...);

      /// ExitCode result to string
      static const string ExitCodeToString(const ExitCodes& result, const bool& colored = false);

      /// Assert for all code, generate exception if something goes wrong
      inline static void Assert(const bool& logicResult)
      {
        if (!logicResult)
          abort();
      }
      /// Assert for all code, generate exception if something goes wrong
      inline static void Assert(const ExitCodes& logicResult)
      {
        if (logicResult != Output::Success)
          abort();
      }

      /// Assert for all code, generate exception if something goes wrong with a message
      static void Assert(const bool& logicResult, const string& message, ...);
      /// Assert for all code, generate exception if something goes wrong with a message
      static void Assert(const ExitCodes& logicResult, const string& message, ...);

      /// Assert for unit test, generate exception if something goes wrong with a message
      static void AssertTest(const bool& logicResult, const string& message, ...);
      /// Assert for unit test, generate exception if something goes wrong with a message
      static void AssertTest(const ExitCodes& logicResult, const string& message, ...);
  };

  class LogFile
  {
    private:
      static string GetDateTime();
      static void PrintMessage(const string& type, const string& message, va_list args);

    public:
      static int LogMaxFileSize; ///< Max file log size (MB)
      static int LogMaxNumFiles; ///< Max number of backup log files to maintain
      static string LogFolder; ///< Folder name where log files of the program are
      static string LogNameFile; ///< File name of Log

      static void PrintLine(const char& symbol = ' ');
      static void PrintWarningMessage(const string& message, va_list args);
      static void PrintErrorMessage(const string& message, va_list args);
      static void PrintInfoMessage(const string& message, va_list args);
      static void PrintDebugMessage(const string& message, va_list args);

      /// Check if file reach the maximum size and create new file
      static void CheckFileSize(const string& nameFile);
      /// Get the file size in MB
      static double GetFileSize(const string& nameFile);
  };

  class Profiler
  {
    private:
      static map<string, double> times; ///< list of times
      static map<string, double> localTimes; ///< list of times
      static unsigned long totalAvailStartingMemory; ///< Avail Memory on program start
      static int ParseLine(char* line); ///< used in the CheckMemory function

    public:
      static bool ActivateProfiler; ///< Compute time for program profiling
      static string TimeFile; ///< File where time profiling data are stored
      static string MemoryFile; ///< File where memory profiling data are stored

      static double GetTime(const string& nameTime, const bool& localTime = false);

      /// \brief Get Total Physical Memory from file /proc/meminfo (KB)
      /// \param totalMemory will have 3 values: TotMemory, UsedMemory, AvailMemory
      static void GetTotalMemory(vector<unsigned long>& totalMemory);

      /// \brief Get Current Virtual/Physical Memory used by process (KB)
      /// \param memoryUsed will have 2 values: virtual, physical
      static void GetProcessMemory(vector<unsigned long>& memoryUsed);

      /// \brief Compute the virtual memory used by the process
      /// \details Get from http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
      static void CheckMemory(const string& nameMemory, const bool& checkProcessesMemory = false);

      /// Sleep for n milliseconds
      static void Sleep(const unsigned int& milliseconds);

      /// Start time with nameTime if not exist
      static void StartTime(const string& nameTime);
      /// Split time with nameTime if exists and return the global time and local time
      static void SplitTime(const string& nameTime, double& globalTime, double& localTime);
      /// Stop time with nameTime if exists and print to file the result
      static void StopTime(const string& nameTime, const bool& printTime = true);
  };

  /// General print of a vector
  template <typename T>
  ostream& operator<<(ostream& out, const vector<T>& vec)
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
  ostream& operator<<(ostream& out, const vector<T*>& vec)
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
  ostream& operator<<(ostream& out, const list<T>& listToPrint)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)listToPrint.size()) ? listToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)listToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out<< "{";
    for(typename list<T>::const_iterator iterator = listToPrint.begin(); iterator != listToPrint.end(); iterator++)
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
  ostream& operator<<(ostream& out, const list<T*>& listToPrint)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)listToPrint.size()) ? listToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)listToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out<< "{";
    for(typename list<T*>::const_iterator iterator = listToPrint.begin(); iterator != listToPrint.end(); iterator++)
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
  /// General print of a list
  template <typename T>
  ostream& operator<<(ostream& out, const set<T>& setToPrint)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)setToPrint.size()) ? setToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)setToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out<< "{";
    for(typename set<T>::const_iterator iterator = setToPrint.begin(); iterator != setToPrint.end(); iterator++)
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
  ostream& operator<<(ostream& out, const set<T*>& setToPrint)
  {
    const int& maxElementToPrint = Output::MaxElementToPrint;
    const int& startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)setToPrint.size()) ? setToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)setToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out<< "{";
    for(typename set<T*>::const_iterator iterator = setToPrint.begin(); iterator != setToPrint.end(); iterator++)
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

}

#endif
