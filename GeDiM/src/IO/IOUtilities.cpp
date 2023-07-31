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

#include "IOUtilities.hpp"

#include "MpiParallelEnvironment.hpp"
#include "CommonUtilities.hpp"

#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <limits>
#include <dirent.h>

#ifdef _WIN32
#include <windows.h>
#include <ctime>
#elif __linux__
#include <sys/sysinfo.h>
#include <unistd.h>
#elif __APPLE__
#include <sys/sysctl.h>
#include <unistd.h>
#endif // _WIN32

#include <iomanip>

using namespace std;

namespace Gedim
{
  // ***************************************************************************
#ifdef _WIN32
  string Output::BlueColor = "";
  string Output::RedColor = "";
  string Output::GreenColor = "";
  string Output::YellowColor = "";
  string Output::EndColor = "";
#elif __linux__
  string Output::BlueColor = "\033[1;34m";
  string Output::RedColor = "\033[1;31m";
  string Output::GreenColor = "\033[1;32m";
  string Output::YellowColor = "\033[1;33m";
  string Output::EndColor = "\033[0m";
#elif __APPLE__
  string Output::BlueColor = "";
  string Output::RedColor = "";
  string Output::GreenColor = "";
  string Output::YellowColor = "";
  string Output::EndColor = "";
#endif // _WIN32
  int Output::MaxElementToPrint = 1000;
  int Output::StartingIndexToPrint = 0;
  // ***************************************************************************
  void Output::CreateFolder(const string& nameFolder)
  {
    vector<string> folders = Output::StringSplit(nameFolder, '/');
    ostringstream nameFolderPath;

    for(unsigned int i = 0; i < folders.size(); i++)
    {
      nameFolderPath<< folders[i]<< "/";

      if(folders[i] == ".")
        continue;

#ifdef _WIN32
      mkdir(nameFolderPath.str().c_str());
#elif __linux__
      mkdir(nameFolderPath.str().c_str(), 0777);
#elif __APPLE__
      mkdir(nameFolderPath.str().c_str(), 0777);
#endif
    }
  }
  // *************************************************************************
  vector<string> Output::StringSplit(const string& stringToSplit, const char character)
  {
    stringstream splitter(stringToSplit);
    string tempString;
    vector<string> strings;

    while(getline(splitter, tempString, character))
      strings.push_back(tempString);

    return strings;
  }
  // ***************************************************************************
  bool Output::FileExists(const string& nameFile)
  {
    bool fileExists = false;

#ifdef _WIN32
    fileExists = (_access(nameFile.c_str(), 0) != -1);
#elif __linux__
    fileExists = (access(nameFile.c_str(), F_OK) != -1);
#elif __APPLE__
    fileExists = (access(nameFile.c_str(), F_OK) != -1);
#endif // _WIN32

    return fileExists;
  }
  // ***************************************************************************
  void Output::GetFilePath(const string& nameFilePath, string& filePath, string& fileName, string& fileExtension)
  {
    size_t foundPath = nameFilePath.find_last_of("/\\");
    size_t foundExtension = nameFilePath.find_last_of(".");

    if (foundPath != string::npos)
      filePath = nameFilePath.substr(0, foundPath + 1);

    if (foundPath != string::npos && foundExtension != string::npos)
      fileName = nameFilePath.substr(foundPath + 1, foundExtension - foundPath - 1);
    else if (foundPath != string::npos)
      fileName = nameFilePath.substr(foundPath + 1);
    else if (foundExtension != string::npos)
      fileName = nameFilePath.substr(0, foundExtension);
    else
      fileName = nameFilePath;

    if (foundExtension != string::npos && ((foundPath != string::npos && foundPath < foundExtension) || foundPath == string::npos) )
      fileExtension = nameFilePath.substr(foundExtension + 1);
  }
  // ***************************************************************************
  void Output::GetFileList(const string& mainDirectory, vector<string>& out, const FileFilter& filter, const bool& hiddenElements)
  {
    out.clear();

#ifdef WINDOWS
    HANDLE dir;
    WIN32_FIND_DATA file_data;

    if ((dir = FindFirstFile((mainDirectory + "/*").c_str(), &file_data)) == INVALID_HANDLE_VALUE)
      return; /* No files found */

    do {
      const string file_name = file_data.cFileName;
      const string full_file_name = mainDirectory + "/" + file_name;
      const bool is_directory = (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;

      if (file_name[0] == '.')
        continue;

      if (filter == DFNResult::Files && is_directory)
        continue;

      out.push_back(file_name);
    } while (FindNextFile(dir, &file_data));

    FindClose(dir);
#else
    DIR *dir;
    class dirent *ent;
    class stat st;

    dir = opendir(mainDirectory.c_str());

    if (dir == NULL)
    {
      Output::PrintErrorMessage("Cannot reach folder '%s'", false, mainDirectory.c_str());
      return;
    }

    while ((ent = readdir(dir)) != NULL) {
      const string file_name = ent->d_name;
      const string full_file_name = mainDirectory + "/" + file_name;

      if (!hiddenElements && file_name[0] == '.')
        continue;

      if (stat(full_file_name.c_str(), &st) == -1)
        continue;

      const bool is_directory = (st.st_mode & S_IFDIR) != 0;

      switch (filter)
      {
        case Output::Files:
          if (is_directory)
            continue;
          break;
        case Output::Directories:
          if (!is_directory)
            continue;
          break;
        default:
          break;
      }

      out.push_back(file_name);
    }
    closedir(dir);
#endif
  }
  // *****************************************************************
  void Output::FindPaths(const string& mainDirectory, const string& objectNameToFind, vector<string>& paths, const FileFilter& filter, const bool& hiddenElements)
  {
#ifdef WINDOWS

#else
    DIR *dir;
    class dirent *ent;
    class stat st;

    dir = opendir(mainDirectory.c_str());
    while ((ent = readdir(dir)) != NULL)
    {
      const string file_name = ent->d_name;
      const string full_file_name = mainDirectory + "/" + file_name;

      if (!hiddenElements && file_name[0] == '.')
        continue;

      if (stat(full_file_name.c_str(), &st) == -1)
        continue;

      const bool is_directory = (st.st_mode & S_IFDIR) != 0;

      if (is_directory)
        FindPaths(full_file_name, objectNameToFind, paths, filter, hiddenElements);

      switch (filter)
      {
        case Output::Files:
          if (is_directory)
            continue;
          break;
        case Output::Directories:
          if (!is_directory)
            continue;
          break;
        default:
          break;
      }

      if (file_name != objectNameToFind)
        continue;

      paths.push_back(full_file_name);
    }
    closedir(dir);
#endif
  }
  // ***************************************************************************
  void Output::GetBinaryFileSize(const string& nameFile,
                                 unsigned int& fileSize,
                                 const unsigned int& sizeOfSingleDataToWrite,
                                 const unsigned int& startingPosition)
  {
    /// <ul>

    /// <li> Open file
    ifstream file;
    file.open(nameFile.c_str(), ios::binary | ios::ate);

    if (file.fail())
      throw runtime_error("File '" + nameFile +"' cannot be opened");

    /// <li> Get file size
    fileSize = 0;

    fileSize += file.tellg();

    if (fileSize < (startingPosition * sizeOfSingleDataToWrite))
      throw runtime_error("Error in file '" + nameFile +
                          "': uncorrect starting position. Expected " +
                          to_string(startingPosition) +", obtained " +
                          to_string(fileSize / sizeof(double)));

    fileSize -= startingPosition * sizeOfSingleDataToWrite;

    file.close();

    /// </ul>
  }
  // ***************************************************************************
  void Output::ReadBinaryFile(const string& nameFile, void* dataToRead, const unsigned int& sizeOfSingleDataToWrite, const unsigned int& dataSizeToRead, const unsigned int& startingPosition)
  {
    /// <ul>

    /// <li> Get file dimension
    unsigned int fileSize = 0;
    GetBinaryFileSize(nameFile, fileSize, sizeOfSingleDataToWrite, startingPosition);

    if (fileSize == 0)
      return;

    /// <li> Open file
    ifstream file;
    file.open(nameFile.c_str(), ios::binary | ios::ate);
    if (file.fail())
      throw runtime_error("File '" + nameFile + "' cannot be opened");

    file.seekg(startingPosition * sizeOfSingleDataToWrite);

    if (fileSize < (dataSizeToRead * sizeOfSingleDataToWrite))
      throw runtime_error("Error in file '" + nameFile +
                          "': uncorrect size. Expected " +
                          to_string(dataSizeToRead) + ", obtained " +
                          to_string(fileSize / sizeof(double)));

    /// <li> Read from file
    unsigned int fileSizeToRead = dataSizeToRead == 0 ? fileSize : dataSizeToRead * sizeOfSingleDataToWrite;
    int fileContentSize = dataSizeToRead == 0 ? fileSize / sizeOfSingleDataToWrite : dataSizeToRead;

    if (fileContentSize == 0)
      throw runtime_error("Error in file '" + nameFile +
                          "': uncorrect content size. Obtained " +
                          to_string(fileContentSize));

    file.read((char*)dataToRead, fileSizeToRead);

    file.close();

    /// </ul>
  }
  // ***************************************************************************
  bool Output::ReadBinaryFile(const string& nameFile, vector<double>& dataToRead, const unsigned int& dataSizeToRead, const unsigned int& startingPosition)
  {
    /// <ul>

    /// <li> Get file dimension
    unsigned int fileSize = 0;

    if (dataSizeToRead == 0)
    {
      Output::GetBinaryFileSize(nameFile, fileSize, sizeof(double), startingPosition);

      if (fileSize == 0)
        return true;
    }

    int fileContentSize = dataSizeToRead == 0 ? fileSize / sizeof(double) : dataSizeToRead;

    if (fileContentSize == 0)
    {
      Output::PrintErrorMessage("Error in file '%s': uncorrect content size. Obtained %d", false, nameFile.c_str(), fileContentSize);
      return false;
    }

    dataToRead.resize(fileContentSize, 0.0);
    ReadBinaryFile(nameFile, dataToRead.data(), sizeof(double), fileContentSize, startingPosition);

    return true;
    /// </ul>
  }
  // ***************************************************************************
  void Output::WriteBinaryFile(const string& nameFile,
                               const void* dataToWrite,
                               const unsigned int& sizeOfSingleDataToWrite,
                               const unsigned int& dataSizeToWrite,
                               const bool& append)
  {
    /// <ul>

    if (dataToWrite == NULL)
      throw runtime_error("empty data to write");

    unsigned int fileDataSize = dataSizeToWrite;

    if (fileDataSize == 0)
      return;

    /// <li> Exporting file
    ofstream file;

    if (append)
      file.open(nameFile.c_str(), ios::binary | ofstream::app);
    else
      file.open(nameFile.c_str(), ios::binary);

    if(file.fail())
      throw runtime_error("File '" + nameFile + "' cannot be opened");

    file.write((char *)dataToWrite, fileDataSize * sizeOfSingleDataToWrite);
    file.close();

    /// </ul>
  }
  // ***************************************************************************
  bool Output::WriteBinaryFile(const string& nameFile,
                               const vector<double>& dataToWrite,
                               const unsigned int& dataSizeToWrite,
                               const unsigned int& dataStartingPositionToWrite,
                               const bool& append)
  {
    /// <ul>

    if (dataSizeToWrite > dataToWrite.size())
    {
      Output::PrintErrorMessage("Error in write file '%s': uncorrect size. Required %d, given %d", false, nameFile.c_str(), dataSizeToWrite, dataToWrite.size());
      return false;
    }

    if (dataStartingPositionToWrite > dataToWrite.size())
    {
      Output::PrintErrorMessage("Error in write file '%s': uncorrect starting position. Required %d, given %d", false, nameFile.c_str(), dataStartingPositionToWrite, dataToWrite.size());
      return false;
    }

    unsigned int fileDataSize = (dataSizeToWrite == 0) ? dataToWrite.size() : dataSizeToWrite;

    if (fileDataSize == 0)
      return true;

    /// <li> Exporting file
    WriteBinaryFile(nameFile, dataToWrite.data() + dataStartingPositionToWrite, sizeof(double), fileDataSize, append);

    return true;
    /// </ul>
  }
  // ***************************************************************************
  void Output::PrintLine(char symbol, bool onlyMaster)
  {
#if LOGGING > 0

    if (onlyMaster && MpiParallelEnvironment::Process().Rank() != 0)
      return;

#if LOGGING==1 || LOGGING==3
    cout << " >> ";
    for (int i = 1; i < 72; i++)
      cout<< symbol;
    cout<< endl;
#endif

#if LOGGING==2 || LOGGING==3
    LogFile::PrintLine(symbol);
#endif

#endif // LOGGING
  }
  // ***************************************************************************
  ostream& Output::PrintStars(ostream& out)
  {
    for (int i = 1; i < 76; i++)
      out<< '*';
    out<< endl;

    return out;
  }
  // **********************************
  ostream& Output::PrintLines(ostream& out)
  {
    for (int i = 1; i < 76; i++)
      out<< '-';
    out<< endl;

    return out;
  }
  // ***************************************************************************
  void Output::PrintStatusProgram(const string& message, ...)
  {
#if LOGGING > 0

    if (MpiParallelEnvironment::Process().Rank() != 0)
      return;

    va_list args, dupArgs;
    va_start(args, message);
    va_copy(dupArgs, args);

    // Printing
    ostringstream newMessage;
#if LOGGING==1 || LOGGING==3
    newMessage<< " >> "<< message<< " "<< Output::GreenColor<< "SUCCESS"<< Output::EndColor;
    newMessage<< endl;

    vprintf(newMessage.str().c_str(), args);
#endif

#if LOGGING==2 || LOGGING==3
    newMessage.str("");
    newMessage.clear();
    newMessage<< message<< " "<< "SUCCESS";
    LogFile::PrintInfoMessage(newMessage.str(), dupArgs);
#endif

    va_end(args);

#endif // LOGGING
  }
  // ***************************************************************************
  void Output::PrintGenericMessageOnLine(const string& message, ...)
  {
#if LOGGING > 0

    if (MpiParallelEnvironment::Process().Rank() != 0)
      return;

    va_list args;
    va_start(args, message);

    // Printing
#if LOGGING==1 || LOGGING==3
    ostringstream newMessage;
    newMessage<< "\r"<< message;
    newMessage<< flush;

    vprintf(newMessage.str().c_str(), args);
    cout<< flush;
#endif

    va_end(args);

#endif // LOGGING
  }
  // ***************************************************************************
  void Output::PrintErrorMessage(const string& message, const bool& onlyMaster, ...)
  {
#if VERBOSE >= 1
#if LOGGING > 0

    if (onlyMaster && MpiParallelEnvironment::Process().Rank() != 0)
      return;

    va_list args, dupArgs;
    va_start(args, onlyMaster);
    va_copy(dupArgs, args);

    ostringstream newMessage;
    newMessage<< Output::RedColor<< " >> "<< message<< Output::EndColor;
    if(!onlyMaster)
      newMessage<< " on process "<< MpiParallelEnvironment::Process().Rank();
    newMessage<< endl;

    vprintf(newMessage.str().c_str(), args);

#if LOGGING==2 || LOGGING==3
    LogFile::PrintErrorMessage(message, dupArgs);
#endif

    va_end(args);

#endif // LOGGING
#endif // VERBOSE
  }
  // ***************************************************************************
  void Output::PrintWarningMessage(const string& message, const bool& onlyMaster, ...)
  {
#if VERBOSE >= 2
#if LOGGING > 0

    if (onlyMaster && MpiParallelEnvironment::Process().Rank() != 0)
      return;

    va_list args, dupArgs;
    va_start(args, onlyMaster);
    va_copy(dupArgs, args);

    ostringstream newMessage;
    newMessage<< Output::YellowColor<< " >> "<< message<< Output::EndColor;
    if(!onlyMaster)
      newMessage<< " on process "<< MpiParallelEnvironment::Process().Rank();
    newMessage<< endl;

    vprintf(newMessage.str().c_str(), args);

#if LOGGING==2 || LOGGING==3
    LogFile::PrintWarningMessage(message, dupArgs);
#endif

    va_end(args);

#endif // LOGGING
#endif // VERBOSE
  }
  // ***************************************************************************
  void Output::PrintGenericMessage(const string& message, const bool& onlyMaster, ...)
  {
#if VERBOSE >= 3
#if LOGGING > 0

    if (onlyMaster && MpiParallelEnvironment::Process().Rank() != 0)
      return;

    va_list args, dupArgs;
    va_start(args, onlyMaster);
    va_copy(dupArgs, args);

    // Printing
#if LOGGING==1 || LOGGING==3
    ostringstream newMessage;
    newMessage<< " >> "<< message;
    if(!onlyMaster)
      newMessage<< " on process "<< MpiParallelEnvironment::Process().Rank();
    newMessage<< endl;

    vprintf(newMessage.str().c_str(), args);
#endif

#if LOGGING==2 || LOGGING==3
    LogFile::PrintInfoMessage(message, dupArgs);
#endif

    va_end(args);

#endif // LOGGING
#endif // VERBOSE
  }
  // ***************************************************************************
  void Output::PrintDebugMessage(const string& message, const bool& onlyMaster, ...)
  {
    Utilities::Unused(message); //< Variable may be unused depending on macro
    Utilities::Unused(onlyMaster); //< Variable may be unused depending on macro

#if VERBOSE >= 4
#if LOGGING > 0
    if (onlyMaster && MpiParallelEnvironment::Process().Rank() != 0)
      return;

    va_list args, dupArgs;
    va_start(args, onlyMaster);
    va_copy(dupArgs, args);

    ostringstream newMessage;
    newMessage<< Output::BlueColor<< " >> "<< message<< Output::EndColor;
    if(!onlyMaster)
      newMessage<< " on process "<< MpiParallelEnvironment::Process().Rank();
    newMessage<< endl;

    vprintf(newMessage.str().c_str(), args);

#if LOGGING==2 || LOGGING==3
    LogFile::PrintDebugMessage(message, dupArgs);
#endif

    va_end(args);

#endif // LOGGING
#endif // VERBOSE
  }
  // ***************************************************************************
  void Output::PrintSuccessMessage(const string& message, const bool& onlyMaster, ...)
  {
#if LOGGING > 0

    if (onlyMaster && MpiParallelEnvironment::Process().Rank() != 0)
      return;

    va_list args, dupArgs;
    va_start(args, onlyMaster);
    va_copy(dupArgs, args);

#if LOGGING==1 || LOGGING==3
    ostringstream newMessage;
    newMessage<< Output::GreenColor<< " >> "<< message<< Output::EndColor;
    if(!onlyMaster)
      newMessage<< " on process "<< MpiParallelEnvironment::Process().Rank();
    newMessage<< endl;

    vprintf(newMessage.str().c_str(), args);
#endif

#if LOGGING==2 || LOGGING==3
    LogFile::PrintInfoMessage(message, dupArgs);
#endif

    va_end(args);

#endif // LOGGING
  }
  // ***************************************************************************
  void Output::Assert(const bool& logicResult, const string& message, ...)
  {
    if (logicResult)
      return;

    if (MpiParallelEnvironment::Process().Rank() != 0 || message.empty())
      abort();

#if LOGGING > 0
    va_list args, dupArgs;
    va_start(args, message);
    va_copy(dupArgs, args);

    // Printing
    ostringstream newMessage;
#if LOGGING==1 || LOGGING==3
    newMessage<< " >> "<< message;
    newMessage<< " "<< Output::RedColor<< "FAILED"<< Output::EndColor;

    newMessage<< endl;

    vprintf(newMessage.str().c_str(), args);
#endif

#if LOGGING==2 || LOGGING==3
    newMessage.str("");
    newMessage.clear();
    newMessage<< message;
    newMessage<< " "<< "FAILED";

    LogFile::PrintInfoMessage(newMessage.str(), dupArgs);
#endif

    va_end(args);
#endif // LOGGING

    abort();
  }
  // ***************************************************************************
  int LogFile::LogMaxFileSize = 2;
  int LogFile::LogMaxNumFiles = 5;
  string LogFile::LogFolder = "Log";
  string LogFile::LogNameFile = "LogFile.log";
  // ***************************************************************************
  string LogFile::GetDateTime()
  {
    ostringstream timeBuilder;

    time_t now = time(0);   // get time now
    tm* ltm = localtime(&now);

    // print various components of tm structure.
    timeBuilder<< setfill('0') << setw(4)<< 1900 + ltm->tm_year<< "/";
    timeBuilder<< setfill('0') << setw(2)<< (1 + ltm->tm_mon)<< "/";
    timeBuilder<< setfill('0') << setw(2)<< ltm->tm_mday<< " ";
    timeBuilder<< setfill('0') << setw(2)<< ltm->tm_hour<< ":";
    timeBuilder<< setfill('0') << setw(2)<< ltm->tm_min<< ":";
    timeBuilder<< setfill('0') << setw(2)<< ltm->tm_sec;

    return timeBuilder.str();
  }
  // ***************************************************************************
  double LogFile::GetFileSize(const string& nameFile)
  {
    streampos begin,end;
    ifstream file(nameFile.c_str());
    begin = file.tellg();
    file.seekg (0, ios::end);
    end = file.tellg();
    file.close();
    double fileSize = (end - begin) / 100000.0; // MegaBytes

    return fileSize;
  }
  // ***************************************************************************
  void LogFile::CheckFileSize(const string& nameFile)
  {
    if(!Output::FileExists(nameFile))
      return;

    // Check file dimension
    double actualFileSize = GetFileSize(nameFile);
    if(actualFileSize > LogMaxFileSize)
    {
      // Create new log file
      ostringstream newNameFile, oldNameFile;
      int maxFiles = LogMaxNumFiles, numFiles = maxFiles;

      newNameFile<< nameFile<< "_"<< numFiles;
      if(Output::FileExists(newNameFile.str()))
        remove(newNameFile.str().c_str());

      do
      {
        if(numFiles > 1)
          oldNameFile<< nameFile<< "_"<< (numFiles - 1);
        else
          oldNameFile<< nameFile;

        if(Output::FileExists(oldNameFile.str()))
          rename(oldNameFile.str().c_str(), newNameFile.str().c_str());

        newNameFile.str("");
        newNameFile.clear();
        oldNameFile.str("");
        oldNameFile.clear();

        numFiles--;
        newNameFile<< nameFile<< "_"<< numFiles;
      }
      while(numFiles > 0);
    }
  }
  // ***************************************************************************
  void LogFile::PrintMessage(const string& type, const string& message, va_list args)
  {
    if (MpiParallelEnvironment::Process().Rank() != 0)
      return;

    ostringstream nameFolder, nameFile;

    nameFolder<< LogFolder<< "/";
    Output::CreateFolder(nameFolder.str());

    nameFile<< nameFolder.str()<< LogNameFile;
    CheckFileSize(nameFile.str());

    if(args != NULL)
    {
      FILE * pFile = fopen(nameFile.str().c_str(),"a");

      ostringstream logMessage;
      logMessage<< GetDateTime()<< " ";
      if(!type.empty())
        logMessage<< type<< ": ";
      logMessage<< message<< endl;

      vfprintf(pFile, logMessage.str().c_str(), args);
      va_end(args);

      fclose(pFile);
    }
    else
    {
      ofstream outFile(nameFile.str().c_str(), std::ios_base::app);

      outFile<< GetDateTime()<< " ";
      if(!type.empty())
        outFile<< type<< ": ";
      outFile<< message<< endl;

      outFile.close();
    }
  }
  // ***************************************************************************
  void LogFile::PrintLine(const char& symbol)
  {
    ostringstream separatorBuilder;

    for (int i = 1; i < 61; i++)
      separatorBuilder<< symbol;

    LogFile::PrintMessage("", separatorBuilder.str(), NULL);
  }
  // ***************************************************************************
  void LogFile::PrintWarningMessage(const string& message, va_list args)
  {
    LogFile::PrintMessage("WARN", message, args);
  }
  // ***************************************************************************
  void LogFile::PrintErrorMessage(const string& message, va_list args)
  {
    LogFile::PrintMessage("ERR ", message, args);
  }
  // ***************************************************************************
  void LogFile::PrintInfoMessage(const string& message, va_list args)
  {
    LogFile::PrintMessage("INFO", message, args);
  }
  // ***************************************************************************
  void LogFile::PrintDebugMessage(const string& message, va_list args)
  {
    LogFile::PrintMessage("DBG ", message, args);
  }
  // ***************************************************************************
  map<string, double> Profiler::times;
  map<string, double> Profiler::localTimes;
  unsigned long Profiler::totalAvailStartingMemory = 0;
  bool Profiler::ActivateProfiler = false;
  string Profiler::TimeFile = "Time";
  string Profiler::MemoryFile = "Memory";
  // ***************************************************************************
  double Profiler::GetTime(const string& nameTime, const bool& localTime)
  {
    if(!ActivateProfiler)
      return 0.0;

    double time = 0.0;

    if(localTime)
    {
      map<string, double>::iterator timeFound = localTimes.find(nameTime);
      if(timeFound != localTimes.end())
        time = timeFound->second;
    }
    else
    {
      map<string, double>::iterator timeFound = times.find(nameTime);
      if(timeFound != times.end())
        time = timeFound->second;
    }

    return time;
  }
  // ***************************************************************************
  void Profiler::Sleep(const unsigned int& milliseconds)
  {
#ifdef _WIN32
    Sleep(milliseconds);
#else
    usleep(milliseconds * 1000); // takes microseconds
#endif
  }
  // ***************************************************************************
  void Profiler::StartTime(const string& nameTime)
  {
    if (!ActivateProfiler)
      return;

    double time = 0.0;

    map<string, double>::iterator timeFound = times.find(nameTime);
    map<string, double>::iterator localTimeFound = localTimes.find(nameTime);

#if USE_MPI == 1
    MPI::COMM_WORLD.Barrier();
    time = MPI::Wtime();
    MPI::COMM_WORLD.Bcast(&time, 1, MPI::DOUBLE, 0);
#elif USE_MPI == 0
    time = clock();
#endif // USE_MPI

    if(timeFound != times.end())
      timeFound->second = time;
    else
      times.insert(pair<string, double>(nameTime, time));

    if(localTimeFound != localTimes.end())
      localTimeFound->second = time;
    else
      localTimes.insert(pair<string, double>(nameTime, time));
  }
  // *************************************************************************
  void Profiler::SplitTime(const string& nameTime, double& globalTime, double& localTime)
  {
    if (!ActivateProfiler)
      return;

    double endTime = 0.0, localEndTime = 0.0;

#if USE_MPI == 1
    localEndTime = MPI::Wtime();
    MPI::COMM_WORLD.Allreduce(&localEndTime, &endTime, 1, MPI::DOUBLE, MPI::MAX);
#elif USE_MPI == 0
    localEndTime = clock();
    endTime = localEndTime;
#endif // USE_MPI

    map<string, double>::iterator timeFound = times.find(nameTime);
    map<string, double>::iterator localTimeFound = localTimes.find(nameTime);
    if(timeFound == times.end())
    {
      Output::PrintErrorMessage(" -- Profiling: No time '" + nameTime + "' exists --", true);
      return;
    }
    if(localTimeFound == localTimes.end())
    {
      Output::PrintErrorMessage(" -- Profiling: No local time '" + nameTime + "' exists --", true);
      return;
    }

    globalTime = ComputeTime(timeFound->second, endTime);
    localTime = ComputeTime(localTimeFound->second, localEndTime);
  }
  // *************************************************************************
  double Profiler::StopTime(const string& nameTime, const bool& printTime)
  {
    if (!ActivateProfiler)
      return 0.0;

    double endTime = 0.0, localEndTime = 0.0;

#if USE_MPI == 1
    localEndTime = MPI::Wtime();
    MPI::COMM_WORLD.Allreduce(&localEndTime, &endTime, 1, MPI::DOUBLE, MPI::MAX);
#elif USE_MPI == 0
    localEndTime = clock();
    endTime = localEndTime;
#endif // USE_MPI

    map<string, double>::iterator timeFound = times.find(nameTime);
    map<string, double>::iterator localTimeFound = localTimes.find(nameTime);
    if(timeFound == times.end())
    {
      Output::PrintErrorMessage(" -- Profiling: No time '" + nameTime + "' exists --", true);
      return 0.0;
    }
    if(localTimeFound == localTimes.end())
    {
      Output::PrintErrorMessage(" -- Profiling: No local time '" + nameTime + "' exists --", true);
      return 0.0;
    }

    const double timeDifference = ComputeTime(timeFound->second, endTime);
    const double localTimeDifference = ComputeTime(localTimeFound->second, localEndTime);

    times[nameTime] = timeDifference;
    localTimes[nameTime] = localTimeDifference;

    if (!printTime)
      return timeDifference;

    WriteTime(nameTime,
              timeDifference,
              localTimeDifference);

    return timeDifference;
  }
  // *************************************************************************
  void Profiler::WriteTime(const std::string& nameTime,
                           const double& globalTime,
                           const double& localTime)
  {
    ostringstream nameFolder;

    nameFolder<< LogFile::LogFolder<< "/";
    Output::CreateFolder(nameFolder.str());

    // Print global time file
    if (MpiParallelEnvironment::Process().Rank() == 0)
    {
#if LOGGING==1 || LOGGING==3
      int hours = globalTime / 3600;
      int minutes = (globalTime - hours * 3600) / 60;
      int seconds = (globalTime - hours * 3600 - minutes * 60);
      double milliseconds = (globalTime - hours * 3600 - minutes * 60 - seconds) * 1000;
      cout<< " -- "<< Output::BlueColor<< "Profiling"<< Output::EndColor<< ": Time '"<< nameTime<< "': "<< hours<< " h "<< minutes<< " m "<< seconds<< " s "<< milliseconds<< " ms --"<< endl;
#endif

      ostringstream nameFile;

      nameFile<< nameFolder.str()<< TimeFile<< ".csv";
      bool fileExists = Output::FileExists(nameFile.str());

      LogFile::CheckFileSize(nameFile.str());

      ofstream outFile(nameFile.str().c_str(), ios_base::app);
      outFile.precision(8);

      if(!fileExists)
      {
        outFile<< "TotProc"<< " ";
        outFile<< "Type"<< " ";
        outFile<< "Time"<< endl;
      }

      outFile<< scientific<< MpiParallelEnvironment::Process().NumberProcesses()<< " "<< nameTime<< " "<< globalTime<< endl;
    }

    // Print local time files
    ostringstream nameFile;

    nameFolder<< "Processes/";
    Output::CreateFolder(nameFolder.str());

    nameFile<< nameFolder.str()<< TimeFile<< "_"<< MpiParallelEnvironment::Process().Rank()<< ".csv";
    bool fileExists = Output::FileExists(nameFile.str());

    LogFile::CheckFileSize(nameFile.str());

    ofstream outFile(nameFile.str().c_str(), ios_base::app);
    outFile.precision(8);

    if(!fileExists)
    {
      outFile<< "TotProc"<< " ";
      outFile<< "Type"<< " ";
      outFile<< "Time"<< endl;
    }

    outFile<< scientific<< MpiParallelEnvironment::Process().NumberProcesses()<< " "<< nameTime<< " "<< localTime<< endl;
  }
  // *************************************************************************
  int Profiler::ParseLine(char* line)
  {
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);

    return i;
  }
  // *************************************************************************
  void Profiler::GetProcessMemory(vector<unsigned long>& memoryUsed)
  {
    memoryUsed.resize(2, 0);

    // Get Virtual Memory
    string token;
    ifstream file("/proc/self/status");

    if (!file.is_open())
      return;

    while(file >> token)
    {
      if(token == "VmSize:")
        file >> memoryUsed[0];
      if(token == "VmRSS:")
        file >> memoryUsed[1];

      file.ignore(numeric_limits<streamsize>::max(), '\n');
    }
    file.close();
  }
  // *************************************************************************
  void Profiler::GetTotalMemory(vector<unsigned long>& totalMemory)
  {
    // NOTE: It is possible also to get Virtual memory from the file using the SwapTotal and SwapFree
    totalMemory.resize(3, 0);

    string token;
    ifstream file("/proc/meminfo");

    if (!file.is_open())
      return;

    while(file >> token)
    {
      if(token == "MemTotal:")
        file >> totalMemory[0];
      if(token == "MemAvailable:")
        file >> totalMemory[2];

      file.ignore(numeric_limits<streamsize>::max(), '\n');
    }
    file.close();

    // Get Used Physical Memory
    totalMemory[1] = totalMemory[0] - totalMemory[2];
  }
  // *************************************************************************
  void Profiler::CheckMemory(const string& nameMemory, const bool& checkProcessesMemory)
  {
    if(!ActivateProfiler)
      return;

#if USE_MPI == 1
    MPI::COMM_WORLD.Barrier();
#endif // USE_MPI

    // Get information on Virtual and Physical RAM
    vector<unsigned long> totalMemory;
    GetTotalMemory(totalMemory);
    vector<unsigned long> processMemory;
    GetProcessMemory(processMemory);

    vector<int> gigabytes(2, 0);
    vector<int> megabytes(2, 0);
    vector<int> kilobytes(2, 0);

    for(int i=0; i< 2; i++)
    {
      gigabytes[i] = (processMemory[i] / 1000000);
      megabytes[i] = (processMemory[i] - gigabytes[i] * 1000000) / 1000;
      kilobytes[i] = (processMemory[i] - gigabytes[i] * 1000000 - megabytes[i] * 1000);
    }

    // Communicate total memory of processes- Only for Physical RAM
    int totalProcessesMemoryGigabytes = 0;
    int totalProcessesMemoryMegabytes = 0;
    int totalProcessesMemoryKilobytes = 0;

#if USE_MPI == 1
    MPI::COMM_WORLD.Allreduce(&gigabytes[1], &totalProcessesMemoryGigabytes, 1, MPI::INT, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&megabytes[1], &totalProcessesMemoryMegabytes, 1, MPI::INT, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&kilobytes[1], &totalProcessesMemoryKilobytes, 1, MPI::INT, MPI::SUM);
#elif USE_MPI == 0
    totalProcessesMemoryGigabytes = gigabytes[1];
    totalProcessesMemoryMegabytes = megabytes[1];
    totalProcessesMemoryKilobytes = kilobytes[1];
#endif // USE_MPI

    unsigned long totalProcessMemory = (unsigned long)totalProcessesMemoryGigabytes * 1000000 + (unsigned long)totalProcessesMemoryMegabytes * 1000 + (unsigned long)totalProcessesMemoryKilobytes;
    totalProcessesMemoryGigabytes = totalProcessMemory / 1000000;
    totalProcessesMemoryMegabytes = (totalProcessMemory - totalProcessesMemoryGigabytes * 1000000) / 1000;
    totalProcessesMemoryKilobytes = (totalProcessMemory - totalProcessesMemoryGigabytes * 1000000 - totalProcessesMemoryMegabytes * 1000);

    if(totalAvailStartingMemory == 0)
      totalAvailStartingMemory = totalMemory[2] + totalProcessMemory;
    unsigned long totalAvailMemory = totalMemory[2];

#if LOGGING==1 || LOGGING==3
    int availableMemoryGygabytes = 0;
    int availableMemoryMegabytes = 0;
    int availableMemoryKilobytes = 0;

    availableMemoryGygabytes = totalAvailMemory / 1000000;
    availableMemoryMegabytes = (totalAvailMemory - availableMemoryGygabytes * 1000000) / 1000;
    availableMemoryKilobytes = (totalAvailMemory - availableMemoryGygabytes * 1000000 - availableMemoryMegabytes * 1000);
#endif // LOGGING

#if USE_MPI == 1
    MPI::COMM_WORLD.Barrier();
#endif // USE_MPI

    // Get Information on PHYSICAL RAM
    if(MpiParallelEnvironment::Process().Rank() == 0)
    {
#if LOGGING==1 || LOGGING==3
      if(checkProcessesMemory)
      {
        cout<< " -- "<< Output::BlueColor<< "Profiling"<< Output::EndColor<< ": RAM Memory '"<< nameMemory<<
               "': Total"<< "      "<< " Used / Free: "<<
               setfill(' ') << setw(5)<< totalProcessesMemoryGigabytes<< " / "<<
               setfill(' ') << setw(5)<< availableMemoryGygabytes<< " GB "<<
               setfill(' ') << setw(5)<< totalProcessesMemoryMegabytes<< " / "<<
               setfill(' ') << setw(5)<< availableMemoryMegabytes<< " MB "<<
               setfill(' ') << setw(5)<< totalProcessesMemoryKilobytes<< " / "<<
               setfill(' ') << setw(5)<< availableMemoryKilobytes<< " KB "<< " --"<< endl;
      }
      else
      {
        cout<< " -- "<< Output::BlueColor<< "Profiling"<< Output::EndColor<< ": RAM Memory '"<< nameMemory<<
               "': Used / Free: "<<
               totalProcessesMemoryGigabytes<< " / "<<
               availableMemoryGygabytes<< " GB "<<
               totalProcessesMemoryMegabytes<< " / "<<
               availableMemoryMegabytes<< " MB "<<
               totalProcessesMemoryKilobytes<< " / "<<
               availableMemoryKilobytes<< " KB "<< " --"<< endl;
      }
#endif

      ostringstream nameFile, nameFolder;

      nameFolder<< LogFile::LogFolder<< "/";
      Output::CreateFolder(nameFolder.str());

      nameFile<< nameFolder.str()<< MemoryFile<< ".csv";
      bool fileExists = Output::FileExists(nameFile.str());

      LogFile::CheckFileSize(nameFile.str());

      ofstream outFile(nameFile.str().c_str(), ios_base::app);
      outFile.precision(8);

      if(!fileExists)
      {
        outFile<< "TotProc"<< " ";
        outFile<< "Type"<< " ";
        outFile<< "MemoryUsed(KB)"<< " ";
        outFile<< "MemoryAvail(KB)"<< " ";
        outFile<< "PercUsed(%)"<< endl;
      }
      outFile<< MpiParallelEnvironment::Process().NumberProcesses()<< " "<<
                nameMemory<< " "<<
                totalProcessMemory<< " "<<
                totalAvailMemory<< " "<<
                (totalProcessMemory == 0 ? 0 : (unsigned int)(1.0 / (1.0 + totalAvailStartingMemory / totalProcessMemory) * 100.0))<< endl;
    }
#if USE_MPI == 1
    MPI::COMM_WORLD.Barrier();
#endif // USE_MPI

    if(!checkProcessesMemory)
      return;

#if LOGGING==1 || LOGGING==3
    cout<< " -- "<< Output::BlueColor<< "Profiling"<< Output::EndColor<< ": RAM Memory '"<< nameMemory<<
           "': Process "<< setfill(' ') << setw(3)<< MpiParallelEnvironment::Process().Rank() << " Used / Free: "<<
           setfill(' ') << setw(5)<< gigabytes[1]<< " / "<<
           setfill(' ') << setw(5)<< availableMemoryGygabytes<< " GB "<<
           setfill(' ') << setw(5)<< megabytes[1]<< " / "<<
           setfill(' ') << setw(5)<< availableMemoryMegabytes<< " MB "<<
           setfill(' ') << setw(5)<< kilobytes[1]<< " / "<<
           setfill(' ') << setw(5)<< availableMemoryKilobytes<< " KB "<< " --"<< endl;
#endif

#if USE_MPI == 1
    MPI::COMM_WORLD.Barrier();
#endif // USE_MPI
  }
  // ***************************************************************************
}
