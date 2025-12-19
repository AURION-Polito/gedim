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

#ifndef __GEDIM_IOStream_H
#define __GEDIM_IOStream_H

#include "Gedim_Macro.hpp"
#include <unordered_set>

#if USE_MPI == 1
#include <mpi.h>
#endif

#include "IOUtilities.hpp"
#include <cstdarg>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace Gedim
{
template <typename matrixType>
std::string MatrixToString(const matrixType &mat, const std::string &matrixTypeStr, const std::string &matrixName)
{
    std::ostringstream str;
    str.precision(16);
    str << matrixName << " = " << matrixTypeStr << "(" << mat.rows() << ", " << mat.cols() << ");" << std::endl;
    for (unsigned int c = 0; c < mat.cols(); c++)
    {
        str << matrixName << ".col(" << c << ")<< ";
        for (unsigned int r = 0; r < mat.rows(); r++)
            str << std::scientific << (r == 0 ? "" : ",") << mat(r, c);
        str << ";" << std::endl;
    }

    return str.str();
}

template <typename matrixType>
std::string MatrixCollectionToString(const std::vector<matrixType> &matCollection,
                                     const std::string &matrixTypeStr,
                                     const std::string &matrixName)
{
    std::ostringstream str;
    str.precision(16);

    str << matrixName << " = std::vector<" << matrixTypeStr << ">(" << matCollection.size() << ");" << std::endl;
    for (unsigned int v = 0; v < matCollection.size(); v++)
    {
        str << MatrixToString<matrixType>(matCollection.at(v), matrixTypeStr, matrixName + "[" + std::to_string(v) + "]");
    }

    return str.str();
}

template <typename matrixType>
std::string MatrixCollectionToString(const std::vector<std::vector<matrixType>> &matCollection,
                                     const std::string &matrixTypeStr,
                                     const std::string &matrixName)
{
    std::ostringstream str;
    str.precision(16);

    str << matrixName << " = std::vector<std::vector<" << matrixTypeStr << ">>(" << matCollection.size() << ");" << std::endl;
    for (unsigned int v = 0; v < matCollection.size(); v++)
    {
        str << MatrixCollectionToString<matrixType>(matCollection.at(v), matrixTypeStr, matrixName + "[" + std::to_string(v) + "]");
    }

    return str.str();
}

template <typename matrixType>
std::string MatrixCollectionToString(const std::vector<std::vector<std::vector<matrixType>>> &matCollection,
                                     const std::string &matrixTypeStr,
                                     const std::string &matrixName)
{
    std::ostringstream str;
    str.precision(16);

    str << matrixName << " = std::vector<std::vector<std::vector<" << matrixTypeStr << ">>>(" << matCollection.size()
        << ");" << std::endl;
    for (unsigned int v = 0; v < matCollection.size(); v++)
    {
        str << MatrixCollectionToString<matrixType>(matCollection.at(v), matrixTypeStr, matrixName + "[" + std::to_string(v) + "]");
    }

    return str.str();
}

/// General print of a vector
template <typename T> std::ostream &operator<<(std::ostream &out, const std::vector<T> &vec)
{
    const int &maxElementToPrint = Output::MaxElementToPrint;
    const int &startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeVector = (maxElementToPrint == 0 || maxElementToPrint > (int)vec.size()) ? vec.size() : maxElementToPrint;
    unsigned int startIndexVector = (startingIndex >= (int)vec.size()) ? 0 : startingIndex;

    out << "[";
    for (unsigned int i = startIndexVector; i < startIndexVector + sizeVector && i < vec.size(); i++)
        out << (i != startIndexVector ? "," : "") << vec.at(i);
    out << "]";

    return out;
}
/// General print of a vector
template <typename T> std::ostream &operator<<(std::ostream &out, const std::vector<T *> &vec)
{
    const int &maxElementToPrint = Output::MaxElementToPrint;
    const int &startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeVector = (maxElementToPrint == 0 || maxElementToPrint > (int)vec.size()) ? vec.size() : maxElementToPrint;
    unsigned int startIndexVector = (startingIndex >= (int)vec.size()) ? 0 : startingIndex;

    out << "[";
    for (unsigned int i = startIndexVector; i < startIndexVector + sizeVector && i < vec.size(); i++)
        out << (i != startIndexVector ? "," : "") << *vec.at(i);
    out << "]";

    return out;
}
/// General print of an array
template <typename T, size_t s> std::ostream &operator<<(std::ostream &out, const std::array<T, s> &vec)
{
    const int &maxElementToPrint = Output::MaxElementToPrint;
    const int &startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeVector = (maxElementToPrint == 0 || maxElementToPrint > (int)vec.size()) ? vec.size() : maxElementToPrint;
    unsigned int startIndexVector = (startingIndex >= (int)vec.size()) ? 0 : startingIndex;

    out << "[";
    for (unsigned int i = startIndexVector; i < startIndexVector + sizeVector && i < vec.size(); i++)
        out << (i != startIndexVector ? "," : "") << vec.at(i);
    out << "]";

    return out;
}
/// General print of an array
template <typename T, size_t s> std::ostream &operator<<(std::ostream &out, const std::array<T *, s> &vec)
{
    const int &maxElementToPrint = Output::MaxElementToPrint;
    const int &startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeVector = (maxElementToPrint == 0 || maxElementToPrint > (int)vec.size()) ? vec.size() : maxElementToPrint;
    unsigned int startIndexVector = (startingIndex >= (int)vec.size()) ? 0 : startingIndex;

    out << "[";
    for (unsigned int i = startIndexVector; i < startIndexVector + sizeVector && i < vec.size(); i++)
        out << (i != startIndexVector ? "," : "") << *vec.at(i);
    out << "]";

    return out;
}
/// General print of a list
template <typename T> std::ostream &operator<<(std::ostream &out, const std::list<T> &listToPrint)
{
    const int &maxElementToPrint = Output::MaxElementToPrint;
    const int &startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)listToPrint.size()) ? listToPrint.size()
                                                                                                    : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)listToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out << "{";
    for (typename std::list<T>::const_iterator iterator = listToPrint.begin(); iterator != listToPrint.end(); iterator++)
    {
        if (counter >= startIndexList)
        {
            out << (counter != startIndexList ? "," : "") << *iterator;
            elementPrinted++;
        }

        if (elementPrinted >= sizeList)
            break;

        counter++;
    }
    out << "}";

    return out;
}
/// General print of a list
template <typename T> std::ostream &operator<<(std::ostream &out, const std::list<T *> &listToPrint)
{
    const int &maxElementToPrint = Output::MaxElementToPrint;
    const int &startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)listToPrint.size()) ? listToPrint.size()
                                                                                                    : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)listToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out << "{";
    for (typename std::list<T *>::const_iterator iterator = listToPrint.begin(); iterator != listToPrint.end(); iterator++)
    {
        if (counter >= startIndexList)
        {
            out << (counter != startIndexList ? "," : "") << **iterator;
            elementPrinted++;
        }

        if (elementPrinted >= sizeList)
            break;

        counter++;
    }
    out << "}";

    return out;
}
/// General print of a unordered_set
template <typename T> std::ostream &operator<<(std::ostream &out, const std::unordered_set<T> &setToPrint)
{
    const int &maxElementToPrint = Output::MaxElementToPrint;
    const int &startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)setToPrint.size()) ? setToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)setToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out << "{";
    for (typename std::unordered_set<T>::const_iterator iterator = setToPrint.begin(); iterator != setToPrint.end(); iterator++)
    {
        if (counter >= startIndexList)
        {
            out << (counter != startIndexList ? "," : "") << *iterator;
            elementPrinted++;
        }

        if (elementPrinted >= sizeList)
            break;

        counter++;
    }
    out << "}";

    return out;
}
/// General print of a unordered_set
template <typename T> std::ostream &operator<<(std::ostream &out, const std::unordered_set<T *> &setToPrint)
{
    const int &maxElementToPrint = Output::MaxElementToPrint;
    const int &startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)setToPrint.size()) ? setToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)setToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out << "{";
    for (typename std::unordered_set<T *>::const_iterator iterator = setToPrint.begin(); iterator != setToPrint.end(); iterator++)
    {
        if (counter >= startIndexList)
        {
            out << (counter != startIndexList ? "," : "") << **iterator;
            elementPrinted++;
        }

        if (elementPrinted >= sizeList)
            break;

        counter++;
    }
    out << "}";

    return out;
}
/// General print of a set
template <typename T> std::ostream &operator<<(std::ostream &out, const std::set<T> &setToPrint)
{
    const int &maxElementToPrint = Output::MaxElementToPrint;
    const int &startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)setToPrint.size()) ? setToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)setToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out << "{";
    for (typename std::set<T>::const_iterator iterator = setToPrint.begin(); iterator != setToPrint.end(); iterator++)
    {
        if (counter >= startIndexList)
        {
            out << (counter != startIndexList ? "," : "") << *iterator;
            elementPrinted++;
        }

        if (elementPrinted >= sizeList)
            break;

        counter++;
    }
    out << "}";

    return out;
}
/// General print of a set
template <typename T> std::ostream &operator<<(std::ostream &out, const std::set<T *> &setToPrint)
{
    const int &maxElementToPrint = Output::MaxElementToPrint;
    const int &startingIndex = Output::StartingIndexToPrint;

    unsigned int sizeList = (maxElementToPrint == 0 || maxElementToPrint > (int)setToPrint.size()) ? setToPrint.size() : maxElementToPrint;
    unsigned int startIndexList = (startingIndex >= (int)setToPrint.size()) ? 0 : startingIndex;
    unsigned int counter = 0, elementPrinted = 0;

    out << "{";
    for (typename std::set<T *>::const_iterator iterator = setToPrint.begin(); iterator != setToPrint.end(); iterator++)
    {
        if (counter >= startIndexList)
        {
            out << (counter != startIndexList ? "," : "") << **iterator;
            elementPrinted++;
        }

        if (elementPrinted >= sizeList)
            break;

        counter++;
    }
    out << "}";

    return out;
}
/// General print of a map
template <typename map_key, typename map_val>
std::ostream &operator<<(std::ostream &out, const std::map<map_key, map_val> &mapToPrint)
{
    unsigned int counter = 0;

    out << "{";
    for (typename std::map<map_key, map_val>::const_iterator it = mapToPrint.begin(); it != mapToPrint.end(); ++it)
    {
        out << (counter != 0 ? "," : "") << "{" << it->first << ", " << it->second << "}";
        counter++;
    }
    out << "}";

    return out;
}
/// General print of a map
template <typename map_key, typename map_val>
std::ostream &operator<<(std::ostream &out, const std::map<map_key, map_val *> &mapToPrint)
{
    unsigned int counter = 0;

    out << "{";
    for (typename std::map<map_key, map_val *>::const_iterator it = mapToPrint.begin(); it != mapToPrint.end(); ++it)
    {
        out << (counter != 0 ? "," : "") << "{" << it->first << ", " << (it->second == nullptr ? "NULL" : (*it->second)) << "}";
        counter++;
    }
    out << "}";

    return out;
}
/// General print of a unordered map
template <typename map_key, typename map_val>
std::ostream &operator<<(std::ostream &out, const std::unordered_map<map_key, map_val> &mapToPrint)
{
    unsigned int counter = 0;

    out << "{";
    for (typename std::unordered_map<map_key, map_val>::const_iterator it = mapToPrint.begin(); it != mapToPrint.end(); ++it)
    {
        out << (counter != 0 ? "," : "") << "{" << it->first << ", " << it->second << "}";
        counter++;
    }
    out << "}";

    return out;
}

} // namespace Gedim

#endif
