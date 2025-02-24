#ifndef __ImportExportUtilities_H
#define __ImportExportUtilities_H

#include "Eigen/Eigen"

namespace Gedim_ImportExport_Utilities
{

template <class T, int Rows, int Cols>
std::ostream &operator<<(std::ostream &out, const Eigen::Matrix<T, Rows, Cols> &matrix)
{
    out << matrix.rows() << "," << matrix.cols() << ",";
    out << "{";
    for (unsigned int c = 0; c < matrix.cols(); c++)
    {
        out << (c != 0 ? ",{" : "{");
        for (unsigned int r = 0; r < matrix.rows(); r++)
            out << (r != 0 ? ", " : "") << matrix(r, c);
        out << "}";
    }
    out << "}";

    return out;
}

template <class T, std::size_t s> std::ostream &operator<<(std::ostream &out, const std::array<T, s> &elements)
{
    out << elements.size() << ",";
    out << "{";
    unsigned int i = 0;
    for (const auto &element : elements)
    {
        out << (i != 0 ? "," : "") << element;
        i++;
    }
    out << "}";

    return out;
}

template <class T> std::ostream &operator<<(std::ostream &out, const std::vector<T> &elements)
{
    out << elements.size() << ",";
    out << "{";
    unsigned int i = 0;
    for (const auto &element : elements)
    {
        out << (i != 0 ? "," : "") << element;
        i++;
    }
    out << "}";

    return out;
}

template <class T, int Rows, int Cols> std::istream &operator>>(std::istream &in, Eigen::Matrix<T, Rows, Cols> &matrix)
{
    char separator;
    unsigned int rows = 0, cols = 0;
    in >> rows;
    in >> separator;
    in >> cols;
    matrix.resize(rows, cols);
    in >> separator;

    in >> separator;
    for (unsigned int c = 0; c < matrix.cols(); c++)
    {
        in >> separator;
        for (unsigned int r = 0; r < matrix.rows(); r++)
        {
            in >> matrix(r, c);
            in >> separator;
        }
        in >> separator;
    }

    return in;
}

template <class T, std::size_t s> std::istream &operator>>(std::istream &in, std::array<T, s> &elements)
{
    char separator;
    unsigned int size = 0;
    in >> size;
    in >> separator;

    in >> separator;
    for (unsigned int v = 0; v < size; v++)
    {
        T element;
        in >> element;
        in >> separator;
        elements[v] = element;
    }

    return in;
}

template <class T> std::istream &operator>>(std::istream &in, std::vector<T> &elements)
{
    char separator;
    unsigned int size = 0;
    in >> size;
    elements.resize(size);
    in >> separator;

    in >> separator;
    for (unsigned int v = 0; v < size; v++)
    {
        T element;
        in >> element;
        in >> separator;
        elements[v] = element;
    }

    return in;
}
} // namespace Gedim_ImportExport_Utilities

#endif // __IOUtilities_H
