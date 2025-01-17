#ifndef __CONFIGURATIONPROPERTYSUPPORTEDTYPES_H
#define __CONFIGURATIONPROPERTYSUPPORTEDTYPES_H

#include "Eigen/Eigen"
#include "StringsUtilities.hpp"
#include <set>
#include <string>
#include <typeinfo>

namespace Gedim
{
class ConfigurationPropertySupportedTypes
{
  public:
    enum SupportedTypes
    {
        Unknown = 0,
        Int = 1,
        Char = 2,
        String = 3,
        Double = 4,
        VectorInt = 5,
        VectorDouble = 6,
        Bool = 7,
        SetUInt = 8,
        UInt = 9,
        VectorUInt = 10,
        Matrix = 11
    };

    /// Get Supported Type from the type
    template <class K> static SupportedTypes GetSupportedType()
    {
        const char *propertyType = typeid(K).name();

        if (propertyType == typeid(int).name())
            return Int;

        if (propertyType == typeid(char).name())
            return Char;

        if (propertyType == typeid(double).name())
            return Double;

        if (propertyType == typeid(std::string).name())
            return String;

        if (propertyType == typeid(std::vector<int>).name())
            return VectorInt;

        if (propertyType == typeid(std::vector<double>).name())
            return VectorDouble;

        if (propertyType == typeid(bool).name())
            return Bool;

        if (propertyType == typeid(unsigned int).name())
            return UInt;

        if (propertyType == typeid(std::vector<unsigned int>).name())
            return VectorUInt;

        if (propertyType == typeid(std::set<unsigned int>).name())
            return SetUInt;

        if (propertyType == typeid(Eigen::MatrixXd).name())
            return Matrix;

        return Unknown;
    }

    template <class K> static std::string TypeToString()
    {
        return TypeToString(GetSupportedType<K>());
    }

    /// Convert the SupportedType to string
    static std::string TypeToString(const SupportedTypes &type)
    {
        switch (type)
        {
        case Int:
            return "int";
        case Char:
            return "char";
        case String:
            return "string";
        case Double:
            return "double";
        case VectorInt:
            return "vector_int";
        case VectorDouble:
            return "vector_double";
        case Bool:
            return "bool";
        case UInt:
            return "uint";
        case SetUInt:
            return "set_uint";
        case VectorUInt:
            return "vector_uint";
        case Matrix:
            return "Matrix";
        default:
            return "unknown";
        }
    }

    /// Convert a string to SupportedType
    static SupportedTypes StringToType(const std::string &stringType)
    {
        std::string stringTypeToLow = StringsUtilities::ToLower(stringType);

        if (stringTypeToLow == "int")
            return Int;
        if (stringTypeToLow == "char")
            return Char;
        if (stringTypeToLow == "string")
            return String;
        if (stringTypeToLow == "double")
            return Double;
        if (stringTypeToLow == "vector_int")
            return VectorInt;
        if (stringTypeToLow == "vector_double")
            return VectorDouble;
        if (stringTypeToLow == "bool")
            return Bool;
        if (stringTypeToLow == "set_uint")
            return SetUInt;
        if (stringTypeToLow == "uint")
            return UInt;
        if (stringTypeToLow == "vector_uint")
            return VectorUInt;
        if (stringTypeToLow == "Matrix")
            return Matrix;

        return Unknown;
    }

    template <class K> static bool IsSupported()
    {
        return GetSupportedType<K>() != ConfigurationPropertySupportedTypes::Unknown;
    }
};
} // namespace Gedim

#endif // __CONFIGURATIONPROPERTYSUPPORTEDTYPES_H
