#ifndef __CONFIGURATIONS_H
#define __CONFIGURATIONS_H

#include <string>
#include <unordered_map>

#include "ConfigurationPropertySupportedTypes.hpp"

#include "StringsUtilities.hpp"

using namespace std;

namespace Gedim
{
class Configurations final
{
  private:
    template <class T> struct ConfigurationProperty
    {
        string Id;
        string Description = "";
        T Value = T();
    };

    struct ExportProperty
    {
        string Description = "";
        string Value;
    };

    static unsigned int numberProperties;
    static unordered_map<string, ConfigurationPropertySupportedTypes::SupportedTypes> propertyTypes;
    static unordered_map<string, void *> properties;

    /// Remove Single Property by Id and Type
    static void RemoveProperty(const string &id, const ConfigurationPropertySupportedTypes::SupportedTypes &type);

    /// Get Property Type by Id
    static ConfigurationPropertySupportedTypes::SupportedTypes &GetPropertyType(const string &id)
    {
        return propertyTypes[id];
    }

    /// Get Property by Id
    template <class T> static ConfigurationProperty<T> &GetProperty(const string &id)
    {
        return *static_cast<ConfigurationProperty<T> *>(properties[id]);
    }

    /// Get Property Value To string
    static ExportProperty GetPropertyForExport(const string &id, const ConfigurationPropertySupportedTypes::SupportedTypes &type);

  public:
    ~Configurations()
    {
        Reset();
    }

    static const unsigned int &NumberProperties()
    {
        return numberProperties;
    }

    /// Read the configuration from csv file
    static void InitializeFromCsv(const string &ConfigurationsFile, const bool &hasHeader = true, const char &separator = ';');

    /// Read the configuration from ini file
    static void InitializeFromIni(const string &ConfigurationsFile);

    /// Read the configuration from command arguments
    /// @note format shall contain keyword id + char + type + char + value (example: id:type=value)
    static void Initialize(const int &argc, char **argv, const string &format = "id:type=value");

    /// Export the configuration to csv file
    static void ExportToCsv(const string &filePath, const bool &append = false, const char &separator = ';');

    /// Export the configuration to ini file
    static void ExportToIni(const string &filePath, const bool &append = false, const string &section = "");

    /// Reset the Configurations class
    static void Reset();

    /// Remove Single Property
    static void RemoveProperty(const string &id);

    /// Check if Property exists
    static bool ExistsProperty(const string &id)
    {
        return !(propertyTypes.find(id) == propertyTypes.end());
    }

    template <class T> static bool CheckPropertyType(const string &id)
    {
        return GetPropertyType(id) == ConfigurationPropertySupportedTypes::GetSupportedType<T>();
    }

    /// Convert Property by string Type and string value
    static void ConvertPropertyFromString(const string &id, const string &type, const string &value, const string &description = "");

    /// Add or Overwrite Property
    template <class T> static void AddProperty(const string &id, const T &value = T(), const string &description = "")
    {
        if (!ConfigurationPropertySupportedTypes::IsSupported<T>())
            throw invalid_argument("Property '" + id + "': type not supported");

        RemoveProperty(id);

        properties.insert(pair<string, void *>(id, new ConfigurationProperty<T>()));
        propertyTypes.insert(pair<string, ConfigurationPropertySupportedTypes::SupportedTypes>(
            id,
            ConfigurationPropertySupportedTypes::GetSupportedType<T>()));
        numberProperties++;

        ConfigurationProperty<T> &property = GetProperty<T>(id);

        property.Id = id;
        property.Value = value;
        property.Description = description;
    }

    /// Add or Overwrite Property converting const char* to string
    static void AddProperty(const string &id, const char *value = "", const string &description = "");

    /// Set Property Value
    template <class T> static void SetPropertyValue(const string &id, const T value)
    {
        if (!ExistsProperty(id))
            throw invalid_argument("");

        if (!CheckPropertyType<T>(id))
        {
            string propertyType = ConfigurationPropertySupportedTypes::TypeToString(GetPropertyType(id));
            string valueType = ConfigurationPropertySupportedTypes::TypeToString<T>();
            throw invalid_argument("Property '" + id + "' has type " + propertyType + " and not " + valueType);
        }

        ConfigurationProperty<T> &property = GetProperty<T>(id);
        property.Value = value;
    }

    /// Set Property Value converting const char* to string
    static void SetPropertyValue(const string &id, const char *value);

    /// Get Property Value
    template <class T> static const T &GetPropertyValue(const string &id)
    {
        if (!ExistsProperty(id))
            throw invalid_argument("");

        if (!CheckPropertyType<T>(id))
        {
            string propertyType = ConfigurationPropertySupportedTypes::TypeToString(GetPropertyType(id));
            string valueType = ConfigurationPropertySupportedTypes::TypeToString<T>();
            throw invalid_argument("Property '" + id + "' has type " + propertyType + " and not " + valueType);
        }

        const ConfigurationProperty<T> &property = GetProperty<T>(id);
        return property.Value;
    }
};
} // namespace Gedim

#endif // __CONFIGURATIONS_H
