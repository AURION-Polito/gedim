#ifndef __CONFIGURATIONS_H
#define __CONFIGURATIONS_H

#include <string>
#include <unordered_map>

#include "ConfigurationPropertySupportedTypes.hpp"

#include "StringsUtilities.hpp"

namespace Gedim
{
class Configurations final
{
  private:
    template <class T> struct ConfigurationProperty
    {
        std::string Id;
        std::string Description = "";
        T Value = T();
    };

    struct ExportProperty
    {
        std::string Description = "";
        std::string Value;
    };

    static unsigned int numberProperties;
    static std::unordered_map<std::string, ConfigurationPropertySupportedTypes::SupportedTypes> propertyTypes;
    static std::unordered_map<std::string, void *> properties;

    /// Remove Single Property by Id and Type
    static void RemoveProperty(const std::string &id, const ConfigurationPropertySupportedTypes::SupportedTypes &type);

    /// Get Property Type by Id
    static ConfigurationPropertySupportedTypes::SupportedTypes &GetPropertyType(const std::string &id)
    {
        return propertyTypes[id];
    }

    /// Get Property by Id
    template <class T> static ConfigurationProperty<T> &GetProperty(const std::string &id)
    {
        return *static_cast<ConfigurationProperty<T> *>(properties[id]);
    }

    /// Get Property Value To string
    static ExportProperty GetPropertyForExport(const std::string &id,
                                               const ConfigurationPropertySupportedTypes::SupportedTypes &type);

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
    static void InitializeFromCsv(const std::string &ConfigurationsFile, const bool &hasHeader = true, const char &separator = ';');

    /// Read the configuration from ini file
    static void InitializeFromIni(const std::string &ConfigurationsFile);

    /// Read the configuration from command arguments
    /// @note format shall contain keyword id + char + type + char + value (example: id:type=value)
    static void Initialize(const int &argc, char **argv, const std::string &format = "id:type=value");

    /// Export the configuration to csv file
    static void ExportToCsv(const std::string &filePath, const bool &append = false, const char &separator = ';');

    /// Export the configuration to ini file
    static void ExportToIni(const std::string &filePath, const bool &append = false, const std::string &section = "");

    /// Reset the Configurations class
    static void Reset();

    /// Remove Single Property
    static void RemoveProperty(const std::string &id);

    /// Check if Property exists
    static bool ExistsProperty(const std::string &id)
    {
        return !(propertyTypes.find(id) == propertyTypes.end());
    }

    template <class T> static bool CheckPropertyType(const std::string &id)
    {
        return GetPropertyType(id) == ConfigurationPropertySupportedTypes::GetSupportedType<T>();
    }

    /// Convert Property by string Type and string value
    static void ConvertPropertyFromString(const std::string &id,
                                          const std::string &type,
                                          const std::string &value,
                                          const std::string &description = "");

    /// Add or Overwrite Property
    template <class T>
    static void AddProperty(const std::string &id, const T &value = T(), const std::string &description = "")
    {
        if (!ConfigurationPropertySupportedTypes::IsSupported<T>())
            throw std::invalid_argument("Property '" + id + "': type not supported");

        RemoveProperty(id);

        properties.insert(std::pair<std::string, void *>(id, new ConfigurationProperty<T>()));
        propertyTypes.insert(std::pair<std::string, ConfigurationPropertySupportedTypes::SupportedTypes>(
            id,
            ConfigurationPropertySupportedTypes::GetSupportedType<T>()));
        numberProperties++;

        ConfigurationProperty<T> &property = GetProperty<T>(id);

        property.Id = id;
        property.Value = value;
        property.Description = description;
    }

    /// Add or Overwrite Property converting const char* to string
    static void AddProperty(const std::string &id, const char *value = "", const std::string &description = "");

    /// Set Property Value
    template <class T> static void SetPropertyValue(const std::string &id, const T value)
    {
        if (!ExistsProperty(id))
            throw std::invalid_argument("");

        if (!CheckPropertyType<T>(id))
        {
            std::string propertyType = ConfigurationPropertySupportedTypes::TypeToString(GetPropertyType(id));
            std::string valueType = ConfigurationPropertySupportedTypes::TypeToString<T>();
            throw std::invalid_argument("Property '" + id + "' has type " + propertyType + " and not " + valueType);
        }

        ConfigurationProperty<T> &property = GetProperty<T>(id);
        property.Value = value;
    }

    /// Set Property Value converting const char* to string
    static void SetPropertyValue(const std::string &id, const char *value);

    /// Get Property Value
    template <class T> static const T &GetPropertyValue(const std::string &id)
    {
        if (!ExistsProperty(id))
            throw std::invalid_argument("");

        if (!CheckPropertyType<T>(id))
        {
            std::string propertyType = ConfigurationPropertySupportedTypes::TypeToString(GetPropertyType(id));
            std::string valueType = ConfigurationPropertySupportedTypes::TypeToString<T>();
            throw std::invalid_argument("Property '" + id + "' has type " + propertyType + " and not " + valueType);
        }

        const ConfigurationProperty<T> &property = GetProperty<T>(id);
        return property.Value;
    }
};
} // namespace Gedim

#endif // __CONFIGURATIONS_H
