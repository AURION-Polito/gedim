#include "Configurations.hpp"
#include "IOUtilities.hpp"
#include "StringsUtilities.hpp"
// #include "CsvExporter.hpp"
#include <set>

using namespace Eigen;

namespace Gedim
{
// ***************************************************************************
std::unordered_map<std::string, ConfigurationPropertySupportedTypes::SupportedTypes> Configurations::propertyTypes;
std::unordered_map<std::string, void *> Configurations::properties;
unsigned int Configurations::numberProperties = 0;
typedef unsigned int uint;
// ***************************************************************************
void Configurations::RemoveProperty(const std::string &id, const ConfigurationPropertySupportedTypes::SupportedTypes &type)
{
    switch (type)
    {
    case ConfigurationPropertySupportedTypes::Int:
        delete static_cast<const ConfigurationProperty<int> *>(properties[id]);
        break;
    case ConfigurationPropertySupportedTypes::Char:
        delete static_cast<const ConfigurationProperty<char> *>(properties[id]);
        break;
    case ConfigurationPropertySupportedTypes::String:
        delete static_cast<const ConfigurationProperty<std::string> *>(properties[id]);
        break;
    case ConfigurationPropertySupportedTypes::Double:
        delete static_cast<const ConfigurationProperty<double> *>(properties[id]);
        break;
    case ConfigurationPropertySupportedTypes::VectorInt:
        delete static_cast<const ConfigurationProperty<std::vector<int>> *>(properties[id]);
        break;
    case ConfigurationPropertySupportedTypes::VectorDouble:
        delete static_cast<const ConfigurationProperty<std::vector<double>> *>(properties[id]);
        break;
    case ConfigurationPropertySupportedTypes::Bool:
        delete static_cast<const ConfigurationProperty<bool> *>(properties[id]);
        break;
    case ConfigurationPropertySupportedTypes::UInt:
        delete static_cast<const ConfigurationProperty<unsigned int> *>(properties[id]);
        break;
    case ConfigurationPropertySupportedTypes::VectorUInt:
        delete static_cast<const ConfigurationProperty<std::vector<unsigned int>> *>(properties[id]);
        break;
    case ConfigurationPropertySupportedTypes::SetUInt:
        delete static_cast<const ConfigurationProperty<std::set<unsigned int>> *>(properties[id]);
        break;
    case ConfigurationPropertySupportedTypes::Matrix:
        delete static_cast<const ConfigurationProperty<MatrixXd> *>(properties[id]);
        break;
    default:
        throw std::invalid_argument("Property '%s' cannot be deleted. Type not recognized");
    }

    properties.erase(id);
    numberProperties--;
}
// ***************************************************************************
Configurations::ExportProperty Configurations::GetPropertyForExport(const std::string &id,
                                                                    const ConfigurationPropertySupportedTypes::SupportedTypes &type)
{
    Configurations::ExportProperty exportProperty;

    std::ostringstream propertyStringConverter;
    std::string propertyDescription = "";

    switch (type)
    {
    case ConfigurationPropertySupportedTypes::Int: {
        propertyStringConverter << GetProperty<int>(id).Value;
        propertyDescription = GetProperty<int>(id).Description;
        break;
    }
    case ConfigurationPropertySupportedTypes::Char: {
        propertyStringConverter << GetProperty<char>(id).Value;
        propertyDescription = GetProperty<char>(id).Description;
        break;
    }
    case ConfigurationPropertySupportedTypes::String: {
        propertyStringConverter << GetProperty<std::string>(id).Value;
        propertyDescription = GetProperty<std::string>(id).Description;
        break;
    }
    case ConfigurationPropertySupportedTypes::Double: {
        propertyStringConverter << GetProperty<double>(id).Value;
        propertyDescription = GetProperty<double>(id).Description;
        break;
    }
    case ConfigurationPropertySupportedTypes::VectorInt: {
        propertyStringConverter << GetProperty<std::vector<int>>(id).Value;
        propertyDescription = GetProperty<std::vector<int>>(id).Description;
        break;
    }
    case ConfigurationPropertySupportedTypes::VectorDouble: {
        propertyStringConverter << GetProperty<std::vector<double>>(id).Value;
        propertyDescription = GetProperty<std::vector<double>>(id).Description;
        break;
    }
    case ConfigurationPropertySupportedTypes::Bool: {
        propertyStringConverter << GetProperty<bool>(id).Value;
        propertyDescription = GetProperty<bool>(id).Description;
        break;
    }
    case ConfigurationPropertySupportedTypes::UInt: {
        propertyStringConverter << GetProperty<unsigned int>(id).Value;
        propertyDescription = GetProperty<unsigned int>(id).Description;
        break;
    }
    case ConfigurationPropertySupportedTypes::VectorUInt: {
        propertyStringConverter << GetProperty<std::vector<unsigned int>>(id).Value;
        propertyDescription = GetProperty<std::vector<unsigned int>>(id).Description;
        break;
    }
    case ConfigurationPropertySupportedTypes::SetUInt: {
        propertyStringConverter << GetProperty<std::set<unsigned int>>(id).Value;
        propertyDescription = GetProperty<std::set<unsigned int>>(id).Description;
        break;
    }
    case ConfigurationPropertySupportedTypes::Matrix: {
        propertyStringConverter << GetProperty<MatrixXd>(id).Value;
        propertyDescription = GetProperty<MatrixXd>(id).Description;
        break;
    }
    default:
        throw std::invalid_argument("Property '" + id + "' cannot be converted. Type not recognized");
    }

    exportProperty.Value = propertyStringConverter.str();
    exportProperty.Description = propertyDescription;

    return exportProperty;
}
// ***************************************************************************
void Configurations::ConvertPropertyFromString(const std::string &id, const std::string &type, const std::string &value, const std::string &description)
{
    auto convertedType = ConfigurationPropertySupportedTypes::StringToType(type);

    switch (convertedType)
    {
    case ConfigurationPropertySupportedTypes::UInt: {
        AddProperty(id, uint(), description);

        unsigned int parsedValue = StringsUtilities::Parse<unsigned int>(value);
        SetPropertyValue(id, parsedValue);

        break;
    }
    case ConfigurationPropertySupportedTypes::Int: {
        AddProperty(id, int(), description);

        int parsedValue = StringsUtilities::Parse<int>(value);
        SetPropertyValue(id, parsedValue);

        break;
    }
    case ConfigurationPropertySupportedTypes::Char: {
        AddProperty(id, char(), description);

        char parsedValue = StringsUtilities::Parse<char>(value);
        SetPropertyValue(id, parsedValue);

        break;
    }
    case ConfigurationPropertySupportedTypes::String: {
        AddProperty(id, std::string(), description);

        std::string parsedValue = StringsUtilities::Parse<std::string>(value);
        SetPropertyValue(id, parsedValue);

        break;
    }
    case ConfigurationPropertySupportedTypes::Double: {
        AddProperty(id, double(), description);

        double parsedValue = StringsUtilities::Parse<double>(value);
        SetPropertyValue(id, parsedValue);

        break;
    }
    case ConfigurationPropertySupportedTypes::VectorUInt: {
        AddProperty(id, std::vector<unsigned int>(), description);

        std::vector<unsigned int> parsedValue = StringsUtilities::Parse<std::vector<unsigned int>>(value);
        SetPropertyValue(id, parsedValue);

        break;
    }
    case ConfigurationPropertySupportedTypes::VectorInt: {
        AddProperty(id, std::vector<int>(), description);

        std::vector<int> parsedValue = StringsUtilities::Parse<std::vector<int>>(value);
        SetPropertyValue(id, parsedValue);

        break;
    }
    case ConfigurationPropertySupportedTypes::VectorDouble: {
        AddProperty(id, std::vector<double>(), description);

        std::vector<double> parsedValue = StringsUtilities::Parse<std::vector<double>>(value);
        SetPropertyValue(id, parsedValue);

        break;
    }
    case ConfigurationPropertySupportedTypes::SetUInt: {
        AddProperty(id, std::set<unsigned int>(), description);

        std::set<unsigned int> parsedValue = StringsUtilities::Parse<std::set<unsigned int>>(value);
        SetPropertyValue(id, parsedValue);

        break;
    }
    case ConfigurationPropertySupportedTypes::Bool: {
        AddProperty(id, bool(), description);

        bool parsedValue = StringsUtilities::Parse<bool>(value);
        SetPropertyValue(id, parsedValue);

        break;
    }
    default:
        throw std::invalid_argument("Property '" + id + "': type not supported");
    }
}
// ***************************************************************************
void Configurations::AddProperty(const std::string &id, const char *value, const std::string &description)
{
    RemoveProperty(id);

    properties.insert(make_pair(id, new ConfigurationProperty<std::string>()));
    propertyTypes.insert(make_pair(id, ConfigurationPropertySupportedTypes::String));
    numberProperties++;

    ConfigurationProperty<std::string> &property = GetProperty<std::string>(id);

    property.Id = id;
    property.Value = std::string(value);
    property.Description = description;
}
// ***************************************************************************
void Configurations::SetPropertyValue(const std::string &id, const char *value)
{
    if (!ExistsProperty(id))
        throw std::invalid_argument("");

    ConfigurationProperty<std::string> &property = GetProperty<std::string>(id);
    property.Value = std::string(value);
}
// ***************************************************************************
void Configurations::InitializeFromCsv(const std::string &inputFile, const bool &hasHeader, const char &separator)
{
    if (!Output::FileExists(inputFile))
        throw std::runtime_error("File '" + inputFile + "' not found");

    std::ifstream inFile;
    inFile.open(inputFile.c_str());

    if (!inFile.is_open())
        throw std::runtime_error("Cannot open file '" + inputFile + "'");

    int num_line = 0;
    while (!inFile.eof())
    {
        std::string line;
        getline(inFile, line);
        num_line++;

        // Skip header
        if (num_line == 1 && hasHeader)
            continue;

        std::vector<std::string> columns = StringsUtilities::Split(line, separator);

        if (columns.size() == 0)
            continue;

        switch (columns.size())
        {
        case 2: {
            const std::string &id = columns[0];
            const std::string &type = "string";
            const std::string &value = columns[1];

            ConvertPropertyFromString(id, type, value);
            break;
        }
        case 3: {
            const std::string &id = columns[0];
            const std::string &type = columns[1];
            const std::string &value = columns[2];

            ConvertPropertyFromString(id, type, value);
            break;
        }
        case 4: {
            const std::string &id = columns[0];
            const std::string &type = columns[1];
            const std::string &value = columns[2];
            const std::string &description = columns[3];

            ConvertPropertyFromString(id, type, value, description);
            break;
        }
        default: {
            Gedim::Output::PrintWarningMessage("Skip property at line " + std::to_string(num_line), true);
            break;
        }
        }
    }

    inFile.close();
}
// ***************************************************************************
void Configurations::InitializeFromIni(const std::string &inputFile)
{
    if (!Output::FileExists(inputFile))
        throw std::runtime_error("File '" + inputFile + "' not found");

    std::ifstream inFile;
    inFile.open(inputFile.c_str());

    if (!inFile.is_open())
        throw std::runtime_error("Cannot open file '" + inputFile + "'");

    // Find Property in File
    unsigned int num_line = 0;
    while (!inFile.eof())
    {
        std::string line;

        getline(inFile, line);
        num_line++;

        // Skip Comment Line
        if (line[0] == '#')
            continue;

        std::vector<std::string> columns = StringsUtilities::Split(line, ' ');

        if (columns.size() == 0)
            continue;

        switch (columns.size())
        {
        case 2: {
            const std::string &id = columns[0];
            const std::string &type = "string";
            const std::string &value = columns[1];

            ConvertPropertyFromString(id, type, value);
            break;
        }
        case 3: {
            const std::string &id = columns[0];
            const std::string &type = columns[1];
            const std::string &value = columns[2];

            ConvertPropertyFromString(id, type, value);
            break;
        }
        default:
            Gedim::Output::PrintWarningMessage("Skip property at line " + std::to_string(num_line), true);
            break;
        }
    }

    inFile.close();
}
// ***************************************************************************
void Configurations::Initialize(const int &argc, char **argv, const std::string &format)
{
    char typeSeparator = StringsUtilities::FindSeparator(format, "id", "type");
    char valueSeparator = StringsUtilities::FindSeparator(format, "type", "value");

    for (int i = 0; i < argc; i++)
    {
        std::vector<std::string> strings = StringsUtilities::Split(argv[i], {typeSeparator, valueSeparator});

        switch (strings.size())
        {
        case 2: {
            const std::string &id = strings[0];
            const std::string &type = "string";
            const std::string &value = strings[1];

            ConvertPropertyFromString(id, type, value);
            break;
        }
        case 3: {
            const std::string &id = strings[0];
            const std::string &type = strings[1];
            const std::string &value = strings[2];

            ConvertPropertyFromString(id, type, value);
            break;
        }
        default:
            break;
        }
    }
}
// ***************************************************************************
void Configurations::ExportToCsv(const std::string &filePath, const bool &append, const char &separator)
{
    std::ofstream file;

    if (append)
        file.open(filePath.c_str(), std::ofstream::app);
    else
        file.open(filePath.c_str());

    if (file.fail())
        throw std::runtime_error("File '" + filePath + "' cannot be opened");

    if (!append)
    {
        file << "Id" << separator;
        file << "Type" << separator;
        file << "Value" << separator;
        file << "Description" << std::endl;
    }

    for (auto it = propertyTypes.begin(); it != propertyTypes.end(); it++)
    {
        const std::string &id = it->first;
        const ConfigurationPropertySupportedTypes::SupportedTypes &type = it->second;
        const std::string propertyType = ConfigurationPropertySupportedTypes::TypeToString(type);
        const ExportProperty property = GetPropertyForExport(id, type);

        const std::string &csvDescription = property.Description;
        const std::string &csvKey = id;
        const std::string &csvType = propertyType;
        const std::string &csvValue = property.Value;

        file << csvKey << separator;
        file << csvType << separator;
        file << csvValue << separator;
        file << csvDescription << std::endl;
    }

    file.close();
}
// ***************************************************************************
void Configurations::ExportToIni(const std::string &filePath, const bool &append, const std::string &section)
{
    std::ofstream file;

    if (append)
        file.open(filePath.c_str(), std::ofstream::app);
    else
        file.open(filePath.c_str());

    if (file.fail())
        throw std::runtime_error("File '" + filePath + "' cannot be opened");

    if (!append && !section.empty())
        file << "# " << section << std::endl;

    for (auto it = propertyTypes.begin(); it != propertyTypes.end(); it++)
    {
        const std::string &id = it->first;
        const ConfigurationPropertySupportedTypes::SupportedTypes &type = it->second;
        const std::string propertyType = ConfigurationPropertySupportedTypes::TypeToString(type);
        const ExportProperty property = GetPropertyForExport(id, type);

        const std::string &iniDescription = property.Description;
        const std::string &iniKey = id;
        const std::string &iniType = propertyType;
        const std::string &iniValue = property.Value;

        if (iniKey.empty() || iniValue.empty() || iniType.empty())
            continue;

        file << "#####################################################" << std::endl;

        if (!iniDescription.empty())
            file << "# " << iniDescription << std::endl;

        file << iniKey << " " << iniType << " " << iniValue << std::endl;
    }

    file.close();
}
// ***************************************************************************
void Configurations::Reset()
{
    for (auto it = propertyTypes.begin(); it != propertyTypes.end(); it++)
    {
        const std::string &id = it->first;
        const ConfigurationPropertySupportedTypes::SupportedTypes &type = it->second;

        RemoveProperty(id, type);
    }

    propertyTypes.clear();
    numberProperties = 0;
}
// ***************************************************************************
void Configurations::RemoveProperty(const std::string &id)
{
    if (!ExistsProperty(id))
        return;

    RemoveProperty(id, GetPropertyType(id));
    propertyTypes.erase(id);
}
// ***************************************************************************
} // namespace Gedim
