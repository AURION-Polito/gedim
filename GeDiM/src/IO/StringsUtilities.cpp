#include <algorithm>
#include <cctype>
#include <list>
#include <sstream>

#include "StringsUtilities.hpp"

using namespace std;

namespace Gedim
{
// ***************************************************************************
vector<string> StringsUtilities::Split(const string &stringToSplit, const char &character)
{
    stringstream splitter(stringToSplit);
    string tempString;
    vector<string> strings;

    while (getline(splitter, tempString, character))
        strings.push_back(tempString);

    return strings;
}
// ***************************************************************************
vector<string> StringsUtilities::Split(const string &stringToSplit, const vector<char> &characters)
{
    list<string> stringsTotal;
    stringsTotal.push_back(stringToSplit);

    for (unsigned int c = 0; c < characters.size(); c++)
    {
        const char &character = characters[c];

        list<string> subStrings;

        for (auto it = stringsTotal.begin(); it != stringsTotal.end(); ++it)
        {
            const string &subStringToSplit = *it;
            vector<string> strings = Split(subStringToSplit, character);

            for (unsigned int s = 0; s < strings.size(); s++)
                subStrings.push_back(strings[s]);
        }

        stringsTotal.clear();
        for (auto it = subStrings.begin(); it != subStrings.end(); ++it)
            stringsTotal.push_back(*it);
    }

    return vector<string>(make_move_iterator(stringsTotal.begin()), make_move_iterator(stringsTotal.end()));
}
// ***************************************************************************
char StringsUtilities::FindSeparator(const string &stringToSearch, const string &keyOne, const string &keyTwo)
{
    char separator;

    size_t positionKeyOne = stringToSearch.find(keyOne);
    if (positionKeyOne == std::string::npos)
        throw runtime_error("");

    size_t positionKeyTwo = stringToSearch.find(keyTwo);
    if (positionKeyTwo == std::string::npos)
        throw runtime_error("");

    if (positionKeyTwo == positionKeyOne)
        throw runtime_error("");

    if (positionKeyTwo > positionKeyOne)
    {
        size_t separatorPosition = keyOne.length() + positionKeyOne;

        if (positionKeyTwo - separatorPosition != 1)
            throw runtime_error("");

        separator = stringToSearch[separatorPosition];
    }
    else
    {
        size_t separatorPosition = keyTwo.length() + positionKeyTwo;

        if (positionKeyOne - separatorPosition != 1)
            throw runtime_error("");

        separator = stringToSearch[separatorPosition];
    }

    return separator;
}
// ***************************************************************************
string StringsUtilities::ToLower(const string &input)
{
    string output = input;
    transform(output.begin(), output.end(), output.begin(), [](unsigned char c) { return std::tolower(c); });

    return output;
}
// ***************************************************************************
string StringsUtilities::ToUpper(const string &input)
{
    string output = input;
    transform(output.begin(), output.end(), output.begin(), [](unsigned char c) { return std::toupper(c); });

    return output;
}
// ***************************************************************************
template <> char StringsUtilities::Parse<char>(const string &objectString)
{
    return objectString[0];
}
// ***************************************************************************
template <> int StringsUtilities::Parse<int>(const string &objectString)
{
    return std::stoi(objectString);
}
// ***************************************************************************
template <> unsigned int StringsUtilities::Parse<unsigned int>(const string &objectString)
{
    return std::stoi(objectString);
}
// ***************************************************************************
template <> float StringsUtilities::Parse<float>(const string &objectString)
{
    return std::stof(objectString);
}
// ***************************************************************************
template <> double StringsUtilities::Parse<double>(const string &objectString)
{
    return std::stod(objectString);
}
// ***************************************************************************
template <> string StringsUtilities::Parse<string>(const string &objectString)
{
    return objectString;
}
// ***************************************************************************
template <> bool StringsUtilities::Parse<bool>(const string &objectString)
{
    bool objectConverted;
    const string &objStringToLow = ToLower(objectString);

    if (objStringToLow == "false" || objStringToLow == "f" || objStringToLow == "0")
    {
        objectConverted = false;
        return objectConverted;
    }

    if (objStringToLow == "true" || objStringToLow == "t" || objStringToLow == "1")
    {
        objectConverted = true;
        return objectConverted;
    }

    throw invalid_argument("");
}
// ***************************************************************************
/// @note Supported vector separator are {',' ';' ' '}
template <> vector<int> StringsUtilities::Parse<vector<int>>(const string &objectString)
{
    vector<int> objectConverted;

    vector<string> numberStrings = Split(objectString, {'[', ',', ';', ']', ' '});

    list<int> numbers;
    for (unsigned int s = 0; s < numberStrings.size(); s++)
        numbers.push_back(Parse<int>(numberStrings[s]));

    objectConverted.clear();
    objectConverted.reserve(numbers.size());

    for (auto it = numbers.begin(); it != numbers.end(); it++)
        objectConverted.push_back(*it);

    return objectConverted;
}
// ***************************************************************************
/// @note Supported vector separator are {',' ';' ' '}
template <> vector<unsigned int> StringsUtilities::Parse<vector<unsigned int>>(const string &objectString)
{
    vector<unsigned int> objectConverted;

    vector<string> numberStrings = Split(objectString, {'[', ',', ';', ']', ' '});

    list<unsigned int> numbers;
    for (unsigned int s = 0; s < numberStrings.size(); s++)
        numbers.push_back(Parse<unsigned int>(numberStrings[s]));

    objectConverted.clear();
    objectConverted.reserve(numbers.size());

    for (auto it = numbers.begin(); it != numbers.end(); it++)
        objectConverted.push_back(*it);

    return objectConverted;
}
// ***************************************************************************
/// @note Supported vector separator are {',' ';' ' '}
template <> vector<double> StringsUtilities::Parse<vector<double>>(const string &objectString)
{
    vector<double> objectConverted;

    vector<string> numberStrings = Split(objectString, {'[', ',', ';', ']', ' '});

    list<double> numbers;
    for (unsigned int s = 0; s < numberStrings.size(); s++)
        numbers.push_back(Parse<double>(numberStrings[s]));

    objectConverted.clear();
    objectConverted.reserve(numbers.size());

    for (auto it = numbers.begin(); it != numbers.end(); it++)
        objectConverted.push_back(*it);

    return objectConverted;
}
// ***************************************************************************
/// @note Supported vector separator are {',' ';' ' '}
template <> set<unsigned int> StringsUtilities::Parse<set<unsigned int>>(const string &objectString)
{
    set<unsigned int> objectConverted;

    vector<string> numberStrings = Split(objectString, {'[', ',', ';', ']', ' '});

    list<unsigned int> numbers;
    for (unsigned int s = 0; s < numberStrings.size(); s++)
        numbers.push_back(Parse<unsigned int>(numberStrings[s]));

    objectConverted.clear();

    for (auto it = numbers.begin(); it != numbers.end(); it++)
        objectConverted.insert(*it);

    return objectConverted;
}
// ***************************************************************************
} // namespace Gedim
