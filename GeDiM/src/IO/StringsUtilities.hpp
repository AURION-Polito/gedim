#ifndef __StringsUtilities_H
#define __StringsUtilities_H

#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace Gedim
{
class StringsUtilities
{
  public:
    /// Divide a string by a character in a vector of strings
    /// @example stringToSplit="pippo_pe" character='_' -> result=["pippo", "pe"]
    static std::vector<std::string> Split(const std::string &stringToSplit, const char &character = ' ');

    /// Divide a string by a set of characters in a vector of strings
    /// @example stringToSplit="pippo_pe:pu" characters={'_',':'} -> result=["pippo", "pe", "pu"]
    static std::vector<std::string> Split(const std::string &stringToSplit,
                                          const std::vector<char> &characters = std::vector<char>(' '));

    /// Find inside a string a separator between two keys
    /// @example stringToSearch="id:value" keyOne="id" keyTwo="value" -> separator=':'
    static char FindSeparator(const std::string &stringToSearch, const std::string &keyOne, const std::string &keyTwo);

    /// Convert string to lower
    static std::string ToLower(const std::string &input);

    /// Convert string to upper
    static std::string ToUpper(const std::string &input);

    /// Parse a string to object
    template <class T> static T Parse(const std::string &objectString);
};
} // namespace Gedim

#endif // __StringsUtilities_H
