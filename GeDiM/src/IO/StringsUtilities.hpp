#ifndef __StringsUtilities_H
#define __StringsUtilities_H

#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace Gedim
{
class StringsUtilities
{
  public:
    /// Divide a string by a character in a vector of strings
    /// @example stringToSplit="pippo_pe" character='_' -> result=["pippo", "pe"]
    static vector<string> Split(const string &stringToSplit, const char &character = ' ');

    /// Divide a string by a set of characters in a vector of strings
    /// @example stringToSplit="pippo_pe:pu" characters={'_',':'} -> result=["pippo", "pe", "pu"]
    static vector<string> Split(const string &stringToSplit, const vector<char> &characters = vector<char>(' '));

    /// Find inside a string a separator between two keys
    /// @example stringToSearch="id:value" keyOne="id" keyTwo="value" -> separator=':'
    static char FindSeparator(const string &stringToSearch, const string &keyOne, const string &keyTwo);

    /// Convert string to lower
    static string ToLower(const string &input);

    /// Convert string to upper
    static string ToUpper(const string &input);

    /// Parse a string to object
    template <class T> static T Parse(const string &objectString);
};
} // namespace Gedim

#endif // __StringsUtilities_H
