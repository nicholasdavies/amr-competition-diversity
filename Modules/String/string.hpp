// string.hpp - string utility functions
// (C) 2013 Nicholas G Davies

#ifndef STRING_HPP
#define STRING_HPP

#include <string>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <vector>
#include <sstream>
#include <map>

// Trim left whitespace
inline std::string ltrim(std::string s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// Trim right whitespace
inline std::string rtrim(std::string s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// Trim left and right whitespace
inline std::string trim(std::string s)
{
    return ltrim(rtrim(s));
}

// String tolower
inline std::string lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

// String toupper
inline std::string upper(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

// Serialize a vector of strings with the given separator
inline std::string serialize(std::vector<std::string>& v, std::string sep = " ")
{
    std::string ret;
    for (std::vector<std::string>::iterator i = v.begin(); i != v.end(); ++i)
    {
        ret += *i;
        if (i + 1 != v.end())
            ret += sep;
    }
    return ret;
}

// Unserialize a string into a vector of strings based on the given separator
inline std::vector<std::string> unserialize(std::string s, std::string sep = " ")
{
    std::vector<std::string> ret;
    std::string::size_type n0 = 0, n1;
    while ((n1 = s.find(sep, n0)) != std::string::npos)
    {
        ret.push_back(s.substr(n0, n1 - n0));
        n0 = n1 + sep.length();
    }
    if (n0 != s.length())
        ret.push_back(s.substr(n0));
    return ret;
}

// Unserialize a string into a vector of strings based on any number of intervening whitespace characters being permitted
inline std::vector<std::string> unserialize_ws(std::string s)
{
    std::string sep = " \f\n\r\t\v";
    std::vector<std::string> ret;
    std::string::size_type n0 = s.find_first_not_of(sep), n1;
    while ((n1 = s.find_first_of(sep, n0)) != std::string::npos)
    {
        ret.push_back(s.substr(n0, n1 - n0));
        n0 = s.find_first_not_of(sep, n1);
    }
    if (n0 != std::string::npos)
        ret.push_back(s.substr(n0));
    return ret;
}

// Replace a given token in a string with another.
inline std::string replace(std::string str, std::string from, std::string to)
{
    std::string::size_type p = 0;
    while ((p = str.find(from, p)) != std::string::npos)
    {
        str.replace(p, from.length(), to);
        p += to.length();
    }
    return str;
}

// Make a string safe for storing in a datafile.
inline std::string encode(std::string orig)
{
    return replace(replace(orig, ":", "{[COLON]}"), "\n", "{[NEWLINE]}");
}

// Decode a safeified string.
inline std::string decode(std::string coded)
{
    return replace(replace(coded, "{[COLON]}", ":"), "{[NEWLINE]}", "\n");
}

// Encode and serialize for datafiles
inline std::string encode_and_serialize(std::vector<std::string> v)
{
    for (auto& s : v)
    {
        s = encode(s);
    }
    return serialize(v, ":");
}

// Unserialize and decode for datafiles
inline std::vector<std::string> unserialize_and_decode(std::string s)
{
    std::vector<std::string> ret = unserialize(s, ":");
    for (auto& s : ret)
    {
        s = decode(s);
    }
    return ret;
}

// turn a string of the form "1, 3-5, 10" or "*" into a vector of the form {1, 3, 4, 5, 10}
inline std::vector<int> extract_range(std::string s, int def_min = 0, int def_max = -1)
{
    std::vector<int> range;
    if (def_max >= def_min && s == "*")
    {
        while (def_min <= def_max) range.push_back(def_min++);
    }
    else
    {
        auto elements = unserialize(s, ",");
        for (auto& e : elements)
        {
            std::string::size_type dash;
            if ((dash = e.find('-')) != std::string::npos)
            {
                int first = stoi(e.substr(0, dash)), second = stoi(e.substr(dash + 1));
                while (first <= second) range.push_back(first++);
            }
            else
            {
                range.push_back(stoi(e));
            }
        }
    }
    return range;
}

typedef std::map<std::string, std::string> param_map;

// turn a string of the form "x = y, w = z" into a map of form m["x"] = "y", m["w"] = "z"
inline param_map extract_parameters(std::string s, std::string sep = ",", std::string eq = "=")
{
    param_map m;
    auto bits = unserialize(s, sep);

    for (auto& b : bits)
    {
        auto equals = b.find(eq);
        if (equals != std::string::npos)
            m[trim(b.substr(0, equals))] = trim(b.substr(equals + 1));
    }

    return m;
}

inline bool param_defined(param_map& m, std::string key)
{
    return m.count(key) > 0;
}

inline bool param_defined(param_map& m, std::string key, param_map& s)
{
    return m.count(key) > 0 || s.count(key) > 0;
}

inline std::string param_entry(param_map& m, std::string key, std::string def)
{
    if (m.count(key) > 0) return m[key];
    return def;
}

inline double param_entry(param_map& m, std::string key, double def)
{
    if (m.count(key) > 0) return stod(m[key]);
    return def;
}

inline std::string param_entry(param_map& m, std::string key, std::string def, param_map& s)
{
    if (m.count(key) > 0) return m[key];
    if (s.count(key) > 0) return s[key];
    return def;
}

inline double param_entry(param_map& m, std::string key, double def, param_map& s)
{
    if (m.count(key) > 0) return stod(m[key]);
    if (s.count(key) > 0) return stod(s[key]);
    return def;
}

#endif