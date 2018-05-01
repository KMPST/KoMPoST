// Read an INI file into easy-to-access name/value pairs.

#include <algorithm>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "ini.h"
#include "INIReader.h"

#define FOUND_KEY(ok, value, default_value)   if (ok) { return value; } else { std::cout << "*** Key not found -- " <<section<<"."<<name<<" -- returning default value: "<< default_value << " ***" << std::endl ; return default_value ; }

using std::string;

INIReader::INIReader(string filename)
{
    _error = ini_parse(filename.c_str(), ValueHandler, this);
}

int INIReader::ParseError()
{
    return _error;
}

string INIReader::Get(string section, string name, string default_value)
{
    string key = MakeKey(section, name);
    return (_values.count(key) ? _values[key] : default_value) ;
}
string INIReader::GetString(string section, string name, string default_value)
{
    string key = MakeKey(section, name);
    FOUND_KEY((_values.count(key)),(_values[key]),(default_value))
}

long INIReader::GetInteger(string section, string name, long default_value)
{
    string valstr = Get(section, name, "");
    const char* value = valstr.c_str();
    char* end;
    // This parses "1234" (decimal) and also "0x4D2" (hex)
    long n = strtol(value, &end, 0);
    FOUND_KEY( (end>value), n, default_value) 
}

double INIReader::GetReal(string section, string name, double default_value)
{
    string valstr = Get(section, name, "");
    const char* value = valstr.c_str();
    char* end;
    double n = strtod(value, &end);
    FOUND_KEY( (end>value), n, default_value) 
}

bool INIReader::GetBoolean(string section, string name, bool default_value)
{
    string valstr = Get(section, name, "");
    // Convert to lower case to make string comparisons case-insensitive
    std::transform(valstr.begin(), valstr.end(), valstr.begin(), ::tolower);
    if (valstr == "true" || valstr == "yes" || valstr == "on" || valstr == "1")
        return true;
    else if (valstr == "false" || valstr == "no" || valstr == "off" || valstr == "0")
        return false;
    else
        std::cout << "*** Key not found in " <<section<<"."<<name<<", returning default value: "<< default_value << " ***" << std::endl ; 
    return default_value;
}

string INIReader::MakeKey(string section, string name)
{
    string key = section + "." + name;
    // Convert to lower case to make section/name lookups case-insensitive
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);
    return key;
}

int INIReader::ValueHandler(void* user, const char* section, const char* name,
                            const char* value)
{
    INIReader* reader = (INIReader*)user;
    string key = MakeKey(section, name);
    if (reader->_values[key].size() > 0)
        reader->_values[key] += "\n";
    reader->_values[key] += value;
    return 1;
}
