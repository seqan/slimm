#ifndef __STRINGUTILS_H_
#define __STRINGUTILS_H_

#include <string>
#include <vector>

using namespace std;

class StringUtils
{

public:
static void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

static void split(const string& str, vector<string>*out, const string& delim=",")
{
    string::size_type lpos = 0;
    string::size_type pos = str.find_first_of(delim, lpos);
    while(lpos != string::npos)
    {
        out->insert(out->end(), str.substr(lpos,pos - lpos));
        lpos = ( pos == string::npos ) ?  string::npos : pos + 1;
        pos = str.find_first_of(delim, lpos);
    }
}
};
#endif
