#ifndef EMC_split_hh_included
#define EMC_split_hh_included

#include <vector>
#include <string>
#include <sstream>

namespace EMC {

inline std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
    	if (token != "")
    		tokens.push_back(token);
    }
    return tokens;
}

inline std::vector<std::string> split(const std::vector<std::string>& v, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    for ( auto const& s : v) {
        std::istringstream tokenStream(s);
        while (std::getline(tokenStream, token, delimiter)) {
        	if (token != "")
        		tokens.push_back(token);
        }
    }
    return tokens;
}

} /* namespace EMC */

#endif /* EMC_split_hh_included */
