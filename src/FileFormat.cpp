#include "FileFormat.h"
#include <fstream>

namespace plz2branch {
	
//Format is plz,inhabitants. header is plz,inhabitants
std::unordered_map<uint32_t, uint32_t> plz2inhabitants(std::string const & fileName) {
	std::unordered_map<uint32_t, uint32_t> result;
	std::ifstream file(fileName);
	if (!file.is_open()) {
		throw std::runtime_error("Could not open plz2inhabitants file");
	}
	std::string header;
	std::getline(file, header);
	
	while(!file.eof() && file.good()) {
		uint32_t plz, inhabitants;
		char sep;
		file >> plz >> sep >> inhabitants;
		result[plz] = inhabitants;
	}
	return result;
}
	
} //end namespaceplz
