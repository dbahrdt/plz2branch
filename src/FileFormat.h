#pragma once
#include <stdint.h>
#include <string>
#include <unordered_map>

namespace plz2branch {
	
std::unordered_map<uint32_t, uint32_t> plz2inhabitants(std::string const & fileName);


	
} //end namespace plz2branch
