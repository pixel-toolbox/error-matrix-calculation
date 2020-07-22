#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <stdio.h>

#include "../tasks.hh"
#include "../util/split.hh"
#include "../classes/BinnedTypedMatrix.hh"

namespace fs = std::experimental::filesystem;

using namespace EMC;

int EMC::rebin(bool forceOverwrite, std::string input, std::string output, int count) {

	if (!forceOverwrite && fs::exists(output)) {
		std::cout << "Output file already exists, use -f to force overwrite." << std::endl;
		return 1;
	}

	std::ifstream infile(input);
	BinnedTypedMatrix inputMatrix = BinnedTypedMatrix::readFromFile(infile);



    return 0;
}
