#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <stdio.h>

#include "../tasks.hh"
#include "../util/split.hh"
#include "../classes/BinnedTypedMatrix.hh"

namespace fs = std::experimental::filesystem;

using namespace EMC;

int EMC::assemble(bool forceOverwrite, std::string input, std::string output,
		std::vector<double> inputBinCenters, matrixType matType, int countRate) {

	if (!forceOverwrite && fs::exists(output)) {
		std::cout << "Output file already exists, use -f to force overwrite." << std::endl;
		return 1;
	}

	if (matType == mtDensity) {
		std::cout << "Cannot convert to density" << std::endl;
		return 2;
	}

    std::vector<double> outputBinCenters;

	{ // first get output bins from first file
    	char tmpFilename[256];
    	sprintf(tmpFilename, input.c_str(), inputBinCenters[0]);

		std::string filename(tmpFilename);
		std::ifstream outputBinFile(filename);

		while (!outputBinFile.eof()) {
			std::string line;
			std::getline(outputBinFile, line );

			if (!line.empty()) {
				try {
					double binCntr = std::stod(split(split(split(line, ';'), ','), ' ')[0]);
					outputBinCenters.push_back(binCntr);
				} catch(...) { }
			}
		}

		outputBinFile.close();
	}

	BinnedTypedMatrix thisMatrix(outputBinCenters, inputBinCenters, mtNumber);

	int columnCount = 0;
    for (double inBinCenter : inputBinCenters) {
    	char tmpFilename[256];
    	sprintf(tmpFilename, input.c_str(), inBinCenter);

    	std::string filename(tmpFilename);
    	std::ifstream currentSpectrumFile(filename);

		int rowCount = 0;
		while (!currentSpectrumFile.eof()) {
			std::string line;
			std::getline(currentSpectrumFile, line );

			if (!line.empty()) {
				try {
					double value = std::stod(split(split(split(line, ';'), ','), ' ')[1]);
					thisMatrix.m(rowCount, columnCount) = ValueError(value, value);
					rowCount++;
				} catch(...) { }
			}
		}

		currentSpectrumFile.close();

		columnCount++;
    }

    if (matType == mtProbability) {
    	thisMatrix = thisMatrix.makeProbability(countRate);
    }

    std::ofstream outputF(output);
    thisMatrix.writeToFile(outputF);
    outputF.close();

    return 0;
}
