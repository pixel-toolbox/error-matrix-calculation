#include <fstream>
#include <iostream>
#include <experimental/filesystem>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "../tasks.hh"
#include "../util/split.hh"
#include "../classes/BinnedTypedMatrix.hh"

namespace fs = std::experimental::filesystem;

using namespace EMC;

using namespace boost::numeric;

int EMC::apply(bool forceOverwrite, std::string inputVectorFileName, std::string outputVectorFileName, std::string matrixFileName) {

	if (!forceOverwrite && fs::exists(outputVectorFileName)) {
		std::cout << "Output file already exists, use -f to force overwrite." << std::endl;
		return 1;
	}

	std::ifstream matFile(matrixFileName);
	BinnedTypedMatrix responseMatrix = BinnedTypedMatrix::readFromFile(matFile);
	matFile.close();

	std::vector<double> inputBinCenter;
	std::vector<double> inputBinValue;
	std::ifstream inputVectorFile(inputVectorFileName);
	while (!inputVectorFile.eof()) {
		std::string line;
		std::getline(inputVectorFile, line );

		if (!line.empty()) {
			try {
				auto splitLine = split(split(split(line, ';'), ','), ' ');
				double binCntr = std::stod(splitLine[0]);
				double binValue = std::stod(splitLine[1]);
				inputBinCenter.push_back(binCntr);
				inputBinValue.push_back(binValue);
			} catch(...) { }
		}
	}
	inputVectorFile.close();

	ublas::vector<ValueError> inputVector(responseMatrix.columnCount);
	std::vector<double> inputVectorBinWidth(responseMatrix.columnCount, 0);

	for (unsigned int i = 0; i < inputBinCenter.size() && inputBinCenter[i] < responseMatrix.columnIndex.back() ; i++) {
		int fillPosition = -1;
		while ( (inputBinCenter[i] > responseMatrix.columnIndex.at(fillPosition + 1)) ) {
			fillPosition++;
		}

		//std::cout << inputBinCenter[i] << " --> [" << responseMatrix.columnIndex.at(fillPosition) << ", " << responseMatrix.columnIndex.at(fillPosition + 1) << "]" << std::endl;

		if (fillPosition >= 0) {
			inputVector[fillPosition].value += inputBinValue[i];
			if (i == 0) {
				inputVectorBinWidth[fillPosition] += inputBinCenter[i+1] - inputBinCenter[i];
			} else {
				inputVectorBinWidth[fillPosition] += inputBinCenter[i] - inputBinCenter[i-1];
			}
		}
	}

	for (unsigned int i=0; i<responseMatrix.columnCount; i++) {
		inputVector[i].err_sq = inputVector[i].value;
		inputVector[i].value /= inputVectorBinWidth[i];
		inputVector[i].err_sq /= (inputVectorBinWidth[i]*inputVectorBinWidth[i]);
		std::cout << "inputVectorBinWidth[i] = " << inputVectorBinWidth[i] << std::endl;
	}

	ublas::vector<ValueError> outputVector = ublas::prod(responseMatrix.m, inputVector);

    std::ofstream outputFstep(outputVectorFileName + ".step");

    for (unsigned int i=0; i<responseMatrix.rowCount; i++) {
    	outputFstep << responseMatrix.rowIndex[i]   << " " << outputVector[i].value << std::endl;
    	outputFstep << responseMatrix.rowIndex[i+1] << " " << outputVector[i].value << std::endl;
	}

    outputFstep.close();

    std::ofstream outputFcntr(outputVectorFileName + ".cntr");

    for (unsigned int i=0; i<responseMatrix.rowCount; i++) {
    	outputFcntr << (0.5*(responseMatrix.rowIndex[i] + responseMatrix.rowIndex[i+1]))   << " " << outputVector[i].value << " " << sqrt(outputVector[i].err_sq) << std::endl;
	}

    outputFcntr.close();

    return 0;
}
