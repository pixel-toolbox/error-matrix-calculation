#include <fstream>
#include <iostream>
#include <experimental/filesystem>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
/*#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>*/

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

		std::cout << inputBinCenter[i] << " --> [" << responseMatrix.columnIndex.at(fillPosition) << ", " << responseMatrix.columnIndex.at(fillPosition + 1) << "]" << std::endl;

		if (fillPosition >= 0) {
			inputVector[fillPosition].value += inputBinValue[i];
			inputVectorBinWidth[fillPosition] += 1; //// TODO: change this to measured bin width at correct position
		}
	}

	for (unsigned int i=0; i<responseMatrix.columnCount; i++) {
		inputVector[i].err_sq = inputVector[i].value;
		inputVector[i].value /= inputVectorBinWidth[i];
		inputVector[i].err_sq /= (inputVectorBinWidth[i]*inputVectorBinWidth[i]);
		std::cout << "inputVectorBinWidth[i] = " << inputVectorBinWidth[i] << std::endl;
	}

	for (auto val : inputVector) {
		std::cout << "  (" << val.value << ", " << sqrt(val.err_sq) << ")" << std::endl;
	}

	ublas::vector<ValueError> outputVector = ublas::prod(responseMatrix.m, inputVector);

	for (unsigned int i=0; i<responseMatrix.rowCount; i++) {
		std::cout << responseMatrix.rowIndex[i]   << " " << outputVector[i].value << std::endl;
		std::cout << responseMatrix.rowIndex[i+1] << " " << outputVector[i].value << std::endl;
	}

    return 0;
}


    /*ublas::vector<ValueError> v1(200);

    ublas::matrix<ValueError> m1 = ublas::identity_matrix<ValueError>(200);

    for (int i=0; i<200; i++) {
    	v1[i] = ValueError(i, sqrt(i));
    	m1(i, i) = ValueError(200-i, 42);
    }

    for (ValueError ve : ublas::prod(m1, v1)) {*/
