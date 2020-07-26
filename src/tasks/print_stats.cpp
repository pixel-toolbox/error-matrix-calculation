#include <fstream>
#include <iostream>

#include "../tasks.hh"
#include "../classes/BinnedTypedMatrix.hh"

using namespace EMC;

int EMC::print_stats(std::string matrixFileName) {

	std::ifstream matFile(matrixFileName);
	BinnedTypedMatrix responseMatrix = BinnedTypedMatrix::readFromFile(matFile);
	matFile.close();

	std::cout << "dimension:" << std::endl;
	std::cout << responseMatrix.columnCount << std::endl;

	std::cout << "oscillationSum:" << std::endl;
	double oscillationSum = 0;
	for (unsigned int i=0; i < responseMatrix.rowCount-2; i++) {
		oscillationSum += std::fabs(responseMatrix.m(i, responseMatrix.columnCount-1).value - responseMatrix.m(i+1, responseMatrix.columnCount-1).value);
	}
	std::cout << oscillationSum << std::endl;

    return 0;
}
