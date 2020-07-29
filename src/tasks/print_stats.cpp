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
	for (unsigned int j=2; j < responseMatrix.columnCount; j++) {
		for (unsigned int i=0; i < j-1; i++) {
			oscillationSum += std::fabs(responseMatrix.m(i, j).value - responseMatrix.m(i+1, j).value);
		}
	}
	for (unsigned int j=0; j < responseMatrix.columnCount-2; j++) {
		for (unsigned int i=j+1; i < responseMatrix.rowCount-1; i++) {
			oscillationSum += std::fabs(responseMatrix.m(i, j).value - responseMatrix.m(i+1, j).value);
		}
	}
	std::cout << oscillationSum << std::endl;

    return 0;
}
