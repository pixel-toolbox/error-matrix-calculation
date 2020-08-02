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
	double oscillationSum_Sq = 0;
	for (unsigned int j=2; j < responseMatrix.columnCount; j++) {
		for (unsigned int i=0; i < j-1; i++) {
			oscillationSum += std::fabs(responseMatrix.m(i, j).value - responseMatrix.m(i+1, j).value);
			oscillationSum_Sq += std::fabs(pow(responseMatrix.m(i, j).value - responseMatrix.m(i+1, j).value, 2));
		}
	}
	for (unsigned int j=0; j < responseMatrix.columnCount-2; j++) {
		for (unsigned int i=j+1; i < responseMatrix.rowCount-1; i++) {
			oscillationSum += std::fabs(responseMatrix.m(i, j).value - responseMatrix.m(i+1, j).value);
			oscillationSum_Sq += std::fabs(pow(responseMatrix.m(i, j).value - responseMatrix.m(i+1, j).value, 2));
		}
	}
	std::cout << oscillationSum << std::endl;
	std::cout << oscillationSum_Sq << std::endl;

	double oscillationSumNaive = 0;
	double oscillationSumNaive_Sq = 0;

	for (unsigned int j=0; j < responseMatrix.columnCount-2; j++) {
		for (unsigned int i=0; i < responseMatrix.rowCount-1; i++) {
			oscillationSumNaive += std::fabs(responseMatrix.m(i, j).value - responseMatrix.m(i+1, j).value);
			oscillationSumNaive_Sq += std::fabs(pow(responseMatrix.m(i, j).value - responseMatrix.m(i+1, j).value, 2));
		}
	}

	for (unsigned int i=0; i < responseMatrix.columnCount-2; i++) {
		for (unsigned int j=0; j < responseMatrix.rowCount-1; j++) {
			oscillationSumNaive += std::fabs(responseMatrix.m(j, i).value - responseMatrix.m(j, i+1).value);
			oscillationSumNaive_Sq += std::fabs(pow(responseMatrix.m(j, i).value - responseMatrix.m(j, i+1).value, 2));
		}
	}

	std::cout << oscillationSumNaive << std::endl;
	std::cout << oscillationSumNaive_Sq << std::endl;

    return 0;
}
