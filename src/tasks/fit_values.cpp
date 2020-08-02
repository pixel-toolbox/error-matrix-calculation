#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <stdio.h>
#include <iomanip>

#include "../tasks.hh"
#include "../util/split.hh"
#include "../classes/BinnedTypedMatrix.hh"

namespace fs = std::experimental::filesystem;

using namespace EMC;

const std::string rebinned_ending = ".rebin.mat";
const std::string inverted_ending = ".rebin.inverted.mat";

double getOscillationSumForParameter(std::string input, std::string output,
		int count, double minRow, double minCol, double minRowBin,
		double minColBin, double maxRow, double maxCol) {

	rebin(true, input, output + ".rebin.mat", count, minRow, minCol, minRowBin,
		minColBin, maxRow, maxCol);

	invert(true, output + rebinned_ending, output + inverted_ending);

	std::ifstream matFile(output + inverted_ending);
	BinnedTypedMatrix responseMatrix = BinnedTypedMatrix::readFromFile(matFile);
	matFile.close();

	double oscillationSumNaive = 0;

	for (unsigned int j=0; j < responseMatrix.columnCount-2; j++) {
		for (unsigned int i=0; i < responseMatrix.rowCount-1; i++) {
			oscillationSumNaive += std::fabs(responseMatrix.m(i, j).value - responseMatrix.m(i+1, j).value);
		}
	}

	for (unsigned int i=0; i < responseMatrix.columnCount-2; i++) {
		for (unsigned int j=0; j < responseMatrix.rowCount-1; j++) {
			oscillationSumNaive += std::fabs(responseMatrix.m(j, i).value - responseMatrix.m(j, i+1).value);
		}
	}

	return oscillationSumNaive;
}


int EMC::fit_values(std::string input, std::string output,
		int count, double minRow, double minCol, double minRowBin,
		double minColBin, double maxRow, double maxCol) {

	std::cout << "Column count: " << std::endl;
	std::cout << "   minRow,    minCol, minRowBin, minColBin,    maxRow,    maxCol, err_sum" << std::endl;

    std::streambuf* orig_buf = std::cout.rdbuf();
	bool finished = false;
	while (!finished) {
		finished = true;

	    std::cout.rdbuf(orig_buf);
		std::cout << std::setw(9) << minRow    << ", ";
		std::cout << std::setw(9) << minCol    << ", ";
		std::cout << std::setw(9) << minRowBin << ", ";
		std::cout << std::setw(9) << minColBin << ", ";
		std::cout << std::setw(9) << maxRow    << ", ";
		std::cout << std::setw(9) << maxCol    << ", ";
		std::cout.flush();
	    std::cout.rdbuf(NULL);
	    double err_sum = getOscillationSumForParameter(input, output, count, minRow, minCol, minRowBin, minColBin, maxRow, maxCol);
	    std::cout.rdbuf(orig_buf);
		std::cout << std::setw(11) << err_sum;
		std::cout << std::endl;
	    std::cout.rdbuf(NULL);

	    if        (getOscillationSumForParameter(input, output, count, minRow+.5, minCol  , minRowBin  , minColBin  , maxRow, maxCol  ) < err_sum) {
	    	finished = false;
	    	minRow+=.5;
	    } else if (getOscillationSumForParameter(input, output, count, minRow-.5, minCol  , minRowBin  , minColBin  , maxRow, maxCol  ) < err_sum) {
	    	finished = false;
	    	minRow-=.5;
	    } else if (getOscillationSumForParameter(input, output, count, minRow  , minCol+.25, minRowBin  , minColBin  , maxRow, maxCol  ) < err_sum) {
	    	finished = false;
	    	minCol += .25;
	    } else if (getOscillationSumForParameter(input, output, count, minRow  , minCol-.25, minRowBin  , minColBin  , maxRow, maxCol  ) < err_sum) {
	    	finished = false;
	    	minCol -= .25;
	    } else if (getOscillationSumForParameter(input, output, count, minRow  , minCol  , minRowBin  , minColBin  , maxRow+1, maxCol) < err_sum) {
	    	finished = false;
	    	maxRow++;
	    } else if (getOscillationSumForParameter(input, output, count, minRow  , minCol  , minRowBin  , minColBin  , maxRow-1, maxCol) < err_sum) {
	    	finished = false;
	    	maxRow--;
	    } else if (getOscillationSumForParameter(input, output, count, minRow  , minCol  , minRowBin+.1, minColBin  , maxRow, maxCol  ) < err_sum) {
	    	finished = false;
	    	minRowBin+=.1;
	    } else if (getOscillationSumForParameter(input, output, count, minRow  , minCol  , minRowBin-.1, minColBin  , maxRow, maxCol  ) < err_sum) {
	    	finished = false;
	    	minRowBin-=.1;
	    } else if (getOscillationSumForParameter(input, output, count, minRow  , minCol  , minRowBin  , minColBin+.1, maxRow, maxCol  ) < err_sum) {
	    	finished = false;
	    	minColBin+=.1;
	    } else if (getOscillationSumForParameter(input, output, count, minRow  , minCol  , minRowBin  , minColBin-.1, maxRow, maxCol  ) < err_sum) {
	    	finished = false;
	    	minColBin-=.1;
	    }
	}



    return 0;
}
