#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <stdio.h>
#include <iomanip>
#include <algorithm>

#include "../tasks.hh"
#include "../util/split.hh"
#include "../classes/BinnedTypedMatrix.hh"

namespace fs = std::experimental::filesystem;

using namespace EMC;

double e = 2.718281828459;

double Pois(double k, double l) {
	return pow(l, k) * exp(-l) / tgamma(k + 1.0);
}

double bisect_lambda(double hitProb) {

	double max = 10;
	double min = 0;

	while ((max-min) > (hitProb/10000)) {
		if ( (1-Pois(0, (max+min)/2.0)) > hitProb) {
			max = (max+min)/2.0;
		} else {
			min = (max+min)/2.0;
		}
	}

	return (max+min)/2.0;
}

std::vector<double> convolute(std::vector<double> input, std::vector<double> inputBinBorders, std::vector<double> outputBinCenters, double order_1, double order_2, double totalHits) {
	std::vector<double> retValue(outputBinCenters.size(), 0);

	for (unsigned int i=0; i<input.size(); i++) {
		retValue.at(i) += order_1 * input.at(i);
	}

	for (unsigned int x=0; x<outputBinCenters.size(); x++) {
		for (unsigned int i=1; i<inputBinBorders.size(); i++) {
			for (unsigned int j=1; j<inputBinBorders.size(); j++) {
				retValue.at(x) += order_2 * input.at(i) * input.at(j) * (1/totalHits) *
						std::max(0.0,
							std::min(inputBinBorders.at(j), outputBinCenters.at(x) - inputBinBorders.at(i-1))
							- std::max(outputBinCenters.at(x) - inputBinBorders.at(i),  inputBinBorders.at(j-1))
						);
			}
		}
	}

	return retValue;
}

int EMC::de_pileup(bool forceOverwrite, std::string input, std::string output, int pixelCount, double lowerBound) {

	if (!forceOverwrite && fs::exists(output)) {
		std::cout << "Output file already exists, use -f to force overwrite."
				<< std::endl;
		return 1;
	}

	int totalCount = 0;

	std::vector<double> inputBinCenter;
	std::vector<double> inputBinValue;
	std::ifstream inputVectorFile(input);
	while (!inputVectorFile.eof()) {
		std::string line;
		std::getline(inputVectorFile, line);

		if (!line.empty()) {
			try {
				auto splitLine = split(split(split(line, ';'), ','), ' ');
				double binCntr = std::stod(splitLine[0]);
				double binValue = std::stod(splitLine[1]);
				inputBinCenter.push_back(binCntr);
				if (binCntr <= lowerBound) {
					inputBinValue.push_back(0);
				} else {
					inputBinValue.push_back(binValue);
					totalCount += binValue;
				}
			} catch (...) {
			}
		}
	}
	inputVectorFile.close();

	double p_hit = (double) totalCount / (double) pixelCount;
	std::cout << "p_hit = " << p_hit << std::endl;
	double current_lambda = bisect_lambda(p_hit);
	std::cout << "   ==> estimated lambda = " << current_lambda << std::endl;

	double order_0_count = Pois(0, current_lambda);
	double order_1_count = Pois(1, current_lambda);
	double order_2_count = Pois(2, current_lambda);
	double order_n_count = 1-order_0_count-order_1_count-order_2_count;

	double order_1_relative = order_1_count / (order_1_count+order_2_count);
	double order_2_relative = order_2_count / (order_1_count+order_2_count);
	double error = order_n_count / (order_1_count+order_2_count);

	std::cout << "calculating orders ..." << std::endl;
	std::cout << "0: " << std::setw(10) << order_0_count << " <--> " << std::setw(10) << (order_0_count*pixelCount) << std::endl;
	std::cout << "1: " << std::setw(10) << order_1_count << " <--> " << std::setw(10) << (order_1_count*pixelCount) << std::endl;
	std::cout << "2: " << std::setw(10) << order_2_count << " <--> " << std::setw(10) << (order_2_count*pixelCount) << std::endl;
	std::cout << "n: " << std::setw(10) << order_n_count << " <--> " << std::setw(10) << (order_n_count*pixelCount) << std::endl;

	std::cout << "relative     1: " << order_1_relative << std::endl;
	std::cout << "relative     2: " << order_2_relative << std::endl;
	std::cout << "relative error: " << error << std::endl;

	std::vector<double> inputBinBorders;
	inputBinBorders.push_back(0);
	for (unsigned int i=1; i<inputBinCenter.size(); i++) {
		inputBinBorders.push_back((inputBinCenter.at(i)+inputBinCenter.at(i-1))/2);
	}

	// first assume deconvoluted = convoluted spectrum!
	std::vector<double> deconvSpectrum;
	for (auto v : inputBinValue) {
		// count is scaled up in order to represent the total count correctly !
		deconvSpectrum.push_back(v * (order_1_count+order_2_count) / order_1_count);
	}

	double meanSquareDiff = 1000;
	{
		double totalCount = (order_1_count+order_2_count)*pixelCount;
		unsigned runid = 0;
		while (meanSquareDiff > 0.001 && runid < 20) {

			runid++;
			std::cout << std::endl << "run " << runid << std::endl;


			std::cout << "convoluting ..." << std::endl;
			std::vector<double> convolutedSpectrum = convolute(deconvSpectrum, inputBinBorders, inputBinCenter, order_1_relative, order_2_relative, totalCount);
			std::cout << "            done" << std::endl;

			std::vector<double> delta;

			for (unsigned int j = 0; j<inputBinValue.size(); j++) {
				delta.push_back(inputBinValue.at(j) - convolutedSpectrum.at(j));
			}

			std::cout << "subtracting first order ..." << std::endl;
			for (unsigned int j = 0; j<inputBinValue.size(); j++) {
				deconvSpectrum.at(j) += delta.at(j) * order_1_relative; // * 0.7;
			}

			/*std::cout << "subtracting second order ..." << std::endl;
			for (unsigned int j = 0; j<inputBinValue.size(); j++) {
				for (unsigned int k = 0; k<j; k++){
					deconvSpectrum.at(k) += delta.at(j) * order_2_relative * 0.7 / ((double) j);
				}
			}*/

			double totalCount = 0;
			for (unsigned int i = 0; i<deconvSpectrum.size(); i++) {
				totalCount += deconvSpectrum.at(i);
			}
			current_lambda = (double) totalCount / (double) pixelCount;
			std::cout << "new lambda: " << current_lambda << std::endl;

			order_0_count = Pois(0, current_lambda);
			order_1_count = Pois(1, current_lambda);
			order_2_count = Pois(2, current_lambda);
			order_n_count = 1-order_0_count-order_1_count-order_2_count;

			order_1_relative = order_1_count / (order_1_count+order_2_count);
			order_2_relative = order_2_count / (order_1_count+order_2_count);
			error = order_n_count / (order_1_count+order_2_count);


			meanSquareDiff = 0;
			for (unsigned int j = 0; j<inputBinValue.size(); j++) {
				meanSquareDiff += pow(convolutedSpectrum.at(j) - inputBinValue.at(j), 2);
			}
			meanSquareDiff /= inputBinValue.size();
			std::cout << "meanSquareDiff: " << meanSquareDiff << std::endl;
		}
	}

    std::ofstream outputFcntr(output);

    outputFcntr << "# relativeError: " << error << std::endl;
    outputFcntr << "# meanSquareDiff: " << meanSquareDiff << std::endl;
    outputFcntr << "# fitted lambda: " << current_lambda << std::endl;

    for (unsigned int i=0; i<inputBinCenter.size(); i++) {
    	outputFcntr << inputBinCenter.at(i) << " " << deconvSpectrum.at(i) << std::endl;
	}

    outputFcntr.close();

	return 0;
}
