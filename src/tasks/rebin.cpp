#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <stdio.h>

#include "../tasks.hh"
#include "../util/split.hh"
#include "../classes/BinnedTypedMatrix.hh"

namespace fs = std::experimental::filesystem;

using namespace EMC;

struct Mapping {

	const double min;
	const double max;
	const double minBinSize;
	const double a;

	Mapping(int count, double min, double minBinSize, double max) :
		min(min), max(max), minBinSize(minBinSize),
		a(2.0*((max-min)-minBinSize*count)/(count*count)) {}

	double valueToBinPosition(double value) {
		if (value <= min) {
			return 0;
		}
		return (-minBinSize+sqrt(minBinSize*minBinSize+2.0*a*(value-min)) ) / a;
	}

	double binPositionToValue(double binNumber) {
		if (binNumber == 0) {
			return 0;
		}
		return minBinSize*binNumber+a*binNumber*binNumber/2.0 + min;
	}
};

int EMC::rebin(bool forceOverwrite, std::string input, std::string output,
		int count, double minRow, double minCol, double minRowBin,
		double minColBin, double maxRow, double maxColumn) {

	if (!forceOverwrite && fs::exists(output)) {
		std::cout << "Output file already exists, use -f to force overwrite."
				<< std::endl;
		return 1;
	}

	std::ifstream infile(input);
	BinnedTypedMatrix inputMatrix = BinnedTypedMatrix::readFromFile(infile);
	infile.close();

	if (inputMatrix.matType != mtProbability) {
		std::cout << "Can only rebin probability Matrix" << std::endl;
		return 2;
	}

	Mapping rowMappingHelper(count, minRow, minRowBin, maxRow);
	Mapping colMappingHelper(count, minCol, minColBin, maxColumn);

	std::vector<double> newRowBinBorders;
	std::vector<double> newColumnBinBorders;
	for (int i=0; i<=count; i++) {
		newRowBinBorders.push_back(rowMappingHelper.binPositionToValue(i));
		newColumnBinBorders.push_back(colMappingHelper.binPositionToValue(i));
	}

	std::cout << "newRowBinBorders: ";
	for (double val : newRowBinBorders) {
		std::cout << val << " ";
	}
	std::cout << std::endl;

	std::cout << "newColumnBinBorders: ";
	for (double val : newColumnBinBorders) {
		std::cout << val << " ";
	}
	std::cout << std::endl;

	std::vector<double> rowBinWidths(newRowBinBorders.size()-1, 0);
	std::vector<int> rowMapping;
	{
		double lower = 0;
		double higher = 0;
		for (unsigned int i=0; inputMatrix.rowIndex[i] < maxRow; i++) {
			lower = higher;
			higher = inputMatrix.rowIndex.at(i) + 0.5 * (inputMatrix.rowIndex.at(i+1) - inputMatrix.rowIndex.at(i));
			unsigned int position = rowMappingHelper.valueToBinPosition(inputMatrix.rowIndex.at(i));
			rowBinWidths.at(position) += higher - lower;
			rowMapping.push_back(position);
			//std::cout << "rebinning row: " << inputMatrix.rowIndex.at(i) << " goes into [" << newRowBinBorders.at(position) << "," << newRowBinBorders.at(position+1) << "]" << std::endl;
		}
		for (unsigned int i=0; i<rowBinWidths.size();i++) {
			std::cout << "row " << i << " has width " << rowBinWidths.at(i) << std::endl;
		}
	}
	std::vector<double> columnBinWidths(newColumnBinBorders.size()-1, 0);
	std::vector<int> columnMapping;
	{
		double lower = 0;
		double higher = 0;
		for (unsigned int i=0; inputMatrix.columnIndex[i] < maxColumn; i++) {
			lower = higher;
			higher = inputMatrix.columnIndex.at(i) + 0.5 * (inputMatrix.columnIndex.at(i+1) - inputMatrix.columnIndex.at(i));
			unsigned int position = colMappingHelper.valueToBinPosition(inputMatrix.columnIndex.at(i));
			columnBinWidths.at(position) += higher - lower;
			columnMapping.push_back(position);
			//std::cout << "rebinning column: " << inputMatrix.columnIndex.at(i) << " goes into [" << newColumnBinBorders.at(position) << "," << newColumnBinBorders.at(position+1) << "]" << std::endl;
		}
		for (unsigned int i=0; i<columnBinWidths.size();i++) {
			std::cout << "column " << i << " has width " << columnBinWidths.at(i) << std::endl;
		}
	}

	BinnedTypedMatrix newMat(newRowBinBorders, newColumnBinBorders, mtDensity);
	for (unsigned int oldRow = 0; inputMatrix.rowIndex.at(oldRow) < maxRow; oldRow++) {
		for (unsigned int oldColumn = 0; inputMatrix.columnIndex.at(oldColumn) < maxColumn; oldColumn++) {
			newMat.m(rowMapping.at(oldRow), columnMapping.at(oldColumn)) +=
					inputMatrix.m(oldRow, oldColumn)
							/ (rowBinWidths.at(rowMapping.at(oldRow))
									* columnBinWidths.at(columnMapping.at(oldColumn)));
		}
	}

	std::ofstream outFile(output);
	newMat.writeToFile(outFile);
	outFile.close();

    return 0;
}
