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
	infile.close();

	if (inputMatrix.matType != mtProbability) {
		std::cout << "Can only rebin probability Matrix" << std::endl;
		return 2;
	}

	unsigned int rowCount = inputMatrix.rowIndex.size();
	double maxRow = inputMatrix.rowIndex.at(rowCount-1) +  0.5 * (inputMatrix.rowIndex.at(rowCount-1) - inputMatrix.rowIndex.at(rowCount-2) );
	double aRow = maxRow / log(count);
	std::cout << "rowCount = " << rowCount << ", maxRow = " << maxRow << ", aRow = " << aRow << std::endl;

	unsigned int columnCount = inputMatrix.columnIndex.size();
	double maxColumn = inputMatrix.columnIndex.at(columnCount-1) +  0.5 * (inputMatrix.columnIndex.at(columnCount-1) - inputMatrix.columnIndex.at(columnCount-2) );
	double aColumn = maxColumn / log(count);
	std::cout << "columnCount = " << columnCount << ", maxColumn = " << maxColumn << ", aColumn = " << aColumn << std::endl;

	std::vector<double> newRowBinBorders;
	std::vector<double> newColumnBinBorders;
	for (int i=1; i<=count; i++) {
		newRowBinBorders.push_back(aRow*log(i));
		newColumnBinBorders.push_back(aColumn*log(i));
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

	std::vector<double> rowBinWidths(newRowBinBorders.size(), 0);
	std::vector<int> rowMapping;
	{
		double lower = 0;
		double higher = 0;
		for (unsigned int i=0; i<(inputMatrix.rowIndex.size()-1); i++) {
			lower = higher;
			higher = inputMatrix.rowIndex.at(i) + 0.5 * (inputMatrix.rowIndex.at(i+1) - inputMatrix.rowIndex.at(i));
			unsigned int position = exp(inputMatrix.rowIndex.at(i)/aRow)-1;
			rowBinWidths.at(position) += higher - lower;
			rowMapping.push_back(position);
			std::cout << "rebinning row: " << inputMatrix.rowIndex.at(i) << " goes into [" << newRowBinBorders.at(position) << "," << newRowBinBorders.at(position+1) << "]" << std::endl;
		}
		lower = higher;
		higher = inputMatrix.rowIndex.back() + 0.5 * (inputMatrix.rowIndex.back() - inputMatrix.rowIndex.at(inputMatrix.rowIndex.size()-2));
		unsigned int position = exp(inputMatrix.rowIndex.back()/aRow)-1;
		rowBinWidths.at(position) += higher - lower;
		rowMapping.push_back(position);
		std::cout << "rebinning row: " << inputMatrix.rowIndex.back() << " goes into [" << newRowBinBorders.at(position) << "," << newRowBinBorders.at(position+1) << "]" << std::endl;
		for (unsigned int i=0; i<rowBinWidths.size();i++) {
			std::cout << "row " << i << " has width " << rowBinWidths.at(i) << std::endl;
		}
	}
	std::vector<double> columnBinWidths(newColumnBinBorders.size(), 0);
	std::vector<int> columnMapping;
	{
		double lower = 0;
		double higher = 0;
		for (unsigned int i=0; i<(inputMatrix.columnIndex.size()-1); i++) {
			lower = higher;
			higher = inputMatrix.columnIndex.at(i) + 0.5 * (inputMatrix.columnIndex.at(i+1) - inputMatrix.columnIndex.at(i));
			unsigned int position = exp(inputMatrix.columnIndex.at(i)/aRow)-1;
			columnBinWidths.at(position) += higher - lower;
			columnMapping.push_back(position);
			std::cout << "rebinning column: " << inputMatrix.columnIndex.at(i) << " goes into [" << newColumnBinBorders.at(position) << "," << newColumnBinBorders.at(position+1) << "]" << std::endl;
		}
		lower = higher;
		higher = inputMatrix.columnIndex.back() + 0.5 * (inputMatrix.columnIndex.back() - inputMatrix.columnIndex.at(inputMatrix.columnIndex.size()-2));
		unsigned int position = exp(inputMatrix.columnIndex.back()/aRow)-1;
		columnBinWidths.at(position) += higher - lower;
		columnMapping.push_back(position);
		std::cout << "rebinning column: " << inputMatrix.columnIndex.back() << " goes into [" << newColumnBinBorders.at(position) << "," << newColumnBinBorders.at(position+1) << "]" << std::endl;
		for (unsigned int i=0; i<columnBinWidths.size();i++) {
			std::cout << "row " << i << " has width " << columnBinWidths.at(i) << std::endl;
		}
	}

	BinnedTypedMatrix newMat(newRowBinBorders, newColumnBinBorders, mtDensity);
	for (unsigned int oldRow = 0; oldRow < inputMatrix.rowIndex.size(); oldRow++) {
		for (unsigned int oldColumn = 0; oldColumn < inputMatrix.columnIndex.size(); oldColumn++) {
			newMat.m(rowMapping[oldRow], columnMapping[oldColumn]) +=
					inputMatrix.m(oldRow, oldColumn)
							/ (rowBinWidths[rowMapping[oldRow]]
									* columnBinWidths[columnMapping[oldColumn]]);
		}
	}

	newMat.print();

	std::ofstream outFile(output);
	newMat.writeToFile(outFile);
	outFile.close();

    return 0;
}
