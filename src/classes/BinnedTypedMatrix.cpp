#include "BinnedTypedMatrix.hh"
#include "../util/split.hh"

#include <iostream>
#include <iomanip>

namespace EMC {

std::map<std::string, matrixType> matrixTypeMap = {
    { "number", mtNumber },
    { "probability", mtProbability },
    { "density", mtDensity }
};

std::map<matrixType, std::string> invertedMatrixTypeMap = {
    { mtNumber, "number"},
    { mtProbability, "probability"},
    { mtDensity, "density"}
};

BinnedTypedMatrix BinnedTypedMatrix::readFromFile(std::ifstream& inFile) {
	matrixType matType = mtNone;
	std::vector<double> rowIndex; // in case of density these are borders
	std::vector<double> columnIndex; // in case of density these are borders

	while (!inFile.eof()) {
		std::string line;
		std::getline(inFile, line);
		auto splitLine = split(line, ' ');
		try {
			std::string keyword = splitLine[0];
			if (keyword == "type") {
				matType = matrixTypeMap.at(splitLine[1]);
			} else if (keyword == "rowIndex") {
				for (unsigned int i = 1; i < splitLine.size(); i++) {
					rowIndex.push_back(std::stod(splitLine[i]));
				}
			} else if (keyword == "columnIndex") {
				for (unsigned int i = 1; i < splitLine.size(); i++) {
					columnIndex.push_back(std::stod(splitLine[i]));
				}
			}
		} catch (...) {
			std::cout << "error at line: " << line << std::endl;
		}
		if (matType != mtNone && rowIndex.size() > 0
				&& columnIndex.size() > 0) {
			break;
		}
	}

	BinnedTypedMatrix retMat(rowIndex, columnIndex, matType);

	unsigned int columnCount;
	unsigned int rowCount;
	if (matType == mtDensity) {
		columnCount = columnIndex.size() - 1;
		rowCount = rowIndex.size() - 1;
	} else {
		columnCount = columnIndex.size();
		rowCount = rowIndex.size();
	}

	while (!inFile.eof()) {
		std::string line;
		std::getline(inFile, line);
		if (line != "") {
			try {
				if (line == "value") {
					for (unsigned int i=0; i<rowCount; i++) {
						std::getline(inFile, line);
						auto splitLine = split(line, ' ');
						for (unsigned int j=0; j<columnCount; j++) {
							retMat.m(i, j).value = std::stod(splitLine[j]);
						}
					}
				} else if (line == "error") {
					for (unsigned int i=0; i<rowCount; i++) {
						std::getline(inFile, line);
						auto splitLine = split(line, ' ');
						for (unsigned int j=0; j<columnCount; j++) {
							retMat.m(i, j).err_sq = pow(std::stod(splitLine[j]), 2);
						}
					}
				} else {
					std::cout << "error in BinnedTypedMatrix::readFromFile while reading mat values" << std::endl;
				}
			} catch (...) {
				std::cout << "error at line: " << line << std::endl;
			}
		}
	}

	return retMat;
}

BinnedTypedMatrix BinnedTypedMatrix::makeProbability(int countRate) {
	if (matType != mtNumber) {
		std::cout << "cannot convert to mtProbability in BinnedTypedMatrix::makeProbability" << std::endl;
		exit(1);
	}
	BinnedTypedMatrix retMat(rowIndex, columnIndex, mtProbability);

	for (unsigned int i = 0; i<rowCount; i++) {
		for (unsigned int j = 0; j<columnCount; j++) {
			retMat.m(i, j) = m(i, j) / countRate;
		}
	}

	return retMat;
}

void BinnedTypedMatrix::print() {
	for (unsigned int i = 0; i<rowCount; i++) {
		for (unsigned int j = 0; j<columnCount; j++) {
			std::cout << "thisMatrix.m(" << i << ", " << j << ") = (" << m(i, j).value << ", " << sqrt(m(i, j).err_sq) << ")" << std::endl;
		}
	}
}

void BinnedTypedMatrix::writeToFile(std::ofstream& outFile) {

	outFile << "type " << invertedMatrixTypeMap.at(matType) << std::endl;

	outFile << "rowIndex ";
	for (auto s : rowIndex) {
		outFile << s << " ";
	}
	outFile << std::endl;

	outFile << "columnIndex ";
	for (auto s : columnIndex) {
		outFile << s << " ";
	}
	outFile << std::endl;

	outFile << "value" << std::endl;
	for (unsigned int i = 0; i<rowCount; i++) {
		for (unsigned int j = 0; j<columnCount; j++) {
			outFile << std::setw(12) << m(i, j).value << " ";
		}
		outFile << std::endl;
	}

	outFile << "error" << std::endl;
	for (unsigned int i = 0; i<rowCount; i++) {
		for (unsigned int j = 0; j<columnCount; j++) {
			outFile << std::setw(12) << sqrt(m(i, j).err_sq) << " ";
		}
		outFile << std::endl;
	}
}

} /* namespace EMC */
