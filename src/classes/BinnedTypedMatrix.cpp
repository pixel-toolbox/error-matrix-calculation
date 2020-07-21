#include "BinnedTypedMatrix.hh"

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

BinnedTypedMatrix BinnedTypedMatrix::readFromFile(std::ifstream inFile) {
	return BinnedTypedMatrix(std::vector<double>(), std::vector<double>(), mtNumber);
}

BinnedTypedMatrix BinnedTypedMatrix::makeProbability(int countRate) {
	if (matType != mtNumber) {
		std::cout << "cannot convert to mtProbability in BinnedTypedMatrix::makeProbability" << std::endl;
		exit(1);
	}
	BinnedTypedMatrix retMat(rowIndex, columnIndex, mtProbability);

	for (unsigned int i = 0; i<rowIndex.size(); i++) {
		for (unsigned int j = 0; j<columnIndex.size(); j++) {
			retMat.m(i, j) = m(i, j) / countRate;
		}
	}

	return retMat;
}

void BinnedTypedMatrix::print() {
	for (unsigned int i = 0; i<rowIndex.size(); i++) {
		for (unsigned int j = 0; j<columnIndex.size(); j++) {
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
	for (unsigned int i = 0; i<rowIndex.size(); i++) {
		for (unsigned int j = 0; j<columnIndex.size(); j++) {
			outFile << m(i, j).value << " ";
		}
		outFile << std::endl;
	}

	outFile << "error" << std::endl;
	for (unsigned int i = 0; i<rowIndex.size(); i++) {
		for (unsigned int j = 0; j<columnIndex.size(); j++) {
			outFile << sqrt(m(i, j).err_sq) << " ";
		}
		outFile << std::endl;
	}
}

} /* namespace EMC */
