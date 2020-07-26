#ifndef CLASSES_BINNEDTYPEDMATRIX_HH_
#define CLASSES_BINNEDTYPEDMATRIX_HH_

#include "structs.hh"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <vector>
#include <fstream>

namespace EMC {

extern std::map<std::string, matrixType> matrixTypeMap;

struct BinnedTypedMatrix {
	inline BinnedTypedMatrix(std::vector<double> rowIndex,
			std::vector<double> columnIndex, matrixType matType) :
			matType(matType),
			rowCount(matType == mtDensity ?
				rowIndex.size() - 1 : rowIndex.size()),
			columnCount(matType == mtDensity ?
				columnIndex.size() - 1 : columnIndex.size()),
			rowIndex(rowIndex), columnIndex(columnIndex), m(rowCount, columnCount) {
		std::cout << "creating BinnedTypedMatrix with " << rowCount << "x" << columnCount << std::endl;
	}

	matrixType matType;
	unsigned int rowCount;
	unsigned int columnCount;
	std::vector<double> rowIndex; // in case of density these are borders
	std::vector<double> columnIndex; // in case of density these are borders
	boost::numeric::ublas::matrix<ValueError> m;

	static BinnedTypedMatrix readFromFile(std::ifstream& inFile);
	BinnedTypedMatrix makeProbability(int countRate);
	void print();

	void writeToFile(std::ofstream& outFile);
};

} /* namespace EMC */

#endif /* CLASSES_BINNEDTYPEDMATRIX_HH_ */
