#ifndef SRC_CLASSES_BINNEDTYPEDMATRIX_HH_
#define SRC_CLASSES_BINNEDTYPEDMATRIX_HH_

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
			matType(matType), rowIndex(rowIndex), columnIndex(columnIndex), m(
					rowIndex.size(), columnIndex.size()) {
		std::cout << "creating BinnedTypedMatrix with " << rowIndex.size() << "x" << columnIndex.size() << std::endl;
	}

	matrixType matType;
	std::vector<double> rowIndex; // in case of density these are borders
	std::vector<double> columnIndex; // in case of density these are borders
	boost::numeric::ublas::matrix<ValueError> m;

	static BinnedTypedMatrix readFromFile(std::ifstream& inFile);
	BinnedTypedMatrix makeProbability(int countRate);
	void print();

	void writeToFile(std::ofstream& outFile);
};

} /* namespace EMC */

#endif /* SRC_CLASSES_BINNEDTYPEDMATRIX_HH_ */
