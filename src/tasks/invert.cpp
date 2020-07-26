#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <stdio.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "../tasks.hh"
#include "../classes/BinnedTypedMatrix.hh"

namespace fs = std::experimental::filesystem;

using namespace EMC;

using namespace boost::numeric::ublas;

int EMC::invert(bool forceOverwrite, std::string input, std::string output) {

	if (!forceOverwrite && fs::exists(output)) {
		std::cout << "Output file already exists, use -f to force overwrite." << std::endl;
		return 1;
	}

	std::ifstream infile(input);
	BinnedTypedMatrix inputMatrix = BinnedTypedMatrix::readFromFile(infile);
	infile.close();

	BinnedTypedMatrix newMat(inputMatrix.columnIndex, inputMatrix.rowIndex, inputMatrix.matType);

	typedef permutation_matrix<std::size_t> pmatrix;
	// create a working copy of the input
	matrix<ValueError> A(inputMatrix.m);
	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());
	// perform LU-factorization
	int res = lu_factorize(A,pm);
	if( res != 0 )
		return false;
	// create identity matrix of "inverse"
	newMat.m.assign(identity_matrix<ValueError>(A.size1()));
	// backsubstitute to get the inverse
	lu_substitute(A, pm, newMat.m);

    std::ofstream outputF(output);
    newMat.writeToFile(outputF);
    outputF.close();

    return 0;
}
