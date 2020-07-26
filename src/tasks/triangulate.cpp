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

/*int EMC::triangulate(bool forceOverwrite, std::string input, std::string output) {

	if (!forceOverwrite && fs::exists(output)) {
		std::cout << "Output file already exists, use -f to force overwrite." << std::endl;
		return 1;
	}

	std::ifstream infile(input);
	BinnedTypedMatrix inputMatrix = BinnedTypedMatrix::readFromFile(infile);
	infile.close();



    std::ofstream outputF(output);
    inputMatrix.writeToFile(outputF);
    outputF.close();

    return 0;
}
*/
