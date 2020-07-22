#include <iostream>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

#include "include.hh"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <map>
#include <string>

using namespace EMC;
using namespace std;

namespace EMC {

enum task {
    tNONE = 0,
    ASSEMBLE,
	REBIN,
    UNBINNED_SPECTRUM_ASSEMBLE,
    UNBINNED_ASSEMBLE,
    INVERT,
	APPLY
};

map<string, task> taskNameMap = {
    { "assemble", ASSEMBLE },
    { "rebin", REBIN },
    /*{ "unbinned_spectrum_assemble", UNBINNED_SPECTRUM_ASSEMBLE },
    { "unbinned_assemble", UNBINNED_ASSEMBLE },
    { "invert", INVERT},
    { "apply", APPLY }*/
};

} /* namespace EMC */

using namespace boost::numeric;

int main(int argc, char* argv[]) {
    cout << endl << "error-matrix-calculation v1 (c) Simon Michalke 2020" << endl;


    task todo = tNONE;
    bool forceOverwrite = false;
    string output;
    string input;

    std::vector<double> inputBinCenters;
    matrixType matType = mtNumber;
    int countRate = -1;

    bool fail = false;

	for (int argPos = 1; argPos < argc; argPos++) {
		if (argv[argPos][0] == '-') {
			switch (argv[argPos][1]) {
			case 't': // what task
				argPos++;
				if (taskNameMap.find(argv[argPos]) == taskNameMap.end()) {
					cout << "tasks not recognized!" << endl;
					fail = true;
				} else {
					todo = taskNameMap.at(argv[argPos]);
				}
				break;
			case 'f':
				forceOverwrite = true;
				break;
			case 'i': {
				argPos++;
				input = argv[argPos];
			}
				break;
			case 'o':
				argPos++;
				output = argv[argPos];
				break;
			case 'b': {
				argPos++;
				auto binVec = split(argv[argPos], ',');
				unsigned int i = 0;
				unsigned int binCount = binVec.size();
				double binIncrement = -1;
				for (;i<binVec.size();i++) {
					if (binVec[i] == "...") {
						binIncrement = inputBinCenters[i-1] - inputBinCenters[i-2];
						inputBinCenters.push_back(inputBinCenters[i-1] + binIncrement);
						binCount = std::stod(binVec[i+1]) / (inputBinCenters[i-1] - inputBinCenters[i-2]) - 1;
						break;
					}
					inputBinCenters.push_back(std::stod(binVec[i]));
				}
				for (;i<binCount;i++) {
					inputBinCenters.push_back(inputBinCenters[i] + binIncrement);
				}
			}
				break;
			case 'T':
				argPos++;
				if (matrixTypeMap.find(argv[argPos]) == matrixTypeMap.end()) {
					cout << "matrix type not recognized!" << endl;
					fail = true;
				} else {
					matType = matrixTypeMap.at(argv[argPos]);
				}
				break;
			case 'c':
				argPos++;
				countRate = std::stoi(argv[argPos]);
				break;
			default:
				cout << "unexpected flag in argument number " << argPos << endl;
				fail = true;
			}
		} else {
			cout << "expected \"-\" in argument number " << argPos << endl;
			fail = true;
		}
	}

    if (fail) {
        cout << "aborted" << endl;
        todo = tNONE;
    }

	switch (todo) {
	case ASSEMBLE: {
		cout << "assembling matrix" << endl;
		int retVal = assemble(forceOverwrite, input, output, inputBinCenters,
				matType, countRate);
		cout << "finished with status " << retVal << endl;
		break;
	}
	case REBIN: {
		cout << "rebinning matrix" << endl;
		int retVal = rebin(forceOverwrite, input, output, countRate);
		cout << "finished with status " << retVal << endl;
		break;
	}
	case tNONE:
    default:
        cout << "" << endl << endl;
        cout << " supported general arguments:" << endl;
        cout << "   -t [ assemble | invert | apply ]" << endl;
        cout << "      specifies what task is being executed" << endl;
        cout << "      should be specified first" << endl;
        cout << "   -f overwrite output files / directories, deletes target files before writing" << endl;
        cout << "   -i <filename>" << endl;
        cout << "      input file name" << endl;
        cout << "   -o <filename>" << endl;
        cout << "      output file name" << endl << endl;
        cout << "   In case of a vector the columns are expected to be bin center, value, value error" << endl;
        cout << "   with space character for spacing" << endl;
        cout << " arguments for specific tasks:" << endl;
        cout << "   assemble: assembles a matrix based on putting spectrums into columns, output bins determined by input spectrum" << endl;
        cout << "     -i <filename format string>" << endl;
        cout << "        input file name. String that will be formatted using sprintf with data passed by -b" << endl;
        cout << "        expects 'bin_center counts' in each line" << endl;
        cout << "     -b <bin0>,<bin1>,..." << endl;
        cout << "        list of input bin centers (floating point possible). Passed to -i via sprintf" << endl;
        cout << "        for evently spaced bins e.g. 1,2,...,100 is also possible" << endl;
        cout << "     -T <type>" << endl;
        cout << "        output type; possible values are 'number', 'probability' and 'density'" << endl;
        cout << "        in case of density bin width is determined based on the bin width in input files" << endl;
        cout << "     -c <integer>" << endl;
        cout << "        total events, count rate for normalization" << endl;
        cout << "   rebin: rebins a matrix by merging bins in order to increase statistical significance" << endl;
        cout << "     -c <integer>" << endl;
        cout << "        dimension of output square matrix" << endl;
        /*cout << "   unbinned_spectrum_assemble: assemble a matrix from unbinned spectrum data" << endl;
        cout << "     -i <filename format string>" << endl;
        cout << "        input file name. String that will be formatted using sprintf with data passed by -b" << endl;
        cout << "        expects an energy in each line" << endl;
        cout << "     -b <bin0>,<bin1>,..." << endl;
        cout << "        list of input bin centers (floating point possible). Passed to -i via sprintf" << endl;
        cout << "     -b O <bin0>,<bin1>,..." << endl;
        cout << "        list of output bin separators (floating point possible). Bin center will be placed in between." << endl;
        cout << "     -L" << endl;
        cout << "        Use logarithmic bin center" << endl;
        cout << "     -c <integer>" << endl;
        cout << "        total events, count rate for normalization" << endl;
        cout << "   unbinned_assemble: assemble a matrix from completely unbinned data" << endl;
        cout << "     -i <filename>" << endl;
        cout << "        input file name." << endl;
        cout << "        expects '<input energy> <output energy>' in each line" << endl;
        cout << "     -b I <bin0>,<bin1>,..." << endl;
        cout << "        list of input bin separators (floating point possible). Bin center will be placed in between." << endl;
        cout << "     -b O <bin0>,<bin1>,..." << endl;
        cout << "        list of output bin separators (floating point possible). Bin center will be placed in between." << endl;
        cout << "     -L" << endl;
        cout << "        Use logarithmic bin center" << endl;
        cout << "     -c <integer>" << endl;
        cout << "        total events, count rate for normalization" << endl;*/
        cout << endl;
    }

    /*ublas::vector<ValueError> v1(200);

    ublas::matrix<ValueError> m1 = ublas::identity_matrix<ValueError>(200);

    for (int i=0; i<200; i++) {
    	v1[i] = ValueError(i, sqrt(i));
    	m1(i, i) = ValueError(200-i, 42);
    }

    for (ValueError ve : ublas::prod(m1, v1)) {
    	cout << ve.value << " " << ve.err_sq << std::endl;
    }*/

    return 0;
}
