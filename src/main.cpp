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
    INVERT,
    TRIANGULATE,
    PRINT_STATS,
	APPLY,
	FIT_VALUES
};

map<string, task> taskNameMap = {
    { "assemble", ASSEMBLE },
    { "rebin", REBIN },
    { "invert", INVERT},
    { "triangulate", TRIANGULATE},
    { "print_stats", PRINT_STATS},
    { "apply", APPLY },
    { "fit_values", FIT_VALUES }
};

} /* namespace EMC */

using namespace boost::numeric;

int main(int argc, char* argv[]) {
    cout << endl << "error-matrix-calculation v1 (c) Simon Michalke 2020" << endl;


    task todo = tNONE;
    bool forceOverwrite = false;
    string output;
    string input;

    string matrixFileName;

    std::vector<double> inputBinCenters;
    matrixType matType = mtNumber;
    int countRate = -1;

    double minRow = -1;
    double minCol = -1;
    double minRowBin = -1;
    double minColBin = -1;
    double maxRow = -1;
    double maxCol = -1;

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
			case 'm': {
				argPos++;
				auto maxStrVec = split(argv[argPos], ',');
				minRow = std::stod(maxStrVec[0]);
				minCol = std::stod(maxStrVec[1]);
				minRowBin = std::stod(maxStrVec[2]);
				minColBin = std::stod(maxStrVec[3]);
				maxRow = std::stod(maxStrVec[4]);
				maxCol = std::stod(maxStrVec[5]);
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
			case 'M':
				argPos++;
				matrixFileName = argv[argPos];
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
		int retVal = rebin(forceOverwrite, input, output, countRate, minRow,
				minCol, minRowBin, minColBin, maxRow, maxCol);
		cout << "finished with status " << retVal << endl;
		break;
	}
	case INVERT: {
		cout << "inverting matrix" << endl;
		int retVal = invert(forceOverwrite, input, output);
		cout << "finished with status " << retVal << endl;
		break;
	}
	/*case TRIANGULATE: {
		cout << "triangulating matrix" << endl;
		int retVal = triangulate(forceOverwrite, input, output);
		cout << "finished with status " << retVal << endl;
		break;
	}*/
	case PRINT_STATS: {
		cout << "printing matrix stats matrix" << endl;
		int retVal = print_stats(matrixFileName);
		cout << "finished with status " << retVal << endl;
		break;
	}
	case APPLY: {
		cout << "applying matrix" << endl;
		int retVal = apply(forceOverwrite, input, output, matrixFileName);
		cout << "finished with status " << retVal << endl;
		break;
	}
	case FIT_VALUES: {
		cout << "fitting matrix parameters" << endl;
		int retVal = fit_values(input, output, countRate, minRow,
				minCol, minRowBin, minColBin, maxRow, maxCol);
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
        cout << "     -m <minRow><minCol><minRowBin><minColBin><maxRow><maxCol>" << endl;
        cout << "        row and column bin parameters" << endl;
        cout << "   invert: inverts a matrix" << endl;
        cout << "   print_stats: inverts a matrix" << endl;
        cout << "     -M <filename>" << endl;
        cout << "        filename of the matriy" << endl;
        cout << "   apply: inverts a matrix; input and output are spectrum files" << endl;
        cout << "        first column is expected to be bin center and second the counts in that bin" << endl;
        cout << "        vector will be rebinned and rescaled to matrix dimensions" << endl;
        cout << "        output is bin center and count rate" << endl;
        cout << "     -M <filename>" << endl;
        cout << "        filename of the matriy" << endl;
        cout << "   fit_values: optimizes rebinning based on inverted matrix error metric" << endl;
        cout << "               output is the base filename for the rebinned (.rebin.mat) and inverted (.rebin.inverted.mat) matrix" << endl;
        cout << "               these matrices will be overwritten multiple times!" << endl;
        cout << "     -c <integer>" << endl;
        cout << "        dimension of output square matrix" << endl;
        cout << "     -m <minRow><minCol><minRowBin><minColBin><maxRow><maxCol>" << endl;
        cout << "        initial row and column bin parameters; minCol and maxCol will stay constant" << endl;
        //cout << "   triangulate: triangulate a matrix by deleting off-diagonal values" << endl;
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
