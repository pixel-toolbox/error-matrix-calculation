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
	REBIN_INVERT,
	DE_PILEUP,
    SUBTRACT,
	MULTIPLY_RESPONSE_CURVE,
    PRINT_STATS,
	APPLY,
	FIT_VALUES
};

map<string, task> taskNameMap = {
    { "assemble", ASSEMBLE },
    { "rebin", REBIN },
    { "invert", INVERT},
    { "rebin_invert", REBIN_INVERT},
    { "de_pileup", DE_PILEUP},
    { "subtract", SUBTRACT},
    { "multiply_resp_curve", MULTIPLY_RESPONSE_CURVE},
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

    float rowErrorEst = -1;
    float colErrorEst = -1;

    double lowerBound = 0;

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
			case 'l':
				argPos++;
				lowerBound = std::stof(argv[argPos]);
				break;
			case 'M':
				argPos++;
				matrixFileName = argv[argPos];
				break;
			case 'd': {
				argPos++;
				auto strVec = split(argv[argPos], ',');
				rowErrorEst = std::stod(strVec[0]);
				colErrorEst = std::stod(strVec[1]);
				break;
			}
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
	case REBIN_INVERT: {
		cout << "rebinning and inverting matrix" << endl;
		int retVal = rebin_invert(forceOverwrite, input, output, countRate, minRow,
				minCol, minRowBin, minColBin, maxRow, maxCol, rowErrorEst, colErrorEst);
		cout << "finished with status " << retVal << endl;
		break;
	}
	case DE_PILEUP: {
		cout << "reversing pileup" << endl;
		int retVal = de_pileup(forceOverwrite, input, output, countRate, lowerBound);
		cout << "finished with status " << retVal << endl;
		break;
	}
	case SUBTRACT: {
		cout << "subtracting matrix" << endl;
		int retVal = subtract(forceOverwrite, input, output, matrixFileName);
		cout << "finished with status " << retVal << endl;
		break;
	}
	case MULTIPLY_RESPONSE_CURVE: {
		cout << "multiply-response-curving matrix" << endl;
		int retVal = multiply_resp_curve(forceOverwrite, input, output, matrixFileName);
		cout << "finished with status " << retVal << endl;
		break;
	}
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
        cout << "   -t [ de_pileup | assemble | rebin | invert | rebin_invert | apply | print_stats | fit_values ]" << endl;
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
        cout << "   de_pileup: reconstructs pileup numerically" << endl;
        cout << "     -c <integer>" << endl;
        cout << "        total pixel count" << endl;
        cout << "     -l <float>" << endl;
        cout << "        up to this energy counts are set to 0" << endl;
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
        cout << "        row and column bin parameters, float" << endl;
        cout << "   rebin_invert: rebins a matrix and inverts it. Estimates errors as well." << endl;
        cout << "     -c <integer>" << endl;
        cout << "        dimension of output square matrix" << endl;
        cout << "     -m <minRow><minCol><minRowBin><minColBin><maxRow><maxCol>" << endl;
        cout << "        row and column bin parameters, float" << endl;
        cout << "     -d <row Error, float>,<column Error, float>" << endl;
        cout << "        estimate error based on row and column error." << endl;
        cout << "   invert: inverts a matrix" << endl;
        cout << "   subtract: subtracts spectra" << endl;
        cout << "        uses <filelane>.cntr and <filelane>.step for input, output and filename of subtracted spectrum." << endl;
        cout << "     -M <filename>" << endl;
        cout << "        filename of the spectrum that will be subtracted" << endl;
        cout << "   multiply_resp_curve: scales spectrum with response curve" << endl;
        cout << "        uses <filelane>.cntr and <filelane>.step for input and output" << endl;
        cout << "     -M <filename>" << endl;
        cout << "        filename of the transmission curve" << endl;
        cout << "        this is assumed to be more dense than the input spectrum" << endl;
        cout << "   print_stats: inverts a matrix" << endl;
        cout << "     -M <filename>" << endl;
        cout << "        filename of the matriy" << endl;
        cout << "   apply: inverts a matrix; input and output are spectrum files" << endl;
        cout << "        first column is expected to be bin center and second the counts in that bin" << endl;
        cout << "        vector will be rebinned and rescaled to matrix dimensions" << endl;
        cout << "        output is bin center and count rate" << endl;
        cout << "     -M <filename>" << endl;
        cout << "        filename of the matrix" << endl;
        cout << "   fit_values: optimizes rebinning based on inverted matrix error metric" << endl;
        cout << "               output is the base filename for the rebinned (.rebin.mat) and inverted (.rebin.inverted.mat) matrix" << endl;
        cout << "               these matrices will be overwritten multiple times!" << endl;
        cout << "     -c <integer>" << endl;
        cout << "        dimension of output square matrix" << endl;
        cout << "     -m <minRow><minCol><minRowBin><minColBin><maxRow><maxCol>" << endl;
        cout << "        initial row and column bin parameters; minCol and maxCol will stay constant" << endl;
        cout << endl;
    }

    return 0;
}
