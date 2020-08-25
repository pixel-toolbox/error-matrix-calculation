#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <stdio.h>
#include <random>
#include <string>


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

/// from  https://stackoverflow.com/questions/440133/how-do-i-create-a-random-alpha-numeric-string-in-c
std::string random_string(std::string::size_type length) {
	static auto& chrs = "0123456789"
			"abcdefghijklmnopqrstuvwxyz"
			"ABCDEFGHIJKLMNOPQRSTUVWXYZ";

	thread_local static std::mt19937 rg { std::random_device { }() };
	thread_local static std::uniform_int_distribution<std::string::size_type> pick(
			0, sizeof(chrs) - 2);

	std::string s;

	s.reserve(length);

	while (length--)
		s += chrs[pick(rg)];

	return s;
}

unsigned int modeCount = 5;

void shiftCenters(std::vector<double>& inputVec, unsigned int mode, double error) {
	switch(mode) {
	case 1: {
		for (int i=0; i<inputVec.size(); i++) {
			inputVec.at(i) += error;
		}
		break;
	}
	case 2: {
		for (int i=0; i<inputVec.size(); i++) {
			inputVec.at(i) += error;
			if (inputVec.at(i) < 0) {
				inputVec.at(i) = 0;
			}
		}
		break;
	}
	case 3: {
		double half = inputVec.size() / 2.0;
		for (int i=0; i<inputVec.size(); i++) {
			inputVec.at(i) += error * (i-half) / inputVec.size();
			if (inputVec.at(i) < 0) {
				inputVec.at(i) = 0;
			}
		}
		break;
	}
	case 4: {
		double half = inputVec.size() / 2.0;
		for (int i=0; i<inputVec.size(); i++) {
			inputVec.at(i) -= error * (i-half) / inputVec.size();
			if (inputVec.at(i) < 0) {
				inputVec.at(i) = 0;
			}
		}
		break;
	}
	case 0: {
		/// do nothing, mode 0 does not modify bin centers
		break;
	}
	default:
		std::cout << "this should not happen!!" << std::endl;
		break;
	}
}

int EMC::rebin_invert(bool forceOverwrite, std::string input, std::string output,
		int countRate, double minRow, double minCol, double minRowBin,
		double minColBin, double maxRow, double maxCol, double rowError, double colError) {

	if (!forceOverwrite && fs::exists(output)) {
		std::cout << "Output file already exists, use -f to force overwrite." << std::endl;
		return 1;
	}

    std::streambuf* orig_buf = std::cout.rdbuf();
    std::cout.rdbuf(NULL);

	std::string randRebinName = output + "_" + random_string(5) + "_tmp.rebinned.mat";
	rebin(true, input, randRebinName, countRate, minRow, minCol, minRowBin, minColBin, maxRow, maxCol);

	std::ifstream infile(randRebinName);
	BinnedTypedMatrix inputMatrix = BinnedTypedMatrix::readFromFile(infile);
	infile.close();

	fs::remove(randRebinName);

	BinnedTypedMatrix newMat(inputMatrix.columnIndex, inputMatrix.rowIndex, inputMatrix.matType);
	{
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
	}

	std::vector<std::vector<std::vector<real>>> wiggleResults(newMat.rowCount, std::vector<std::vector<real>>(newMat.columnCount, std::vector<real>()));

	for (unsigned int i=0; i<modeCount; i++) {
		for (unsigned int j=0; j<modeCount; j++) {
			if (i==0 && j==0) {
				continue;
			}

		    std::cout.rdbuf(orig_buf);
			std::cout << "mode " << i << " " << j << std::endl;
			std::cout.flush();
		    std::cout.rdbuf(NULL);

			std::string randKey = random_string(5);

		    std::cout.rdbuf(orig_buf);
			std::cout << "shift, ";
			std::cout.flush();
		    std::cout.rdbuf(NULL);
			std::string tmpShiftedMatrixName = output + "_" + randKey + "_tmp.shifted.mat";
			{ // shift matrix by mode
				std::ifstream infile(input);
				BinnedTypedMatrix inputMat = BinnedTypedMatrix::readFromFile(infile);

				double mult = 0.7071067812;

				if (i==0 || j==0)
					mult=1;

				shiftCenters(inputMat.rowIndex, i, mult*rowError);
				shiftCenters(inputMat.columnIndex, j, mult*colError);

				std::ofstream shiftedMatOutputName(tmpShiftedMatrixName);
				inputMat.writeToFile(shiftedMatOutputName);
			}

		    std::cout.rdbuf(orig_buf);
			std::cout << "rebin, ";
			std::cout.flush();
		    std::cout.rdbuf(NULL);
			std::string tmpRebinnedMatName = output + "_" + randKey + "_tmp.rebin.mat";
			rebin(true, tmpShiftedMatrixName, tmpRebinnedMatName, countRate, minRow, minCol, minRowBin, minColBin, maxRow, maxCol);

		    std::cout.rdbuf(orig_buf);
			std::cout << "invert ";
			std::cout.flush();
		    std::cout.rdbuf(NULL);
			std::string tmpInvertedMatName = output + "_" + randKey + "_tmp.inverted.mat";
			invert(true, tmpRebinnedMatName, tmpInvertedMatName);

			std::ifstream infile(tmpInvertedMatName);
			BinnedTypedMatrix rebinnedInvertedMatrix = BinnedTypedMatrix::readFromFile(infile);
			infile.close();

			fs::remove(tmpShiftedMatrixName);
			fs::remove(tmpRebinnedMatName);
			fs::remove(tmpInvertedMatName);

		    std::cout.rdbuf(orig_buf);
			std::cout << "and done!" << std::endl;
			std::cout.flush();
		    std::cout.rdbuf(NULL);
			for (unsigned int rowPos = 0; rowPos < newMat.rowCount; rowPos++) {
				for (unsigned int colPos = 0; colPos < newMat.columnCount; colPos++) {
					wiggleResults.at(rowPos).at(colPos).push_back(rebinnedInvertedMatrix.m(rowPos, colPos).value);
				}
			}
		}
	}

	for (unsigned int rowPos = 0; rowPos < newMat.rowCount; rowPos++) {
		for (unsigned int colPos = 0; colPos < newMat.columnCount; colPos++) {
			newMat.m(rowPos, colPos).err_sq = 0;
			for (unsigned int i = 0; i<(modeCount*modeCount-1); i++) {

				// use this (quadratic mean)
				newMat.m(rowPos, colPos).err_sq += pow(wiggleResults.at(rowPos).at(colPos).at(i)-newMat.m(rowPos, colPos).value, 2) / (modeCount*modeCount-1);

				// or this (linear mean)
				/*newMat.m(rowPos, colPos).err_sq += std::fabs(wiggleResults.at(rowPos).at(colPos).at(i)-newMat.m(rowPos, colPos).value) / (maxMode*maxMode);
				newMat.m(rowPos, colPos).err_sq = newMat.m(rowPos, colPos).err_sq * newMat.m(rowPos, colPos).err_sq;*/
			}
		}
	}

    std::ofstream outputF(output);
    newMat.writeToFile(outputF);
    outputF.close();

    return 0;
}
