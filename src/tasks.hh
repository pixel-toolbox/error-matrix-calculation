#ifndef EMC_tasks_hh_included
#define EMC_tasks_hh_included

#include "classes/structs.hh"

#include <iostream>
#include <vector>

namespace EMC {

int assemble(bool forceOverwrite, std::string input, std::string output,
		std::vector<double> inputBinCenters, matrixType matType, int countRate);

int rebin(bool forceOverwrite, std::string input, std::string output,
		int countRate, double minRow, double minCol, double minRowBin,
		double minColBin, double maxRow, double maxCol);

int invert(bool forceOverwrite, std::string input, std::string output);

int rebin_invert(bool forceOverwrite, std::string input, std::string output,
		int countRate, double minRow, double minCol, double minRowBin,
		double minColBin, double maxRow, double maxCol, double rowError, double colError);

//int triangulate(bool forceOverwrite, std::string input, std::string output);

int print_stats(std::string matrixFileName);

int apply(bool forceOverwrite, std::string inputVectorFile,
		std::string outputVectorFile, std::string matrixFileName);

int subtract(bool forceOverwrite, std::string inputVectorFile,
		std::string outputVectorFile, std::string spectrumToBeSubtracted);

int multiply_resp_curve(bool forceOverwrite, std::string inputVectorFile,
		std::string outputVectorFile, std::string respCurveFile);

int fit_values(std::string input, std::string output, int countRate,
		double minRow, double minCol, double minRowBin, double minColBin,
		double maxRow, double maxCol);

int de_pileup(bool forceOverwrite, std::string input, std::string output, int pixelCount, double lowerBound);

} /* namespace EMC */

#endif /* EMC_tasks_hh_included */
