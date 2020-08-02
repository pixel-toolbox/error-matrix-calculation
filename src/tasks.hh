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

//int triangulate(bool forceOverwrite, std::string input, std::string output);

int print_stats(std::string matrixFileName);

int apply(bool forceOverwrite, std::string inputVectorFile,
		std::string outputVectorFile, std::string matrixFileName);

int fit_values(std::string input, std::string output, int countRate,
		double minRow, double minCol, double minRowBin, double minColBin,
		double maxRow, double maxCol);

} /* namespace EMC */

#endif /* EMC_tasks_hh_included */
