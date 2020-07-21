#ifndef EMC_tasks_hh_included
#define EMC_tasks_hh_included

#include "classes/structs.hh"

#include <iostream>
#include <vector>

namespace EMC {

int assemble(bool forceOverwrite, std::string input, std::string output,
		std::vector<double> inputBinCenters, matrixType matType, int countRate);

} /* namespace EMC */

#endif /* EMC_tasks_hh_included */
