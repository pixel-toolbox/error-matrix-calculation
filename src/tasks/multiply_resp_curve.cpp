#include <fstream>
#include <iostream>
#include <experimental/filesystem>

#include "../tasks.hh"
#include "../util/split.hh"
#include "../classes/BinnedTypedMatrix.hh"

namespace fs = std::experimental::filesystem;

using namespace EMC;

using namespace boost::numeric;

int EMC::multiply_resp_curve(bool forceOverwrite, std::string inputVectorFile,
		std::string outputVectorFile, std::string respCurveFileName) {

	if (!forceOverwrite && fs::exists(outputVectorFile)) {
		std::cout << "Output file already exists, use -f to force overwrite." << std::endl;
		return 1;
	}

	{
		std::ifstream infile(inputVectorFile + ".step");
		std::ifstream infileCntr(inputVectorFile + ".cntr");
		std::ifstream respCurveFile(respCurveFileName);
		std::ofstream outputFile(outputVectorFile + ".step");
		std::ofstream outputFileCntr(outputVectorFile + ".cntr");

		double lastTransmissionEnergy = 0;
		while (!infile.eof()) {
			std::string lineInput;
			do {
				std::getline(infile, lineInput );
			} while (lineInput.empty() && !infile.eof());
			auto splitLineInput = split(split(split(lineInput, ';'), ','), ' ');
			double x1 = std::stod(splitLineInput.at(0));
			double y = std::stod(splitLineInput.at(1));

			do {
				std::getline(infile, lineInput );
			} while (lineInput.empty() && !infile.eof());
			splitLineInput = split(split(split(lineInput, ';'), ','), ' ');

			double x2 = std::stod(splitLineInput.at(0));

			do {
				std::getline(infileCntr, lineInput );
			} while (lineInput.empty() && !infileCntr.eof());
			splitLineInput = split(split(split(lineInput, ';'), ','), ' ');

			double xcntr = std::stod(splitLineInput.at(0));
			double yerr = std::stod(splitLineInput.at(2));

			std::vector<double> transmissionCurveValues;

			std::string lineResponse;
			double thisTransmissionCenter=0;
			do {
				std::getline(respCurveFile, lineResponse );
				try {
					auto splitLineSubtracted = split(split(split(lineResponse, ';'), ','), ' ');
					double value = std::stod(splitLineSubtracted.at(1));
					thisTransmissionCenter = std::stod(splitLineSubtracted.at(0));
					transmissionCurveValues.push_back(value);
					if ( (2*thisTransmissionCenter - lastTransmissionEnergy) > x2) {
						break;
					}
				} catch(...) {
					std::cout << "error at " << lineResponse << std::endl;
				}
			} while (!respCurveFile.eof());


			if (infile.eof()) {
				break;
			}

			double mult_val = 0;

			for (auto val : transmissionCurveValues) {
				mult_val += 1 / (val * (double)transmissionCurveValues.size());
			}

			try {
				outputFile << x1 << " " << (y*mult_val) << std::endl;
				outputFile << x2 << " " << (y*mult_val) << std::endl;
				outputFileCntr << xcntr << " " << (y*mult_val) << " " << (yerr*mult_val) << std::endl;
			} catch(...) { }
		}

		infile.close();
		infileCntr.close();
		respCurveFile.close();
		outputFile.close();
		outputFileCntr.close();
	}


    return 0;
}
