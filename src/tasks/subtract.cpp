#include <fstream>
#include <iostream>
#include <experimental/filesystem>

#include "../tasks.hh"
#include "../util/split.hh"
#include "../classes/BinnedTypedMatrix.hh"

namespace fs = std::experimental::filesystem;

using namespace EMC;

using namespace boost::numeric;

int EMC::subtract(bool forceOverwrite, std::string inputVectorFile,
		std::string outputVectorFile, std::string spectrumToBeSubtracted) {

	if (!forceOverwrite && fs::exists(outputVectorFile)) {
		std::cout << "Output file already exists, use -f to force overwrite." << std::endl;
		return 1;
	}

	{
		std::ifstream infile(inputVectorFile + ".step");
		std::ifstream subtractedFile(spectrumToBeSubtracted + ".step");
		std::ofstream outputFile(outputVectorFile + ".step");

		while (!infile.eof()) {
			std::string lineInput;
			do {
				std::getline(infile, lineInput );
			} while (lineInput.empty() && !infile.eof());
			auto splitLineInput = split(split(split(lineInput, ';'), ','), ' ');

			std::string lineSubtracted;
			do {
				std::getline(subtractedFile, lineSubtracted );
			} while (lineSubtracted.empty() && !subtractedFile.eof());
			auto splitLineSubtracted = split(split(split(lineSubtracted, ';'), ','), ' ');

			if (infile.eof() || subtractedFile.eof()) {
				break;
			}

			if (std::stod(splitLineInput.at(0)) != std::stod(splitLineSubtracted.at(0))) {
				std::cout << "x-coordinate is not the same: " << splitLineInput.at(0) << "!=" << splitLineSubtracted.at(0) << std::endl;
				return 3;
			}

			double x = std::stod(splitLineInput.at(0));

			double y1 = std::stod(splitLineInput.at(1));
			double y2 = std::stod(splitLineSubtracted.at(1));

			try {
				outputFile << x << " " << (y1-y2) << std::endl;
			} catch(...) { }
		}

		infile.close();
		subtractedFile.close();
		outputFile.close();
	}

	{
		std::ifstream infile(inputVectorFile + ".cntr");
		std::ifstream subtractedFile(spectrumToBeSubtracted + ".cntr");
		std::ofstream outputFile(outputVectorFile + ".cntr");

		while (!infile.eof()) {
			std::string lineInput;
			do {
				std::getline(infile, lineInput );
			} while (lineInput.empty() && !infile.eof());
			auto splitLineInput = split(split(split(lineInput, ';'), ','), ' ');

			std::string lineSubtracted;
			do {
				std::getline(subtractedFile, lineSubtracted );
			} while (lineSubtracted.empty() && !subtractedFile.eof());
			auto splitLineSubtracted = split(split(split(lineSubtracted, ';'), ','), ' ');

			if (infile.eof() || subtractedFile.eof()) {
				break;
			}

			if (std::stod(splitLineInput.at(0)) != std::stod(splitLineSubtracted.at(0))) {
				std::cout << "x-coordinate is not the same: " << splitLineInput.at(0) << "!=" << splitLineSubtracted.at(0) << std::endl;
				return 4;
			}

			double x = std::stod(splitLineInput.at(0));

			double y1 = std::stod(splitLineInput.at(1));
			double y2 = std::stod(splitLineSubtracted.at(1));

			double y1_err = std::stod(splitLineInput.at(2));
			double y2_err = std::stod(splitLineSubtracted.at(2));

			try {
				outputFile << x << " " << (y1-y2) << " " << sqrt(y1_err*y1_err+y2_err*y2_err) << std::endl;
			} catch(...) { }
		}

		infile.close();
		subtractedFile.close();
		outputFile.close();
	}






    return 0;
}
