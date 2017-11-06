#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <Eigen\Dense>
#include <Eigen\Core>

namespace Optimum {
	/**
	* Loads data from a csv file.
	* @param filepath the absolute filepath to the csv file.
	* @param delimiter the delimiter of the csv file, typically a comma.
	* @return MatrixXd an Eigen matrix of unknown dimension.
	*/
	static Eigen::MatrixXd loadCsv(std::string filepath, char delimiter) {
 
		filepath.erase(remove(filepath.begin(), filepath.end(), '\"'), filepath.end());
		std::ifstream classFile(filepath);
		std::string line;
		std::string member;
		int rows = 0, cols = 0;
		std::vector<std::vector<std::string>> values;
		//first get the number of rows and columns. 
		while (getline(classFile, line)) {
			rows += 1;
			std::istringstream s(line);
			std::vector<std::string> vals;
			cols = 0;
			while (getline(s, member, delimiter)) {
				cols += 1;
				vals.push_back(member);
			}
			values.push_back(vals);
		}
		//std::cout << "Reading data: " << rows << " rows, " << cols << " columns." << std::endl;
		//set the data. 
		Eigen::MatrixXd data(rows, cols);
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				data(r, c) = stod(values[r][c]);
			}
		}

		return data;
	}

	/**
	* Writes a csv file given a filePath, a delimeter and the data.
	* @param filePath the file path to the file.
	* @param delimeter the delimeter for it.
	* @para data the data to write.
	*/
	static void writeCsv(std::string filePath, char delimeter, Eigen::MatrixXd data) {
		std::string delim;
		std::stringstream ss;
		ss << delimeter;
		ss >> delim;
		std::string findString(".csv");
		std::size_t found = filePath.find(findString);
		if (found == std::string::npos) {
			filePath.append(".csv");
		}

		std::ofstream outFile(filePath, std::ios::out);
		int rows = data.rows();
		int cols = data.cols();
		for (int r = 0; r < rows; r++) {
			std::ostringstream s;
			for (int c = 0; c < cols; c++) {
				double d = data(r, c);
				s << d;
				if (c < cols - 1){
					s << delim;
				}
			}
			outFile << s.str() << std::endl;
		}
		outFile.close();
		/*std::cout << "File successfully saved." << std::endl;*/
	}
}

#endif