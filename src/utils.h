#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <Eigen\Dense>
#include <Eigen\Core>

using namespace std;
using namespace Eigen;

/**
* Loads data from a csv file. 
* @param filepath the absolute filepath to the csv file.
* @param delimiter the delimiter of the csv file, typically a comma. 
* @return MatrixXd and Eigen matrix of unknown dimension. 
*/
MatrixXd loadCsv(std::string filepath, char delimiter) {
	
	//remove double quotes. 
	filepath.erase(remove(filepath.begin(), filepath.end(), '\"'), filepath.end());
	ifstream classFile(filepath);
	string line;
	string member;
	int rows = 0, cols = 0;
	vector<vector<string>> values;
	//first get the number of rows and columns. 
	while (getline(classFile, line)) {
		rows += 1;
		istringstream s(line);
		vector<string> vals;
		cols = 0;
		while (getline(s, member, delimiter)) {
			cols += 1;
			vals.push_back(member);
		}
		values.push_back(vals);
	}
	cout << "Reading data: " << rows << " rows, " << cols << " columns." << endl;
	//set the data. 
	MatrixXd data(rows, cols);
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
void writeCsv(string filePath, char delimeter, MatrixXd data) {
	string findString(".csv");
	std::size_t found = filePath.find(findString);
	if (found == string::npos) {
		filePath.append(".csv");
	}

	ofstream outFile;
	outFile.open(filePath);
	int rows = data.rows();
	int cols = data.cols();
	for (int r = 0; r < rows; r++) {
		std::ostringstream s;
		for (int c = 0; c < cols; c++) {
			double d = data(r, c);
			s << d;
			if (c < cols - 1){
				s << delimeter;
			}
		}
		outFile << s.str() << endl;
	}
	outFile.close();
	cout << "File successfully saved." << endl;
}
#endif