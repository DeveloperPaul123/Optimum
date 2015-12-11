#include "optimizer.h"
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <algorithm>

using namespace std;
using namespace Eigen;
float f(Vector &vec, MatrixXd &cbct, MatrixXd &bli);
void testTransform();
float rosenbrock(Vector v);
float camelHumps(Vector v);

int main(int argc, const char* argv[]) {

	//test data. 
	float precision = 0.00001;
	int dimension = 2;

	NelderMeadMinimizer minimizer(dimension, precision);
	//really bad start values. 
	Vector v(-10, 30);
	minimizer.initialGuess(v);
	while (!minimizer.done()) {
		float score = rosenbrock(v);
		cout << "Score: " << score << endl;
		v = minimizer.step(v, score);
	}

	cout << "X : " << v[0] << " Y: " << v[1] << endl;

	cout << "Press any key to exit...";
	cin.get();
	return 0;

}

void testTransform() {
	//going to try to map coordinates from one system to the coordinate system of another
	//using the NedlerMead as a minimizing function. 
	double arrA[21][3] = {
		{ -16.6728, 25.9217, 1.9945 },
		{ -1.60935, 26.4446, 1.9945 },
		{ 13.4151, 26.8788, 1.9945 },
		{ -31.0609, 10.4577, 1.9945 },
		{ -16.3166, 11.0408, 1.9945 },
		{ -1.08825, 11.4183, 1.9945 },
		{ 13.8795, 12.1184, 1.9945 },
		{ 28.7055, 12.581, 1.9945 },
		{ -30.5346, -4.51181, 1.9945 },
		{ -15.6519, -4.0138, 1.9945 },
		{ -0.625645, -3.46434, 1.9945 },
		{ 14.4307, -3.0567, 1.9945 },
		{ 29.2584, -2.50372, 1.9945 },
		{ -30.1269, -19.6267, 1.9945 },
		{ -15.2443, -19.0152, 1.9945 },
		{ -0.101031, -18.4339, 1.9945 },
		{ 14.8118, -17.9695, 1.9945 },
		{ 29.8079, -17.53, 1.9945 },
		{ -14.6948, -33.8961, 1.9945 },
		{ 0.421824, -33.4318, 1.9945 },
		{ 15.3027, -32.9089, 1.9945 }
	};

	double arrB[21][3] = {
		{ 848.555, 745.678, 0 },
		{ 1006.32, 744.568, 0 },
		{ 1160.63, 744.322, 0 },
		{ 696.538, 903.447, 0 },
		{ 850.82, 901.093, 0 },
		{ 1007.46, 897.763, 0 },
		{ 1164.11, 897.875, 0 },
		{ 1317.44, 895.364, 0 },
		{ 699.913, 1060.06, 0 },
		{ 853.107, 1057.8, 0 },
		{ 1009.74, 1054.31, 0 },
		{ 1166.4, 1052.07, 0 },
		{ 1318.55, 1050.84, 0 },
		{ 703.31, 1215.68, 0 },
		{ 856.594, 1212.19, 0 },
		{ 1013.23, 1209.88, 0 },
		{ 1168.82, 1207.62, 0 },
		{ 1320.82, 1203.09, 0 },
		{ 859.857, 1366.56, 0 },
		{ 1014.5, 1363.19, 0 },
		{ 1168.67, 1361.94, 0 }
	};

	MatrixXd A(21, 2);
	MatrixXd B(21, 2);

	for (int r = 0; r < 21; r++) {
		A(r, 0) = arrA[r][0];
		A(r, 1) = arrA[r][1];
		double x = arrB[r][0];
		double y = arrB[r][1];

		double correctedX = (x - (0.5 * 2048))*0.099;
		double correctedY = ((0.5 * 2048) - y) * 0.099;
		B(r, 0) = correctedX;
		B(r, 1) = correctedY;
	}

	float precision = 0.0001;
	int dimension = 3;

	NelderMeadMinimizer minimizer(dimension, precision);
	Vector v(0, 0, 0);
	minimizer.initialGuess(v);
	while (!minimizer.done()) {
		float score = f(v, A, B);
		cout << "Score: " << score << endl;
		v = minimizer.step(v, score);
		cout << v[0] << " " << v[1] << " " << v[2] << endl;
	}

	cout << "X: " << v[0] << " mm ";
	cout << "Y: " << v[1] << " mm ";
	cout << "Deg: " << v[2] << " degree " << endl;

}

float f(Vector &vec, MatrixXd &cbct, MatrixXd &bli) {
	Vector temp = Vector(vec);
	float shift_x = vec[0];
	float shift_y = vec[1];
	float deg = vec[2];
	MatrixXd ones = MatrixXd::Ones(21, 1);
	MatrixXd shift(21, 2);
	shift << shift_x*ones, shift_y*ones;

	Matrix2d trans_rot;
	double degreeToRad = 3.14159 / 180.0;
	trans_rot << cos((deg * degreeToRad)), -sin((deg * degreeToRad)), 
	sin((deg*degreeToRad)), cos((deg*degreeToRad));
	

	MatrixXd trans = (bli + shift)*trans_rot;
	//MatrixXd trans = (bli + shift);

	float sumErr = 0.0;
	for (int i = 0; i < 21; i++) {
		double xDiff = cbct(i, 0) - trans(i, 0);
		double yDiff = cbct(i, 1) - trans(i, 1);
		double sqX = xDiff * xDiff;
		double sqY = yDiff * yDiff;
		double sum = sqX + sqY;
		double sqt = sqrt(sum);
		sumErr += sqt;
	}

	return sumErr;
}

/**
* Two local maxima at -1.5, 0 and 1.5 , 0;
* @param v the vector with the x y coordinates. 
*/
float camelHumps(Vector v) {
	float x = v[0];
	float y = v[1];
	return ((-x*x*x*x + 4.5*x*x + 2.0) / pow(2.71828, 2.0 * y*y));
}

/**
* Rosenbrock function. Defined as f(x, y) = (1-x)^2 + alpha * (y-x^2)^2.
* @param v vector with x y point.
* @return float the value of the function at the given point. 
*/
float rosenbrock(Vector v) {
	double x = v[0];
	double y = v[1];

	double alpha = 10.0;
	return (1 - x)*(1 - x) + alpha* ((y - (x*x)) *(y - (x*x)));
}