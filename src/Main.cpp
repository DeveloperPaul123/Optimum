#include "neldermead.h"
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <algorithm>
#include "geneticalgorithm.h"
#include "utils.h"
#include "icp.h"

using namespace std;
using namespace Eigen;
float f(Chromosome<float>, MatrixXd &cbct, MatrixXd &bli);
float f(Vector &vec, MatrixXd &cbct, MatrixXd &bli);
void testTransform();
void testGATransform();
void testNM1D();
void testGA1D();
void testGA2D();
float rosenbrock(Vector v);
float gaRosenbrock(Chromosome<float> c);
float ga1DTest(Chromosome<float> c);
float nm1DTest(Vector v);
float camelHumps(Vector v);
float gaBoothsFunctionTest(Chromosome<float> c);

int main(int argc, const char* argv[]) {

	cout << "Testing ICP..." << endl;
	cout << "Enter reference points path: ";
	char filepath[256];
	cin.getline(filepath, 256);

	std::string fp = std::string(filepath);
	MatrixXd cbct_raw = loadCsv(fp, ',');

	cout << "Enter target data path: ";
	char filepathTarget[256];
	cin.getline(filepathTarget, 256);

	std::string fpTar = std::string(filepathTarget);
	MatrixXd bli_raw = loadCsv(fpTar, ',');

	MatrixXd cbct = cbct_raw.middleCols(1, 2);

	MatrixXd bli = bli_raw.middleCols(1, 2);

	//correct the data. 
	for (int r = 0; r < 21; r++) {
		double x = bli(r, 0);
		double y = bli(r, 1);
		double correctedX = (x - (0.5 * 2048))*0.0971;
		double correctedY = ((0.5 * 2048) - y) * 0.0971;
		bli(r, 0) = correctedX;
		bli(r, 1) = correctedY;
	}


	cout << "CBCT: " << endl;
	cout << cbct << endl;

	cout << "BLI: " << endl;
	cout << bli << endl; cout << endl;

	MatrixXd rotation = OptimalPointMatcher::solveForOptimalRotation(cbct, bli);
	MatrixXd trans = OptimalPointMatcher::solveForOptimalTranslation(cbct, bli, rotation);

	cout << "Rotation: " << endl;
	cout << rotation << endl; cout << endl;
	cout << "Translation: " << endl;
	cout << trans << endl; cout << endl;

	MatrixXd bliNew = OptimalPointMatcher::applyTransformation(bli,(trans * -1.0), rotation);
	cout << "Transformed:" << endl;
	cout << bliNew << endl; cout << endl;
	double rmse = OptimalPointMatcher::RMSE(bliNew, cbct);
	cout << "Error: " << rmse << endl;
	cin.get();

	cout << "Press any key to continue." << endl;
	cin.get();

	testNM1D();

	cout << "Press any key to conitue to GA 1D test..." << endl;
	cin.get();

	testGA1D();

	cout << "Press any key to continue to a transform test..." << endl;
	cin.get();

	testTransform();

	cout << "Press any key to continue to a GA Transform test..." << endl;
	cin.get();

	testGATransform();

	cout << "Press any key to continue to Rosenbrock test..." << endl;
	cin.get();

	//test data. 
	float precision = 0.00001;
	int dimension = 2;

	NelderMeadMinimizer minimizer(dimension, precision);
	//really bad start values. 
	Vector v(5, 8);
	minimizer.initialGuess(v);
	while (!minimizer.done()) {
		float score = rosenbrock(v);
		cout << "Score: " << score << endl;
		v = minimizer.step(v, score);
	}

	cout << "X : " << v[0] << " Y: " << v[1] << " Iterations: " << minimizer.getIterations() << endl;

	cout << "Press any key to continue to the GA Booth's Function Test..." << endl;
	cin.get();

	testGA2D();

	cout << "Press any key to exit...";
	cin.get();

	return 0;

}

void testGA2D() {
	float initialGuess[] = { 3.0f, 5.0f};
	Chromosome<float> chrome(initialGuess, 2);

	GeneticAlgorithm<float> ga;
	Population<float> mPop = ga.generatePopulation(chrome);
	for (int i = 0; i < mPop.size(); i++) {
		Chromosome<float> member = mPop.getMember(i);
		member.setFitness(gaBoothsFunctionTest(member));
		mPop.setMember(i, member);
	}

	//genetic algorithm test. 
	Chromosome<float> best = mPop.getMember(0);
	while (!ga.done()) {
		mPop.set(ga.step(mPop));
		for (int i = 0; i < mPop.size(); i++) {
			Chromosome<float> member = mPop.getMember(i);
			member.setFitness(gaBoothsFunctionTest(member));
			mPop.setMember(i, member);
		}
		mPop.bubbleSort();
		if (mPop.getMember(0).getFitness() < best.getFitness()) {
			best = mPop.getMember(0);
			cout << "Fitness: " << best.getFitness() << " Value: " << best.at(0) << " " << best.at(1) << endl;
		}
	}
	cout << "Best: " << best.at(0) << " " << best.at(1) << " Score: " << best.getFitness() << endl;
}

void testGA1D() {
	float initialGuess[] = { 3.0f };
	Chromosome<float> chrome(initialGuess, 1);

	GeneticAlgorithm<float> ga(150, 300, 0.1f, 0.40f);
	Population<float> mPop = ga.generatePopulation(chrome);
	for (int i = 0; i < mPop.size(); i++) {
		Chromosome<float> member = mPop.getMember(i);
		member.setFitness(ga1DTest(member));
		mPop.setMember(i, member);
	}

	//genetic algorithm test. 
	Chromosome<float> best = mPop.getMember(0);
	while (!ga.done()) {
		mPop.set(ga.step(mPop));
		for (int i = 0; i < mPop.size(); i++) {
			Chromosome<float> member = mPop.getMember(i);
			member.setFitness(ga1DTest(member));
			mPop.setMember(i, member);
		}
		mPop.bubbleSort();
		if (mPop.getMember(0).getFitness() < best.getFitness()) {
			best = mPop.getMember(0);
			cout << "Fitness: " << mPop.getMember(0).getFitness() << " Value: " << mPop.getMember(0).at(0) << endl;
		}
	}
	mPop.bubbleSort();
	best = mPop.getMember(0);
	//true min is x = 3. 
	cout << "Best: " << best.at(0) << " Score: " << best.getFitness() << endl;
	cout << "Error: " << abs(best.at(0) - 3.0f) / 3.0f << endl;
}

void testNM1D() {
	float precision = 0.000001;
	int dimension = 1;

	NelderMeadMinimizer minimizer(dimension, precision);
	Vector v(527.0f);
	minimizer.initialGuess(v);
	while (!minimizer.done()) {
		float score = nm1DTest(v);
		cout << "Score: " << score << endl;
		v = minimizer.step(v, score);
	}
	cout << "Best: " << v.at(0) << " Error: " << abs(v.at(0) - 3.0f) / 3.0f << " Iterations: " << minimizer.getIterations() << endl;
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
	Vector v(1.0f, 1.0f, 1.0f);
	minimizer.initialGuess(v);
	while (!minimizer.done()) {
		float score = f(v, A, B);
		cout << "Score: " << score << endl;
		v = minimizer.step(v, score);
	}

	cout << "X: " << v[0] << " mm ";
	cout << "Y: " << v[1] << " mm ";
	cout << "Deg: " << v[2] << " degree " << endl;

}

void testGATransform() {
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

	float guess[] = { 1.0f, 1.0f, 1.0f };
	Chromosome<float> chrome(guess, 3);

	GeneticAlgorithm<float> ga(500, 1000, 0.09, 0.30f);
	Population<float> pop = ga.generatePopulation(chrome);

	for (int i = 0; i < pop.size(); i++) {
		Chromosome<float> member = pop.getMember(i);
		member.setFitness(f(member, A, B));
		pop.setMember(i, member);
	}

	Chromosome<float> best = pop.getMember(0);
	while (!ga.done()) {
		pop.set(ga.step(pop));
		for (int i = 0; i < pop.size(); i++) {
			Chromosome<float> member = pop.getMember(i);
			member.setFitness(f(member, A, B));
			pop.setMember(i, member);
		}
		pop.bubbleSort();
		
		if (pop.getMember(0).getFitness() < best.getFitness()) {
			best = pop.getMember(0);
			cout << "Fitness: " << best.getFitness() << " " << best.at(0) << " " << best.at(1) << " " << best.at(2) << endl;
		}
	}
	cout << "Best: " << best.at(0) << " " << best.at(1) << " " << best.at(2) << " Score: " << best.getFitness() << endl;
}

/**
* 2D translation and rotation of one set of coordinates to another set of coordinates.
*/
float f(Vector &vec, MatrixXd &cbct, MatrixXd &bli) {
	Vector temp = Vector(vec);
	float shift_x = vec[0];
	float shift_y = vec[1];
	float deg = vec[2];
	MatrixXd ones = MatrixXd::Ones(21, 1);
	MatrixXd shift(2, 1);
	shift(0, 0) = shift_x;
	shift(1, 0) = shift_y;

	Matrix2d trans_rot;
	double degreeToRad = 3.14159 / 180.0;
	trans_rot << cos((deg * degreeToRad)), -sin((deg * degreeToRad)), 
	sin((deg*degreeToRad)), cos((deg*degreeToRad));
	
	MatrixXd trans = OptimalPointMatcher::applyTransformation(bli, shift, trans_rot);
	//MatrixXd trans = (bli + shift);

	return OptimalPointMatcher::RMSE(trans, cbct);
}

/**
* 2D translation and rotation of one set of coordinates to another set of coordinates. 
*/
float f(Chromosome<float> c, MatrixXd &cbct, MatrixXd &bli) {
	float shift_x = c[0];
	float shift_y = c[1];
	float deg = c[2];
	MatrixXd ones = MatrixXd::Ones(21, 1);
	MatrixXd shift(2, 1);
	shift(0, 0) = shift_x;
	shift(1, 0) = shift_y;

	Matrix2d trans_rot;
	double degreeToRad = 3.14159 / 180.0;
	trans_rot << cos((deg * degreeToRad)), -sin((deg * degreeToRad)),
		sin((deg*degreeToRad)), cos((deg*degreeToRad));


	MatrixXd trans = OptimalPointMatcher::applyTransformation(bli, shift, trans_rot);
	//MatrixXd trans = (bli + shift);

	return OptimalPointMatcher::RMSE(trans, cbct);
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

/**
* Rosenbrock function. Defined as f(x, y) = (1-x)^2 + alpha * (y-x^2)^2.
* @param v vector with x y point.
* @return float the value of the function at the given point.
*/
float gaRosenbrock(Chromosome<float> c) {
	double x = c.at(0);
	double y = c.at(1);

	double alpha = 10.0;
	return (1 - x)*(1 - x) + alpha* ((y - (x*x)) *(y - (x*x)));
}

float ga1DTest(Chromosome<float> c) {
	double x = c.at(0);
	//parabola at 3, 5.
	return (pow((x - 3.0f), 2.0f)) + 5.0f;
}

float nm1DTest(Vector v) {
	double x = v[0];
	return (pow((x - 3.0f), 2.0f)) + 5.0f;
}

float gaBoothsFunctionTest(Chromosome<float> c) {
	float x = c.at(0);
	float y = c.at(1);
	return (pow((x + (2 * y) - 7), 2.0f) + pow(((2 * x) + y - 5), 2.0f));
}