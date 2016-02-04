#include <iostream>
#include <math.h>
#include <algorithm>
#include "gtest/gtest.h"
#include "utils.h"
#include "neldermead.h"
#include "icp.h"
#include <Eigen\Dense>
#include <Eigen\Geometry>


/**
* Test fixture for NelderMead
*/
class NelderMeadTest : public ::testing::Test {
protected:
	NelderMeadTest() {
	}

	~NelderMeadTest() {
	}

	void SetUp() {

	}

	void TearDown() {
	}
	Optimum::NelderMeadMinimizer *nmm;
};

/**
* Test fixture for ICP tester. 
*/
class ICPTester : public ::testing::Test {
protected:
	ICPTester() {

	}

	~ICPTester() {

	}

	void SetUp() {

	}

	void TearDown() {

	}
	Optimum::ICP *icp;
};

TEST_F(NelderMeadTest, testPowerFunction) {
	float precision = 0.01;
	int dimension = 1;
	nmm = new Optimum::NelderMeadMinimizer(dimension, precision);
	Optimum::Vector v(5.0f);
	nmm->initialGuess(v);
	while (!nmm->done()) {
		double x = v[0];
		float score = (pow((x - 3.0f), 2.0f)) + 5.0f;
		v = nmm->step(v, score);
	}
	double best = v.at(0);
	EXPECT_NEAR(3.00, best, 0.1);
}

TEST_F(NelderMeadTest, testLinearFunction) {
	float precision = 0.01;
	int dimension = 1;
	nmm = new Optimum::NelderMeadMinimizer(dimension, precision);
	Optimum::Vector v(2.0f);
	nmm->initialGuess(v);
	while (!nmm->done()) {
		double x = v[0];
		float score = 2 * x + 5;
		v = nmm->step(v, score);
	}
	double best = v.at(0);
	EXPECT_NEAR(-2.5, best, 0.1);
}

/**
* Test transform from older test. 
*/
TEST_F(ICPTester, testTransform) {
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

	Eigen::MatrixXd A(21, 2);
	Eigen::MatrixXd B(21, 2);

	for (int r = 0; r < 21; r++) {
		A(r, 0) = arrA[r][0];
		A(r, 1) = arrA[r][1];
		double x = arrB[r][0];
		double y = arrB[r][1];

		double correctedX = (x - (0.5 * 2048))*0.0971;
		double correctedY = ((0.5 * 2048) - y) * 0.0971;
		B(r, 0) = correctedX;
		B(r, 1) = correctedY;
	}

	//mapping B to A.
	Optimum::PointMatcherSettings settings;
	Optimum::PointMatcher *matcher = new Optimum::PointMatcher(settings);
	Optimum::Transform t = matcher->solve(B, A);

	Eigen::MatrixXd transformed = Optimum::PointMatcher::applyTransformation(B, t.translation, t.rotation);
	double error = Optimum::PointMatcher::RMSE(transformed, A);

	EXPECT_NEAR(4.0, error, 1.5);
}

/**
* Test transform from most recent system test. 
*/
TEST_F(ICPTester, testTransformTwo) {
	double arrBLI[12][3] = {
		{ 477.409, 359.402, 0 },
		{ 669.795, 353.769, 0 },
		{ 307.372, 537.054, 0 },
		{ 481.644, 532.634, 0 },
		{ 672.669, 528.672, 0 },
		{ 847.806, 523.96, 0 },
		{ 310.347, 713.711, 0 },
		{ 483.897, 710.377, 0 },
		{ 675.428, 703.24, 0 },
		{ 848.932, 700.778, 0 },
		{489.53, 879.915, 0 },
		{678.346, 876.372, 0 }
	};

	double arrCBCT[12][3] = {
		{ -12.887, 28.2005, -20.875 },
		{ 9.14438, 28.3098, -20.875 },
		{ -33.0168, 8.31108, -20.875 },
		{ -13.0182, 8.18, -20.875 },
		{ 9.05681, 8.00511, -21.125 },
		{28.99, 8.09249, -21.125},
		{ -33.1262, -11.7532, -21.125 },
		{ -13.1711, -11.95, -21.125 },
		{ 9.03499, -12.0373, -21.125 },
		{ 28.8152, -11.9936, -21.125 },
		{ -13.2052, -31.9509, -19.875 },
		{ 8.93298, -31.9168, -20.125 },
	};

	Eigen::MatrixXd cbct_raw(12, 3);
	Eigen::MatrixXd bli_raw(12, 3);

	for (int r = 0; r < cbct_raw.rows(); r++){
		for (int c = 0; c < cbct_raw.cols(); c++) {
			cbct_raw(r, c) = arrCBCT[r][c];
			bli_raw(r, c) = arrBLI[r][c];
		}
	}

	Eigen::MatrixXd cbct = cbct_raw.middleCols(0, 2);
	Eigen::MatrixXd bli = bli_raw.middleCols(0, 2);

	//correct the data. 
	for (int r = 0; r < bli.rows(); r++) {
		double x = bli(r, 0);
		double y = bli(r, 1);
		double correctedX = (x - (0.5 * 1024))*0.1136;
		double correctedY = ((0.5 * 1024) - y) * 0.1136;
		bli(r, 0) = correctedX;
		bli(r, 1) = correctedY;
	}

	Optimum::PointMatcherSettings pSet;
	//perform scaling.
	pSet.doScale = true;
	Optimum::PointMatcher *matcher = new Optimum::PointMatcher(pSet);
	Optimum::Transform transform = matcher->solve(bli, cbct);

	Eigen::MatrixXd transformed = Optimum::PointMatcher::applyTransformation(bli, 
		transform.translation, transform.rotation);
	double rmse = Optimum::PointMatcher::RMSE(transformed, cbct);

	EXPECT_NEAR(rmse, 4.38, 1.0);
}

TEST_F(ICPTester, testICP) {
	Eigen::MatrixXd ref(2, 2);
	Eigen::MatrixXd tar(2, 2);

	ref(0, 0) = 1.0;
	ref(0, 1) = 1.0;
	ref(1, 0) = 5.0;
	ref(1, 1) = 5.0;

	tar(0, 0) = 2.5;
	tar(0, 1) = 2.5;
	tar(1, 0) = 6.5;
	tar(1, 1) = 6.5;

	Optimum::ICPSettings settings;
	settings.maxIterations = 20;
	settings.pointType = Optimum::TWO_D;
	icp = new Optimum::ICP(ref, tar, settings);
	icp->solve();

	Optimum::PointMatcherSettings pset;
	Optimum::PointMatcher *matcher = new Optimum::PointMatcher(pset);
	Optimum::Transform t = matcher->solve(tar, ref);

	std::cout << "Best " << icp->getBestTranslation() << " Optimal: " << t.translation << std::endl;
	std::cout << "Best Rotation " << icp->getBestRotation() << " Optimal Rotation: " << t.rotation << std::endl;
	EXPECT_EQ(t.translation, icp->getBestTranslation());

}
