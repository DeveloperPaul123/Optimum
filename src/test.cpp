#include <iostream>
#include <math.h>
#include <algorithm>
#include "gtest/gtest.h"

#include "neldermead.h"

class NelderMeadTest : public ::testing::Test {
protected:
	NelderMeadTest() {
	}
	
	virtual ~NelderMeadTest() {
	}
	
	virtual void SetUp() {
	}
	
	virtual void TearDown() {
	}
	
	Optimum::NelderMeadMinimizer nmm;
	
};

TEST_F(NelderMeadTest, MethodBarDoesAbc) {

	float precision = 0.000001;
	int dimension = 1;

	nmm = Optimum::NelderMeadMinimizer(dimension, precision);
	Optimum::Vector v(527.0f);
	nmm.initialGuess(v);
	while (!nmm.done()) {
		double x = v[0];
		float score = x - 3.0;
		v = nmm.step(v, score);
	}
	double best = v.at(0);
	EXPECT_NEAR(3.0, v.at(0), 0.1);
}

