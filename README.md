# Optimum
A small, lightweight library with various optimization algorithms including the Nelder-Mead algorithm. 

## What's in it? ##

* **Nelder-Mead Optimization Algorithm:** This algorithm is an optimization algorithm that uses a simplex to converge on a 'best' point. It takes n+1 points for an n-dimensional problem and creates a simplex that is then expanded, contracted, reflected and so on to move towards the minimum. This isn't very robust currently so it will get stuck in local minima. 

* **Genetic Algorithm** (Work in progress): This algorithm is optimization algorithm that introduces random mutations and cross overs (like nature and genetics) to make it more likely to find the true global minimum of some cost function. This is good to use when you know very little about the function you are minimizing or if it has a lot of local minima. 

* **Kabsch Algorithm:** This algorithm provides a way to find the optimal translation and rotation for the mapping of two point clouds to each other. In order to use this you need to know the mapping of your point pairs beforehand. 

* **Iterative Closest Point** (Work in progress): This algorithm provides a means to perform registration between two clouds of points in either 2D or 3D. This can be used when little is know about the mapping between the two point clouds, or if there is no matching points. 

## Why? ##
I made this library because these types of algorithms (for optimization/minimization) are typically only the topic of advanced research papers. I wanted to make this library as a spin-off of a github repository I found [online](https://github.com/blinry/nelder-mead-optimizer). I wanted to improve it and add to it, yet keep it easy to use and incorporate into any c++ project. Future support may include other languages if this does well. 

#### Usage ####

* Nelder-Mead Optimizer

````cpp
float precision = 0.0001;
int dimension = 2;
NelderMeadOptimizer o(dimension, precision);

// request a simplex to start with
Vector v(0.5, 0.5);
o.initialGuess(v);

while (!o.done()) {
    v = o.step(v, f(v));
}
````

* Genetic Algorithm

````cpp
float initialGuess[] = { 3.0f };
	Chromosome<float> chrome(initialGuess, 1);

	GeneticAlgorithm<float> ga(150, 300, 0.1f, 0.40f);
	Population<float> mPop = ga.generatePopulation(chrome);
	for (int i = 0; i < mPop.size(); i++) {
		Chromosome<float> member = mPop.getMember(i);
		//ga1DTest() returns the float values of a parabola with its vertex at (3, 5)
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
````

* Kabsch Algorithm

````cpp
    //pseudocode, load your data points. Could be a CSV file. 
    MatrixXd A = getAData();
    MatrixXd B = getBData();
    //create rotation to map B to A.
    MatrixXd rotation = OptimalPointMatcher::solveForOptimalRotation(A, B);
	MatrixXd trans = OptimalPointMatcher::solveForOptimalTranslation(A, B, rotation);
	//multiply trans by -1 here for corection of coordinate system direction. 
	MatrixXd B_trans = OptimalPointMatcher::applyTransformation(B,(trans * -1.0), rotation);
	//calculate root mean squre error between new and reference points. 
	double rmse = OptimalPointMatcher::RMSE(B_trans, A);
````

#### License ####
-------
This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

See LICENSE file for a copy of the GNU General Public License.

Copyright (C) 2015  Paul T <developer.paul.123@gmail.com>
