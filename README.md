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

* **Note:** All classes are under the `Optimum` namespace. This is to avoid conflicts with classes from other namespaces like std or open cv (which I use frequently). 

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

#### Contributing ####
-------
If you would like to contribute, feel free to fork this repository and create a pull request. This project is configured with Cmake 3.3 or higher and depends on [Eigen] (https://bitbucket.org/eigen/eigen/) which is a header only, linear algebra C++ library. 

#### License ####
-------
Copyright 2016 Paul T

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Copyright (C) 2015  Paul T <developer.paul.123@gmail.com>
