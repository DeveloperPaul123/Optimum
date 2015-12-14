/*
This is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

See LICENSE file for a copy of the GNU General Public License.

Copyright (C) 2015 Paul Tsouchlos

Taken from https://github.com/blinry/nelder-mead-optimizer.

This code has been adapted from the original code by
Sebastian Morr sebastian@morr.cc (Copyright 2013)

Usage:
<p>
float precision = 0.001;
int dimension = 2;
NelderMeadOptimizer o(dimension, precision);

// request a simplex to start with
Vector v(0.5, 0.5);
o.insert(v);
o.insert(Vector(0.1, 0.1));
o.insert(Vector(0.2, 0.7));

while (!o.done()) {
v = o.step(v, f(v));
}

//minimizer

NelderMeadMinimizer r(dimension, precision);

Vector q(2, 1);
r.initialGuess(q);
while(!r.done()) {
float score = f(q);
q = r.step(q, score);
}

</p>
*/

#ifndef NELDERMEAD_H
#define NELDERMEAD_H

#include <ctime>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
using namespace std;

/**
* Float vector with standard operations.
*/
class Vector {
public:

	/**
	* Default constructor.
	*/
	Vector() {
	}

	/**
	* One dimensional vector. 
	*/
	Vector(float c0) {
		coords.push_back(c0);
	}

	/**
	* Two dimensional vector.
	*/
	Vector(float c0, float c1) {
		coords.push_back(c0);
		coords.push_back(c1);
	}

	/**
	* Three dimensional vector.
	*/
	Vector(float c0, float c1, float c2) {
		coords.push_back(c0);
		coords.push_back(c1);
		coords.push_back(c2);
	}

	/**
	* Multidimensional vector constructor, for any size. 
	* @param *values array of float values.
	* @param size the length of the array of values. 
	*/
	Vector(float *values, int size) {
		for (int i = 0; i < size; i++) {
			coords.push_back(values[i]);
		}
	}

	/**
	* [] Operator to make value access easy. 
	* @param i the index position to look up. 
	* @return float the value at the index position. 
	*/
	float& operator[](int i) {
		return coords[i];
	}

	/**
	* Get the value at a given position.
	* @param i the position to lookup.
	* @return float the value.
	*/
	float at(int i) const {
		return coords[i];
	}

	/**
	* Get the dimension of this vector.
	* @return int the dimension of the vector.
	*/
	int dimension() const {
		return coords.size();
	}

	/**
	* Prepare the vector with 0 values.
	* @param size the dimension of this vector.
	*/
	void prepare(int size) {
		for (int i = 0; i<size; i++) {
			coords.push_back(0);
		}
	}

	/**
	* Addition operator for vectors. 
	* @param other the other vector. 
	* @return Vector the result of the addition. 
	*/
	Vector operator+(Vector other) {
		Vector result;
		result.prepare(dimension());
		for (int i = 0; i<dimension(); i++) {
			result[i] = coords[i] + other[i];
		}
		return result;
	}

	/**
	* Plus equals operator for vectors. 
	* @param the other vector. 
	*/
	void operator+=(Vector other) {
		for (int i = 0; i<dimension(); i++) {
			coords[i] += other[i];
		}
	}

	/**
	* Minus operator for vectors. 
	* @param other the other vector. 
	* @return Vector the result of the subtraction. 
	*/
	Vector operator-(Vector other) {
		Vector result;
		result.prepare(dimension());
		for (int i = 0; i<dimension(); i++) {
			result[i] = coords[i] - other[i];
		}
		return result;
	}

	/**
	* Equals (bool) operator. 
	* @param other the other vector to compare to. 
	* @return bool true if the vectors are the same, false otherwise. 
	*/
	bool operator==(Vector other) {
		if (dimension() != other.dimension()) {
			return false;
		}
		for (int i = 0; i<dimension(); i++) {
			if (other[i] != coords[i]) {
				return false;
			}
		}
		return true;
	}

	/**
	* Multiplication operator. 
	* @param factor the factor to multiply this vector by. 
	* @return Vector the new vector. 
	*/
	Vector operator*(float factor) {
		Vector result;
		result.prepare(dimension());
		for (int i = 0; i<dimension(); i++) {
			result[i] = coords[i] * factor;
		}
		return result;
	}

	/**
	* Division operator for vectors. 
	* @param factor the number to divide this vector by. 
	* @return the newly divided vector. 
	*/
	Vector operator/(float factor) {
		Vector result;
		result.prepare(dimension());
		for (int i = 0; i<dimension(); i++) {
			result[i] = coords[i] / factor;
		}
		return result;
	}

	/**
	* Divide equals operator. 
	* @param factor the number to divide this vector by. 
	*/
	void operator/=(float factor) {
		for (int i = 0; i<dimension(); i++) {
			coords[i] /= factor;
		}
	}

	/**
	* Less than operator. 
	* @param other the other vector to compare to. 
	* @return bool true is this vector is less than the other. False otherwise. 
	*/
	bool operator<(const Vector other) const {
		for (int i = 0; i<dimension(); i++) {
			if (at(i) < other.at(i))
				return false;
			else if (at(i) > other.at(i))
				return true;
		}
		return false;
	}

	/**
	* Returns the length of this vector.
	*/
	float length() {
		float sum = 0;
		for (int i = 0; i<dimension(); i++) {
			sum += coords[i] * coords[i];
		}
		return pow(sum, 0.5f);
	}

	/**
	* Set the value of a point in this vector.
	* @param value the new value.
	* @param position the position to set the point at.
	*/
	void set(float value, int position) {
		coords[position] = value;
	}

	/**
	* Calculates the unit vector of this vector.
	* @return Vector the unit vector of this vector. 
	*/
	Vector unitVector() {
		float len = length();
		Vector unitVec;
		unitVec.prepare(dimension());
		for (int i = 0; i < dimension(); i++) {
			float val = 0;
			if (len == 0) {
				val = 1;
			}
			else {
				val = coords[i] / len;
			}
			unitVec.set(val, i);
		}
		return unitVec;
	}

private:
	vector<float> coords;
};

/**
* This class stores known values for vectors.
* It throws unknown vectors.
*/
class ValueDB {
public:
	/**
	* Empty constructor. 
	*/
	ValueDB() {
	}

	/**
	* Look up the stored value for a given vector. 
	* If the vector isn't contained in this DB it is thrown. 
	* @param vec the Vector to look up. 
	* @return float the score associated with the given vector if
	*		it exists in the database. 
	*/
	float lookup(Vector vec) {
		if (!contains(vec)) {
			throw vec;
		}
		else {
			return values[vec];
		}
	}

	/**
	* Inserts a vector into the database with a given value. 
	* @param vec the vector to insert in the database. 
	* @param value the score of the given vector. 
	*/
	void insert(Vector vec, float value) {
		values[vec] = value;
	}
private:

	/**
	* Checks to see if a given vector is contained in this database. 
	* @param vec the vector to check. 
	* @return bool true if the vector exists in the DB, false otherwise. 
	*/
	bool contains(Vector vec) {
		map<Vector, float>::iterator it = values.find(vec); // TODO add tolerance
		return it != values.end();
	}
	map<Vector, float> values;
};

/**
* NelderMeadMinimizer class. Uses the Nelder-Mead method for minimization
* of a given cost function.
*/
class NelderMeadMinimizer {
public:

	NelderMeadMinimizer(int dimension, float termination_distance = 0.001, int maxIterations = 1000000);
	// used in `step` to sort the vectors
	bool operator()(const Vector& a, const Vector& b);
	void initialGuess(Vector v);
	bool done();
	int getIterations();
	Vector step(Vector vec, float score);
	Vector getExtraVector();
	void setExtraVectorScore(float score);
private:

	void insert(Vector vec);
	float f(Vector vec);
	Vector xs;
	int dimension, iterations, maxIter;
	double alpha, gamma, rho, sigma;
	float termination_distance;
	vector<Vector> vectors;
	vector<Vector> gradients;
	ValueDB db;
};


/**
* Nelder Mead maximizer.
* This is the original code from https://github.com/blinry/nelder-mead-optimizer
*/
class NelderMeadOptimizer {
public:
	NelderMeadOptimizer(int dimension, float termination_distance = 0.001) {
		this->dimension = dimension;
		srand(time(NULL));
		alpha = 1;
		gamma = 2;
		rho = -0.5;
		sigma = 0.5;
		this->termination_distance = termination_distance;
	}
	// used in `step` to sort the vectors
	bool operator()(const Vector& a, const Vector& b) {
		return db.lookup(a) < db.lookup(b);
	}
	// termination criteria: each pair of vectors in the simplex has to
	// have a distance of at most `termination_distance`
	bool done() {
		if (vectors.size() < dimension) {
			return false;
		}
		for (int i = 0; i<dimension + 1; i++) {
			for (int j = 0; j<dimension + 1; j++) {
				if (i == j) continue;
				if ((vectors[i] - vectors[j]).length() > termination_distance) {
					return false;
				}
			}
		}
		return true;
	}
	void insert(Vector vec) {
		if (vectors.size() < dimension + 1) {
			vectors.push_back(vec);
		}
	}
	Vector step(Vector vec, float score) {
		db.insert(vec, score);
		try {
			if (vectors.size() < dimension + 1) {
				vectors.push_back(vec);
			}

			// otherwise: optimize!
			if (vectors.size() == dimension + 1) {
				while (!done()) {
					sort(vectors.begin(), vectors.end(), *this);
					Vector cog; // center of gravity
					cog.prepare(dimension);
					for (int i = 1; i <= dimension; i++) {
						cog += vectors[i];
					}
					cog /= dimension;
					Vector best = vectors[dimension];
					Vector worst = vectors[0];
					Vector second_worst = vectors[1];
					// reflect
					Vector reflected = cog + (cog - worst)*alpha;
					if (f(reflected) > f(second_worst) && f(reflected) < f(best)) {
						vectors[0] = reflected;
					}
					else if (f(reflected) > f(best)) {
						// expand
						Vector expanded = cog + (cog - worst)*gamma;
						if (f(expanded) > f(reflected)) {
							vectors[0] = expanded;
						}
						else {
							vectors[0] = reflected;
						}
					}
					else {
						// contract
						Vector contracted = cog + (cog - worst)*rho;
						if (f(contracted) > f(worst)) {
							vectors[0] = contracted;
						}
						else {
							for (int i = 0; i<dimension; i++) {
								vectors[i] = best + (vectors[i] - best)*sigma;
							}
						}
					}
				}

				// algorithm is terminating, output: simplex' center of gravity
				Vector cog;
				for (int i = 0; i <= dimension; i++) {
					cog += vectors[i];
				}
				return cog / (dimension + 1);
			}
			else {
				// as long as we don't have enough vectors, request random ones,
				// with coordinates between 0 and 1. If you want other start vectors,
				// simply ignore these and use `step` on the vectors you want.
				Vector result;
				result.prepare(dimension);
				for (int i = 0; i<dimension; ++i) {
					result[i] = 0.001*(rand() % 1000);
				}
				return result;
			}
		}
		catch (Vector v) {
			return v;
		}
	}
private:
	float f(Vector vec) {
		return db.lookup(vec);
	}
	int dimension;
	float alpha, gamma, rho, sigma;
	float termination_distance;
	vector<Vector> vectors;
	ValueDB db;
};

#endif NELDERMEAD_H