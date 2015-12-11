#include "neldermead.h"

/**
* Constructor for Nelder-Mead Optimizer.
* @param dimension the number of variables we are looking to optimize.
* @param termination_distance the minimum distance that is acceptable as a termination
*			criteria. Default is 0.001.
* @param maxIterations the maximum iterations to run this algorithm before it terminates.
*			Default is 1,000,000.
*/
NelderMeadMinimizer::NelderMeadMinimizer(int dimension, float termination_distance, int maxIterations) {
	this->dimension = dimension;
	srand(time(NULL));
	//adaptive params.
	if (dimension > 0) {
		alpha = 1.0;
		gamma = 1.0 + (2.0 / dimension);
		rho = 0.75 - (1.0 / (2.0 * dimension));
		sigma = 1.0 - (1.0 / dimension);
	}
	else {
		alpha = 1;
		gamma = 2;
		rho = 0.5;
		sigma = 0.5;
	}
	this->maxIter = maxIterations;
	this->iterations = 0;
	this->termination_distance = termination_distance;
}

// used in `step` to sort the vectors
bool NelderMeadMinimizer::operator()(const Vector& a, const Vector& b) {
	return db.lookup(a) < db.lookup(b);
}

/**
* Set the initial guess of the absolute minimum.
* @param v the initial guess.
*/
void NelderMeadMinimizer::initialGuess(Vector v) {
	vectors.push_back(v);
	for (int n = 0; n < dimension; n++) {
		Vector newVec;
		newVec.prepare(dimension);
		for (int i = 0; i < dimension; i++) {
			float tau = 0;
			if (v[i] == 0) {
				tau = 0.00025;
			}
			else {
				tau = 0.05;
			}
			Vector unitVec = v.unitVector();
			float val = v[i] + unitVec[i] * tau;
			newVec.set(val, i);
		}
		vectors.push_back(newVec);
	}

	//set up the extra point xs.  
	xs.prepare(dimension);
	for (int i = 0; i < dimension; i++) {
		Vector v = vectors[i];
		float val = v[i];
		xs.set(val, i);
	}

}

/**
* Termination criteria: each pair of vectors in the simplex has to
* have a distance of at most `termination_distance`
* @return bool true if the criteria has been met.
*/
bool NelderMeadMinimizer::done() {
	if (vectors.size() < dimension) {
		return false;
	}
	if (iterations >= maxIter) {
		return true;
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

/**
* Gets the number of iterations the minimizer has gone through. 
* @return int the number of interations. 
*/
int NelderMeadMinimizer::getIterations() {
	return iterations;
}

Vector NelderMeadMinimizer::getExtraVector() {
	return xs;
}

void NelderMeadMinimizer::setExtraVectorScore(float score) {
	db.insert(xs, score);
}

/**
* Step the algorithm by running through the Nelder-Mead optimization loop.
* @param vec the vector to test.
* @param score the score from the objective function of the current vector.
*/
Vector NelderMeadMinimizer::step(Vector vec, float score) {
	//save vector and score. 
	db.insert(vec, score);
	iterations++;
	try {
		if (vectors.size() < dimension + 1) {
			vectors.push_back(vec);
		}

		// otherwise: optimize!
		if (vectors.size() == dimension + 1) {
			while (!done()) {
				//sorts in ascending order. So first one will be the lowest. 
				sort(vectors.begin(), vectors.end(), *this);
				gradients.clear();

				Vector cog; // center of gravity
				cog.prepare(dimension);
				//calculate centroid for all points except Xn+1
				for (int i = 0; i < dimension; i++) {
					cog += vectors[i];
				}
				cog /= dimension;

				Vector worst = vectors[dimension];
				Vector best = vectors[0];
				Vector second_worst = vectors[1];

				//Testing out gradient method. 		
			/*	float funcDiff = f(best) - f(xs);
				float pointDiff = second_worst[0] - best[0];
				float gX = funcDiff / pointDiff;

				float funcDiffTwo = f(second_worst) - f(xs);
				float pointDiffTwo = best[1] - second_worst[1];
				float gY = funcDiffTwo / pointDiffTwo;*/

				//calculate quasi gradients.
				//Vector gradient(gX, gY);

				// reflect through the COG. 
				Vector reflected = cog + (cog - worst)*alpha;
				//Vector reflected = best - (gradient)*1.0;

				//if the reflected point is better (lower score) than the second worst but
				//has a higher score than the best...replace the worst point with the reflected point. 
				if (f(best) <= f(reflected) && f(reflected) < f(second_worst)) {
					vectors[dimension] = reflected;
				}

				//if the reflected point is the best point so far the compute the 
				//expanded point. 
				else if (f(reflected) < f(best)) {
					// expand
					Vector expanded = cog + (reflected - cog)*gamma;

					//trying out gradient method. 
					//Vector expanded = (best*(1 - gamma)) + (reflected)* gamma;

					//if the expanded point is better than the reflected point 
					//then replace the worst point with the expanded point.
					if (f(expanded) < f(reflected)) {
						vectors[dimension] = expanded;
					}
					//otherwise replace the worst point with the reflected point. 
					else {
						vectors[dimension] = reflected;
					}
				}

				// outside contraction
				else if (f(second_worst) <= f(reflected) && f(reflected) < f(worst)) {
					Vector outContracted = cog + (reflected - cog)*rho;
					if (f(outContracted) <= f(reflected)) {
						//replace worst with outside contraction. 
						vectors[dimension] = outContracted;
					}
					//reduction. For all but the best point replace the point with
					//Xi = X1 + sigma* (Xi - X1);
					else {
						for (int i = 1; i <= dimension; i++) {
							vectors[i] = best + (vectors[i] - best)*sigma;
						}
					}
				}
				//inside contraction. 
				else if (f(reflected) >= f(worst)) {
					Vector inContracted = cog - (reflected - cog)*rho;
					if (f(inContracted) < f(worst)) {
						//replace worst with in contracted. 
						vectors[dimension] = inContracted;
					}
					//reduction. For all but the best point replace the point with
					//Xi = X1 + sigma* (Xi - X1);
					else {
						for (int i = 1; i <= dimension; i++) {
							vectors[i] = best + (vectors[i] - best)*sigma;
						}
					}
				}

				//reduction. For all but the best point replace the point with
				//Xi = X1 + sigma* (Xi - X1);
				else {
					for (int i = 1; i <= dimension; i++) {
						vectors[i] = best + (vectors[i] - best)*sigma;
					}
				}
			}

			// algorithm is terminating, output: simplex' center of gravity
			Vector cog;
			cog.prepare(dimension);
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
				result[i] = 0.001*(rand() % 10000);
			}
			return result;
		}
	}
	catch (Vector v) {
		return v;
	}
}

/**
* Insert a vector into the algorithm.
* @param vec the vector to insert.
*/
void NelderMeadMinimizer::insert(Vector vec) {
	if (vectors.size() < dimension + 1) {
		vectors.push_back(vec);
	}
}

/**
* Look up the score of the given vector.
* @param vec the vector.
* @return float the score of the vector.
*/
float NelderMeadMinimizer::f(Vector vec) {
	return db.lookup(vec);
}