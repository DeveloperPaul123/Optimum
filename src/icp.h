#ifndef ICP_H
#define ICP_H

#include <iostream>
#include <algorithm>
#include <math.h>
#include <Eigen\Dense>
#include <Eigen\Jacobi>

#define PI 3.14159265

using namespace std;
using namespace Eigen;

/**
* Simple point class for ICP algorithm. 
*/
class Point {
	vector<double> coords;
public:

	/**
	* Empty constructor. 
	*/
	Point() {

	}
	
	/**
	* Constructor for a 2D point. 
	* @param x the x coordinate. 
	* @param y the y coordinate. 
	*/
	Point(double x, double y) {
		coords.push_back(x);
		coords.push_back(y);
	}

	/**
	* Constructor for a 3D point.
	* @param x the x coordinate. 
	* @param y the y coordinate. 
	* @param z the z coordinate. 
	*/
	Point(double x, double y, double z) {
		coords.push_back(x);
		coords.push_back(y);
		coords.push_back(z);
	}
	
	/**
	* Prepare the data for a certain number of values. Use
	* this when you initialize a point with the empty constructor
	* before setting data inside the point. 
	* @param size the size this point is to be. 
	*/
	void prepare(int size) {
		for (int i = 0; i < size; i++) {
			coords.push_back(0.0);
		}
	}

	/**
	* Set a value at a given position. 
	* @param value the new value. 
	* @param pos the position to set the new value at. 
	*/
	void setValue(double value, int pos) {
		coords[pos] = value;
	}

	/**
	* Replaces the points in this point up to but not including 'size'
	* @param *points an array of new data. 
	* @param size the size of the new data array. 
	*/
	void setPoints(double * points, int size) {
		for (int i = 0; i < size; i++) {
			coords.push_back(points[i]);
		}
	}
	
	/**
	* Returns the size of this point. 
	* @return int the size of the point. 
	*/
	int size() {
		return coords.size();
	}

	/**
	* Convinience operator for accessing data in this point. 
	* @param i the index to look up. 
	* @return double the data at index point i. 
	*/
	double& operator[] (const int &i) {
		return coords[i];
	}

	/**
	* Minus operator for a point. Note: performs *this - other. 
	* @param other the other point to subtract. 
	* @return Point the result of subtracting other from this point. 
	*/
	Point operator-(Point other) {
		Point p;
		p.prepare(other.size());
		for (int i = 0; i < other.size(); i++) {
			p[i] = coords[i] - other[i];
		}
		return p;
	}

	/**
	* Calculates the Euclidean distance to a given point from this point. 
	* @param other the other point. 
	* @return double the distance. 
	*/
	double distanceTo(Point &other) {
		Point p = Point(*this);
		Point diff = p - other;
		vector<double> sqs;
		for (int i = 0; i < other.size(); i++) {
			double d = diff[i];
			double sq = d*d;
			sqs.push_back(sq);
		}

		double sum = 0.0;
		for (int j = 0; j < sqs.size(); j++) {
			sum += sqs[j];
		}

		return sqrt(sum);
	}

};

/**
* Helper class with static functions to calculate optimal rotation and translation using
* the Kabsch method.
*/
class OptimalPointMatcher {

public:

	/**
	* Solves for the optimal translation between the reference and target points 
	* given an optimal rotation matrix. This is trying to map the reference to the target. 
	* @param ref the reference points.
	* @param target the target points.
	* @optimalRotation the optimal rotation matrix calculated by the function
	*				solveForOptimalRotation.
	*/
	static MatrixXd solveForOptimalTranslation(MatrixXd ref, MatrixXd target, MatrixXd optimalRotation) {
		MatrixXd centroidOne = getCentroid(target);
		MatrixXd centroidTwo = getCentroid(ref);
		MatrixXd translation = ((optimalRotation*-1.0)*centroidOne) + centroidTwo;
		return translation;
	}

	/**
	* Solves for the optimal rotation matrix using singular value decomposition. This is
	* trying to map the reference to the target. 
	* @param ref the reference points
	* @param target the target points.
	* @return MatrixXd the optimal rotation matrix.
	*/
	static MatrixXd solveForOptimalRotation(MatrixXd ref, MatrixXd target) {

		MatrixXd centroidOne = getCentroid(target);
		MatrixXd centroidTwo = getCentroid(ref);

		//center points at centroid
		int rows = ref.rows();
		for (int r = 0; r < rows; r++) {
			ref.row(r) = ref.row(r) - centroidTwo.transpose();
			target.row(r) = target.row(r) - centroidOne.transpose();
		}

		//calculate covarience matrix.
		MatrixXd cov = target.transpose() * ref;

		//calculate svd
		JacobiSVD<MatrixXd> svd(cov, ComputeFullU | ComputeFullV);

		//if svd(H) = [U, S, V], then rotation R = V*U(transpose)
		MatrixXd rotation = svd.matrixV() * svd.matrixU().transpose();

		//check for correction. 
		if (rotation.determinant() < 0) {
			rotation.row(rotation.rows()-1) *= -1;
		}

		return rotation;
	}

	/**
	* Calculates the root mean square error between two matricies.
	* @param dataOne the first matrix in the comparison.
	* @param dataTwo the second matrix in the comparison.
	* @return double the root mean square error or 0 if the matrix sizes don't match.
	*/
	static double RMSE(MatrixXd dataOne, MatrixXd dataTwo) {
		if (dataOne.size() != dataTwo.size()) {
			cout << "Error, matrix size mismatch." << endl;
			return 0.0f;
		}
		else {
			MatrixXd result = dataOne - dataTwo;
			for (int i = 0; i < result.rows(); i++) {
				for (int j = 0; j < result.cols(); j++) {
					double val = result(i, j);
					double newVal = val*val;
					result(i, j) = newVal;
				}
			}

			vector<double> sqrts;
			for (int r = 0; r < result.rows(); r++) {
				double sumSq = result.row(r).sum();
				sqrts.push_back(sqrt(sumSq));
			}

			double rmse = 0.0;
			for (int x = 0; x < sqrts.size(); x++) {
				rmse += sqrts[x];
			}
			return rmse;
		}
	}
	/**
	* Applies a transformation to a given array of data with a translation and a rotation.
	* @param data the original data.
	* @param translation the translation matrix, typically a row or column vector.
	* @param rotation a rotation matrix.
	* @return MatrixXd the transformed data set of data.
	*/
	static MatrixXd applyTransformation(MatrixXd data, MatrixXd translation, MatrixXd rotation) {
		MatrixXd trans(data.rows(), data.cols());
		MatrixXd ones = MatrixXd::Ones(data.rows(), 1);
		if (translation.rows() != data.rows() || translation.cols() != data.cols()) {
			if (translation.cols() == 1) {
				//column vector. 
				for (int i = 0; i < data.cols(); i++) {
					trans.col(i) = translation(i, 0) * ones;
				}
			}
			else if (translation.rows() == 1) {
				//row vector. 
				for (int i = 0; i < data.cols(); i++) {
					trans.col(i) = translation(0, i) * ones;
				}
			}
		}
		return (data * rotation) + trans;
	}

	/**
	* Returns the centroid of a given matrix. The result will be a column vector with the average of each
	* column in the data corresponding to each row of the column vector. So if your input is of the size n X i, then
	* this function will return a column vector of size i x 1 with the number of each row equal to the mean of the ith column in your
	* data.
	* @param points the data.
	* @return MatrixXd a column vector with the centroid points.
	*/
	static MatrixXd getCentroid(MatrixXd points) {

		int columns = points.cols();
		int rows = points.rows();
		//create a column vector. 
		MatrixXd centroid(columns, 1);
		for (int i = 0; i < columns; i++) {
			double sum = points.col(i).sum();
			double cent = sum / rows;
			centroid(i, 0) = cent;
		}
		return centroid;
	}

};

/**
* Enumeration for if we're in 3D or 2D.
*/
enum PointType { TWO_D = 0, THREE_D };

/**
* Structure used to setup parameters for the ICP algorithm. 
*/
struct ICPSettings{
	PointType pointType;
	int maxIterations;
};

/**
* Iterative closest point algorithm. 
*/
class ICP {

	MatrixXd mRef;
	MatrixXd mTarget;
	MatrixXd curRef;
	MatrixXd solvedTransform, solvedRotation;
	MatrixXd transform, lastTransform;
	MatrixXd rotation, lastRotation;
	vector<Point> refPoints;
	vector<Point> tarPoints;
	vector<Point> closestPoints;
	int iterations;

public:
	/**
	* Iterative closest point algorithm. It is assumed that a row in both the target and 
	* reference data is a point. 
	* @param reference the ref data. 
	* @param target the target points. 
	* @param set ICPSettings struct. 
	*/
	ICP(MatrixXd reference, MatrixXd target, ICPSettings set) {
		this->mRef = reference;
		this->mTarget = target;
		this->iterations = set.maxIterations;
		this->curRef = reference;
	}

	void solve() {
		
		rotation = MatrixXd::Identity(mRef.cols(), mRef.cols());
		transform = MatrixXd::Zero(mRef.cols(), 1);

		cout << rotation << endl; cout << endl;
		cout << transform << endl; cout << endl;

		while (iterations > 0) {
			
			closestPoints.clear();
			tarPoints.clear();
			refPoints.clear();
		
			lastTransform = transform;
			lastRotation = rotation;

			MatrixXd closest = match();

			//get rotation and translation. 
			MatrixXd newRot = OptimalPointMatcher::solveForOptimalRotation(closest, mTarget);
			MatrixXd newTrans = OptimalPointMatcher::solveForOptimalTranslation(closest, mTarget, newRot);
			
			double lastAngle = rotMatrixToDegrees(lastRotation);
			double angle = rotMatrixToDegrees(newRot);
			double newAngle = lastAngle + angle;

			rotation = angleToRotMatrix(newAngle, lastRotation.rows());

			transform = newTrans;
			
			cout << rotation << endl; cout << endl;
			cout << transform << endl; cout << endl;
			
			//update the target. 
			MatrixXd newRef = OptimalPointMatcher::applyTransformation(curRef, transform, rotation);
			double error = OptimalPointMatcher::RMSE(mTarget, newRef);

			cout << "Error: " << error << endl;

			curRef = newRef;

			//work towards ending of loop. 
			iterations -= 1;
		}
	}

	MatrixXd getBestTranslation() {
		return transform;
	}

	MatrixXd getBestRotation() {
		return rotation;
	}

private:
	MatrixXd match() {
		//set the points. 
		int rows = mRef.rows();
		int cols = mRef.cols();
		for (int r = 0; r < rows; r++) {
			Point ref;
			Point tar;
			ref.prepare(cols);
			tar.prepare(cols);
			for (int c = 0; c < cols; c++) {
				ref.setValue(curRef(r, c), c);
				tar.setValue(mTarget(r, c), c);
			}
			tarPoints.push_back(tar);
			refPoints.push_back(ref);
		}

		//get the closest points between target and reference. 
		for (int i = 0; i < refPoints.size(); i++) {
			double bestDistance = 1000000.0;
			Point tar = tarPoints.at(i);
			Point best;
			best.prepare(refPoints.at(0).size());
			//linear search. 
			for (int j = 0; j < tarPoints.size(); j++) {
				Point ref = refPoints.at(j);
				double dist = tar.distanceTo(ref);
				if (dist < bestDistance) {
					bestDistance = dist;
					best = ref;
				}
			}
			closestPoints.push_back(best);
		}

		//create matrix of closest points. 
		MatrixXd closest(closestPoints.size(), closestPoints.at(0).size());
		for (int x = 0; x < closestPoints.size(); x++) {
			Point p = closestPoints.at(x);
			for (int y = 0; y < p.size(); y++) {
				closest(x, y) = p[y];
			}
		}

		return closest;
	}

	double rotMatrixToDegrees(MatrixXd rot) {
		return asin(rot(1, 0)) * 180.0 / PI;
	}

	MatrixXd angleToRotMatrix(double angle, int size) {
		MatrixXd rot(size, size);
		rot(0, 0) = cos(angle * PI / 180.0);
		rot(0, 1) = -sin(angle * PI / 180.0);
		rot(1, 0) = sin(angle * PI / 180.0);
		rot(1, 1) = cos(angle * PI / 180.0);

		return rot;
	}
};


#endif ICP_H