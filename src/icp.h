#ifndef ICP_H
#define ICP_H

#include <iostream>
#include <algorithm>
#include <math.h>
#include <Eigen\Dense>
#include <Eigen\Jacobi>
#include <Eigen\Geometry>

#define PI 3.14159265

namespace Optimum {

	/**
	* Enumeration for if we're in 3D or 2D.
	*/
	enum PointType { TWO_D = 0, THREE_D };

	/**
	* Simple point class for ICP algorithm.
	*/
	class Point {
		std::vector<double> coords;
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
			std::vector<double> sqs;
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

	struct PointMatcherSettings {
		PointType dimension = TWO_D;
		bool doScale = false;
	};

	struct Transform {
		Eigen::MatrixXd rotation;
		Eigen::VectorXd translation;
		double scale;
	};
	/**
	* Helper class with static functions to calculate optimal rotation and translation using
	* the Kabsch method.
	*/
	class PointMatcher {

	public:

		PointMatcher(PointMatcherSettings settings) {
			this->mSettings = settings;
		}

		/**
		* Solves for the optimal transformation and rotation (and scale) to map the target data point
		* map to the reference data point map. It is assumed that the mapping of this points are one
		* to one. I.e. the first point in the target is known to map to the first point of the reference. 
		* @param target the target points to align.
		* @param ref the reference points to align to.
		* @return Transform a transformation that contains:
		*		Transform.rotation: the rotation matrix
		*		Transform.translation: the translation vector.
		*		Transform.scale: the scale factor. 
		*/
		Transform solve(Eigen::MatrixXd target, Eigen::MatrixXd ref) {
			//The transform we will return.
			Transform transform;

			Eigen::MatrixXd centroidOne = getCentroid(target);
			Eigen::MatrixXd centroidTwo = getCentroid(ref);

			//center points at centroid
			int rows = ref.rows();
			Eigen::MatrixXd centeredRef(ref.rows(), ref.cols());
			Eigen::MatrixXd centeredTar(target.rows(), target.cols());
			double f_sd2 = 0;
			for (int r = 0; r < rows; r++) {
				//transpose the centroids because they're column vectors.
				//each row is a point.
				centeredRef.row(r) = ref.row(r) - centroidTwo.transpose();
				centeredTar.row(r) = target.row(r) - centroidOne.transpose();

				f_sd2 += centeredTar.row(r).squaredNorm();
			}

			//calculate covarience matrix to align A (target) to B (reference)
			// H = sum ((Pa - centroidA)*(Pb - centroidB)^T)
			//H should be 3x3
			Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(ref.cols(), ref.cols());
			for (int i = 0; i < rows; i++) {
				cov.noalias() += centeredTar.row(i).transpose() * centeredRef.row(i);
			}

			//calculate svd
			Eigen::JacobiSVD<Eigen::MatrixXd> svd(cov, Eigen::QRPreconditioners::FullPivHouseholderQRPreconditioner |
				Eigen::ComputeFullU | Eigen::ComputeFullV);

			//if svd(H) = [U, S, V], then rotation R = V*U(transpose)
			//eigen returns SVD(A) = USV*

			Eigen::MatrixXd rotation = svd.matrixV() * svd.matrixU().transpose();

			//check for correction.
			double f_det = rotation.determinant();
			Eigen::VectorXd e = Eigen::VectorXd::Ones(rotation.rows());
			if (f_det < 0) {
				int rows = e.rows() - 1;
				int cols = e.cols() - 1;
				e(rows, cols) = -1;
			}

			if (f_det < 0) {
				//recompute the rotation part if the determinant was negative
				rotation.noalias() = svd.matrixV() * e.asDiagonal() * svd.matrixU().transpose();
			}

			//renormalizes the rotation matrix, more orthogonal results. 
			//only works in 3d. 
			if (mSettings.dimension == THREE_D) {
				rotation = Eigen::Quaterniond((Eigen::Matrix3d)rotation)
					.normalized().toRotationMatrix();
			}
			
			double f_inv_scale = svd.singularValues().dot(e) / f_sd2;
			if (mSettings.doScale) {
				rotation *= f_inv_scale;
			}
			//centroid one is the target centroid
			//centroid two is the ref centroid
			Eigen::MatrixXd translation = centroidTwo - (rotation * centroidOne);
			transform.rotation = rotation;
			transform.translation = translation;
			transform.scale = f_inv_scale;
			return transform;
		}

		/**
		* Calculates the root mean square error between two matricies.
		* @param dataOne the first matrix in the comparison.
		* @param dataTwo the second matrix in the comparison.
		* @return double the root mean square error or 0 if the matrix sizes don't match.
		*/
		static double RMSE(Eigen::MatrixXd dataOne, Eigen::MatrixXd dataTwo) {
			if (dataOne.size() != dataTwo.size()) {
				std::cout << "Error, matrix size mismatch." << std::endl;
				return 0.0f;
			}
			else {
				//first take difference of the two arrays. 
				Eigen::MatrixXd result = dataOne - dataTwo;
				//create a container for the squares. 
				Eigen::MatrixXd squares = Eigen::MatrixXd::Zero(result.rows(), result.cols());
				for (int i = 0; i < result.rows(); i++) {
					for (int j = 0; j < result.cols(); j++) {
						double val = result(i, j);
						double newVal = val*val;
						squares(i, j) = newVal;
					}
				}

				//vector for the sum of the rows. 
				std::vector<double> sumSqrs;
				for (int r = 0; r < squares.rows(); r++) {
					double sum = squares.row(r).sum();
					//push back the square root of each row sum. 
					sumSqrs.push_back(sum);
				}

				//get the square roots. 
				std::vector<double> sqrts;
				for (int q = 0; q < sumSqrs.size(); q++) {
					sqrts.push_back(std::sqrt(sumSqrs.at(q)));
				}

				//calcuate rmse as the sum of square roots. 
				double rmse = 0.0;
				for (int x = 0; x < sqrts.size(); x++) {
					rmse += sqrts.at(x);;
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
		static Eigen::MatrixXd applyTransformation(Eigen::MatrixXd data, Eigen::MatrixXd translation, Eigen::MatrixXd rotation) {
			Eigen::MatrixXd trans(data.rows(), data.cols());
			Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(data.rows(), 1);
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

			Eigen::MatrixXd transData = data + trans;
			Eigen::MatrixXd transformed(transData.rows(), transData.cols());
			for (int i = 0; i < transData.rows(); i++) {
				transformed.row(i) = rotation * transData.row(i).transpose();
			}
			return transformed;
		}

		/**
		* Returns the centroid of a given matrix. The result will be a column vector with the average of each
		* column in the data corresponding to each row of the column vector. So if your input is of the size n X i, then
		* this function will return a column vector of size i x 1 with the number of each row equal to the mean of the ith column in your
		* data.
		* @param points the data.
		* @return MatrixXd a column vector with the centroid points.
		*/
		static Eigen::MatrixXd getCentroid(Eigen::MatrixXd points) {

			int columns = points.cols();
			int rows = points.rows();
			//create a column vector. 
			Eigen::MatrixXd centroid(columns, 1);
			for (int i = 0; i < columns; i++) {
				double sum = points.col(i).sum();
				double cent = sum / rows;
				centroid(i, 0) = cent;
			}
			return centroid;
		}

		private:
			PointMatcherSettings mSettings;

	};

	

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

		Eigen::MatrixXd mRef;
		Eigen::MatrixXd mTarget;
		Eigen::MatrixXd curRef;
		Eigen::MatrixXd pose;
		std::vector<Point> refPoints;
		std::vector<Point> tarPoints;
		std::vector<Point> closestPoints;
		int iterations;
		PointMatcher *matcher;

	public:
		/**
		* Iterative closest point algorithm. It is assumed that a row in both the target and
		* reference data is a point.
		* @param reference the ref data.
		* @param target the target points.
		* @param set ICPSettings struct.
		*/
		ICP(Eigen::MatrixXd target, Eigen::MatrixXd reference, ICPSettings set) {
			this->mRef = reference;
			this->mTarget = target;
			this->iterations = set.maxIterations;
			this->curRef = reference;
			PointMatcherSettings pSets;
			pSets.dimension = set.pointType;
			this->matcher = new PointMatcher(pSets);
		}

		void solve() {

			Eigen::MatrixXd rotation = Eigen::MatrixXd::Identity(mRef.cols(), mRef.cols());
			Eigen::MatrixXd transform = Eigen::MatrixXd::Zero(1, mRef.cols());
			
			//pose is a dim+1 x dim matrix (so for 3D it's 4x3) that holds the rotation
			//and the translation. 
			pose = Eigen::MatrixXd(rotation.rows() + transform.rows(), rotation.rows());
			pose << rotation, transform;

			while (iterations > 0) {

				closestPoints.clear();
				tarPoints.clear();
				refPoints.clear();

				Eigen::MatrixXd closest = match();

				//get rotation and translation. 
				Transform t = matcher->solve(closest, mTarget);
				Eigen::MatrixXd newRot = t.rotation;
				Eigen::MatrixXd newTrans = t.translation;

				Eigen::MatrixXd newPose(newRot.rows() + newTrans.transpose().rows(), newRot.rows());
				newPose << newRot, newTrans.transpose();

				Eigen::MatrixXd lastPose = pose;

				//update our pose.
				pose = lastPose.cwiseProduct(newPose);

				Eigen::MatrixXd applyTrans = pose.row(pose.cols());
				Eigen::MatrixXd applyRotation = pose.block(0, 0, pose.cols(), pose.cols());

				//update the target. 
				Eigen::MatrixXd newRef = PointMatcher::applyTransformation(curRef, applyTrans, applyRotation);
				double error = PointMatcher::RMSE(mTarget, newRef);

				//check for termination.
				if (error < 0.0001) {
					break;
				}
				curRef = newRef;

				//work towards ending of loop. 
				iterations -= 1;
			}
		}

		Eigen::MatrixXd getBestTranslation() {
			return pose.row(pose.cols());
		}

		Eigen::MatrixXd getBestRotation() {
			return pose.block(0, 0, pose.cols(),pose.cols());
		}

	private:
		Eigen::MatrixXd match() {
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
			Eigen::MatrixXd closest(closestPoints.size(), closestPoints.at(0).size());
			for (int x = 0; x < closestPoints.size(); x++) {
				Point p = closestPoints.at(x);
				for (int y = 0; y < p.size(); y++) {
					closest(x, y) = p[y];
				}
			}

			return closest;
		}

		double rotMatrixToDegrees(Eigen::MatrixXd rot) {
			return asin(rot(1, 0)) * 180.0 / PI;
		}

		Eigen::MatrixXd angleToRotMatrix(double angle, int size) {
			Eigen::MatrixXd rot(size, size);
			rot(0, 0) = cos(angle * PI / 180.0);
			rot(0, 1) = -sin(angle * PI / 180.0);
			rot(1, 0) = sin(angle * PI / 180.0);
			rot(1, 1) = cos(angle * PI / 180.0);

			return rot;
		}
	};
}
#endif ICP_H