#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include <cstring>

#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace std;
using namespace Eigen;

typedef Matrix<double, 2, Dynamic> VertexList;

#define N 1000

typedef bool(*CollisionFunction)(const Vector3d &a, const Vector3d &b);

/* Generate N-dimensional normal distributed samples */
class MultivariateNormalDistribution : normal_distribution<double> {

protected:
	mt19937 &gen;
	
	Matrix<double, Dynamic, Dynamic> covar;
	Matrix<double, Dynamic, Dynamic> transform;
	Matrix<double, Dynamic, 1> mean;

	// drawback: this creates a useless eigenSolver when using Cholesky decomposition, but it yields access to eigenvalues and vectors
	

public:
	MultivariateNormalDistribution(const Matrix<double, Dynamic, 1> &m, const Matrix<double, Dynamic, Dynamic> &c, mt19937 &g) :
		gen(g),
		mean(m),
		covar(c),
		normal_distribution<double>(0, 1)
	{
		LLT<Matrix<double, Dynamic, Dynamic> > cholSolver(covar);
		// We can only use the cholesky decomposition if 
		// the covariance matrix is symmetric, pos-definite.
		// But a covariance matrix might be pos-semi-definite.
		// In that case, we'll go to an EigenSolver
		if (cholSolver.info() == Success) {
			// Use cholesky solver
			transform = cholSolver.matrixL();
		}
		else {
			SelfAdjointEigenSolver<Matrix<double, Dynamic, Dynamic> > eigenSolver(covar);
			transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
		}
	}

	Vector3d randc()
	{
		Vector3d p, q;

		/* Get vector of standard normal distributed samples */
		for (int i = 0; i < p.rows(); i++)
			p[i] = operator()(gen);
		
		return transform * p + mean;
	}
};

/* Is point "t" in polygon "v"?
 * From: https://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html */
bool pnpoly(const VertexList &v, const Vector2d &t)
{
	int i, j;
	bool c = false;

	for (i = 0, j = v.cols()-1; i < v.cols(); j = i++) {
		if (
			((v(1,i) > t(1)) != (v(1,j) > t(1))) &&
			(t(0) < (v(0,j) - v(0,i)) * (t(1) - v(1,i)) / (v(1,j) - v(1,i)) + v(0,i))
		)
			c = !c;
	}

	return c;
}

/* Calculate vertices of car polygon
 * vertx, verty will be filled */
VertexList poly_vehicle(const Vector3d &car)
{
	Rotation2D<double> rot(car[2]);

	VertexList vertices(2, 4);
	Vector2d offset = car.head(2);
	
	double scale = 1.5;
	
	vertices.col(0) <<  2*scale,  1*scale;
	vertices.col(1) << -2*scale,  1*scale;
	vertices.col(2) << -2*scale, -1*scale;
	vertices.col(3) <<  2*scale, -1*scale;
	
	for (int i = 0; i < vertices.cols(); i++)
		vertices.col(i) = rot * vertices.col(i) + offset;
	return vertices;
}

bool collision_vehicle(const Vector3d &vehicle, const Vector3d &obstacle)
{
	/* Create list of X/Y vertex points of car */
	VertexList vertices_vehicle, vertices_obstacle;
	
	vertices_vehicle = poly_vehicle(vehicle);
	vertices_obstacle = poly_vehicle(obstacle);
	
	/* Check if any of the vertices of the obstacle car is inside the vehicle polygon */
	for (int i = 0; i < vertices_obstacle.cols(); i++) {
		if (pnpoly(vertices_vehicle, vertices_obstacle.col(i)) > 0)
			return true;
	}

	return false;
}

bool collision_circle(const Vector3d &vehicle, const Vector3d &obstacle)
{
	const double minDistance = 5;
	
	Vector3d delta = vehicle - obstacle;

	return delta.head(2).norm() < minDistance; /* Collision occured if distance (r) is smaller than 1 */
}

double collision_probability(CollisionFunction cfunc, MultivariateNormalDistribution &vehicle, MultivariateNormalDistribution &obstacle)
{
	unsigned collisions = 0;
	
	for (unsigned j = 0; j < N; j++) {
		Vector3d sample_vehicle, sample_obstacle;
		VertexList vertices_vehicle, vertices_obstacle;

		sample_vehicle = vehicle.randc();
		sample_obstacle = obstacle.randc();
		
		vertices_vehicle = poly_vehicle(sample_vehicle);
		vertices_obstacle = poly_vehicle(sample_obstacle);
		
		bool coll = cfunc(sample_vehicle, sample_obstacle);
		
		/* Print vehicle centers */
		cout << ((coll) ? 1 : 0) << " "
			<< sample_vehicle.transpose() << " "
			<< sample_obstacle.transpose();
		
		/* Print vertices of vehicles */
		for (int i = 0; i < 4; i++)
			cout << " " << vertices_vehicle.col(i).transpose();
		for (int i = 0; i < 4; i++)
			cout << " " << vertices_obstacle.col(i).transpose();
		
		cout << endl;
		
		if (coll)
			collisions++;
	}
	
	cout << "# collisions = " << collisions << endl;
	
	return (double) collisions / N;
}

int main(int argc, char *argv[]) {
	/* Initialize PRNG */
	random_device rd;
	mt19937 gen(rd());

	Vector3d mean_vehicle, mean_obstacle;
	Matrix3d covar_vehicle, covar_obstacle;
		
	mean_vehicle  <<  5, 0, 3.0/4.0 * M_PI; // x, y, heading
	mean_obstacle <<  0, 0, 1.0/4.0 * M_PI; // x, y, heading
	
	covar_vehicle  << 5, -4,  0,
	                 -4,  9,  0,
	                  0,  0,  1.0/16.0*M_PI;

	covar_obstacle << 3,  2,  0,
	                  2,  3,  0,
	                  0,  0,  1.0/16.0*M_PI;
	
	CollisionFunction cfunc = collision_vehicle;
	
	MultivariateNormalDistribution vehicle(mean_vehicle, covar_vehicle, gen);
	MultivariateNormalDistribution obstacle(mean_obstacle, covar_obstacle, gen);
	
	cout << "# Collision Vehicle_x Vehicle_y Obstacle_x Obstacle_y" << endl;
	
	double prob = collision_probability(cfunc, vehicle, obstacle);
	
	cout << "# collision probability is: " << prob << endl;
}