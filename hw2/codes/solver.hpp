#include <stdio.h>
#include <iostream>
#include <Fade_2D.h>
#include <cmath>
#include <tuple>
#include <armadillo>
#include <assert.h>


using namespace std;
using namespace GEOM_FADE2D;

class Neighbour_Volume
{
	public:
		int idx; // index of the Voronoi cell of the neighbour volume in the Voronoi object
		double x;
		double y; // (x,y) position of the midway point between the volume and its neighbour
		double h; // distance with the neighbour's point
		double l; // length of the edge making the interface between the two volumes

		Neighbour_Volume(int idx_, double x_, double y_, double dist_, double edge_length_) : idx(idx_), x(x_), y(y_), h(dist_), l(edge_length_)
		{}
};


/*
    Solve the defined problem using a finite volume method
        with scalar function a, f and g
        with center of the control volume in v_points

        returns the solution sol[i] at points (x[i], y[i])
*/
void solve(const std::function<double(double,double)>& a, const std::function<double(double,double)>& f, const std::function<double(double,double)>& g, std::vector<Point2> v_points, vector<double> &x, vector<double> &y, vector<double> &sol);

/*
	Compute the finite volumes with center in v_points using Delaunay triangulation and Voronoi cells
        return the results in p_voro, v_voronoi_cells and delaunay triangulation
*/
void compute_volumes(std::vector<Point2> v_points, Voronoi2* &p_voro, std::vector<VoroCell2*> &v_voronoi_cells, GEOM_FADE2D::Fade_2D* &dt);


/*
    Get the neighbours of each volume and compute relevant quantities
        p_voro and v_voronoi_cells are the output of the compute columes function
        return adjacent volumes an adjacency list of voronoi cell -> (cell index, x, y, distance, common edge length)
*/
void treat_volumes(Voronoi2 *p_voro, std::vector<VoroCell2*> &v_voronoi_cells, std::vector<std::vector<Neighbour_Volume*>*> &adjacent_volumes);

/*
    Setup system of equations
*/
void setup_equations(const std::function<double(double,double)>& a, const std::function<double(double,double)>& f, const std::function<double(double,double)>& g, int n_points, std::vector<VoroCell2*> &v_voronoi_cells, std::vector<std::vector<Neighbour_Volume*>*> &adjacent_volumes, arma::SpMat<double> &A, arma::Col<double> &b);


/*
    computes the classical euclidean distance between two points
*/
double distance(Point2 *p0, Point2 *p1);
