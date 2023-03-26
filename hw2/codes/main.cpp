#include <stdio.h>
#include <iostream>
#include <Fade_2D.h>
#include <cmath>
#include <tuple>
#include <armadillo>
#include <assert.h>

#include "solver.hpp"


void write_solution(const char fname[], vector<double> x, vector<double> y, vector<double> u) {
    FILE* f = fopen(fname, "w");
    if (f == NULL) {
        perror("fopen");
    }
	int err;
	for (size_t i=0; i<u.size(); i++)
	{
		err = fprintf(f, "%lf %lf %lf\n", x[i], y[i], u[i]);
		if (err < 0) fprintf(stderr, "Error when writing. Do you have enough space left ?");
	}
    err = fclose(f);
    if (err != 0) perror("fclose");
}



int main()
{

	/*
		Define problem functions
	*/
	// diffusivity coefficient
	auto a = [](double x, double y)
	{
		return sin(0.1*x)+cos(0.1*y);
	};
	// source
	auto f = [](double x, double y)
	{
		return x+0.3*y;
	};
	// dirichlet border condition
	auto g = [](double x, double y)
	{
		return sqrt(pow(x,2)+0.1*pow(y,2))+x;
	};

	/*
		Setup points to define volumes
	*/
	// in or out the problem?
	auto problem_region = [](double x, double y)
	{
		return sqrt(pow(x,2)+pow(y,2))<10-1e-2;
	};
	// grid of points inside the domain
	std::vector<Point2> v_points;
	double delta = 1.;
	double length = 20.;
	for(double x=-length; x<=length+1e-6; x+=delta)
	{
		for(double y=-length; y<=length+1e-6; y+=delta)
		{
			if (problem_region(x,y))
				v_points.push_back(Point2(x,y));
		}
	}
	// points on the border
	for (double alpha=0.; alpha<2*M_PI; alpha+=2*M_PI/1000)
		v_points.push_back(Point2(10.*cos(alpha), 10*sin(alpha)));


	/*
		solve the problem	
	*/
	vector<double> x;
	vector<double> y;
	vector<double> sol;
	solve(a, f, g, v_points, x, y, sol);

	/*
		write it
	*/
	write_solution("sol.txt", x, y ,sol);
}

