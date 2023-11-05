#include <math.h>
#include <iostream>
#include "FunctionPointer.hpp"
#include "Triangulation.hpp"
#include "Integration_2d.hpp"

double Integral_2d(function_pointer f, double **triangle, double S) {
	double Ax = 0.5 * (triangle[0][0] + triangle[1][0]);
	double Ay = 0.5 * (triangle[0][1] + triangle[1][1]);

	double Bx = 0.5 * (triangle[1][0] + triangle[2][0]);
	double By = 0.5 * (triangle[1][1] + triangle[2][1]);

	double Cx = 0.5 * (triangle[0][0] + triangle[2][0]);
	double Cy = 0.5 * (triangle[0][1] + triangle[2][1]);

        
//	std::cout  << " Ax = " << Ax << " Ay = " << Ay << " Bx = " << Bx << " By = " << By << " Cx = " << Cx << " Cy = " << Cy << std::endl;
	return (S/3) * (f(Ax, Ay) + f(Bx, By) + f(Cx, Cy));
}

double Integral_2d_N(function_pointer f, int N, double ***triangles) {
        double sum = 0;
        double S = 1/(2 * (double)N * (double)N);
//	std::cout << "------------" << S << std::endl;

	for(int i = 0; i < 2 * N * N; i++) {
		sum += Integral_2d(f, triangles[i], S);
//		std::cout << Integral_2d(f, triangles[i], S) << std::endl;
	}
//        std::cout << sum << std::endl; 
	return sum;

}


