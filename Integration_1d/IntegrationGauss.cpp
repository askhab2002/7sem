#include <math.h>
#include "function_pointer.hpp"
#include "IntegrationGauss.hpp"


double IntegralGauss(double a, double b, function_pointer f) {
	double x1 = (a + b)/2 + ((b - a)/2) * sqrt((double)5/(double)3);
	double x3 = (a + b)/2 - ((b - a)/2) * sqrt((double)5/(double)3);
	double x2 = (a + b)/2;

	return ((b - a)/18) * (5 * f(x1) + 8 * f(x2) + 5 * f(x3));
}	

double IntegralGauss_N(double a, double b, function_pointer f, int N) {
	double sum = 0;

        for(int i = 0; i < N; i++) {
                sum += IntegralGauss(a + i * ((b - a)/N), a + (i + 1) * ((b - a)/N), f);
        }

        return sum;
}
