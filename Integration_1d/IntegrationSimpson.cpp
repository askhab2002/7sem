#include "function_pointer.hpp"
#include "IntegrationSimpson.hpp"
#include <iostream>

double IntegralSimpson(double a, double b, function_pointer f) { 
//	std::cout << f(a) << std::endl;
	return ((b - a)/6) * (f(a) + 4 * f((a + b)/2) + f(b));
}

double IntegralSimpson_N(double a, double b, function_pointer f, int N) {
	double sum = 0;

	for(int i = 0; i < N; i++) {
		sum += IntegralSimpson(a + i * ((b - a)/N), a + (i + 1) * ((b - a)/N), f);
	}

	return sum;
}

