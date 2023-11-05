#include <iostream>
#include "functions.hpp"

double function2(double x1, double x2) {
	return 3 * x1 * x1  + 10 * x1 * x2 + x2 * x2;
}

double function1(double x1, double x2) {
        return x1 * x1 * x1 * x1 + x1 * x1 * x2 * x2 + x2 * x2 * x2 * x2;
}

double function3(double x1, double x2) {
	
	x2 = 0.;
        return x1 + x2;
}
