#include <iostream>
#include "functions.hpp"

double function1(double x1, double x2) {
	return x1 * x1  + x1 * x2 + x2 * x2;
}

double function2(double x1, double x2) {
        return x1 * x1 * x1 * x1 + x1 * x1 * x2 * x2 + x2 * x2 * x2 * x2;
}
