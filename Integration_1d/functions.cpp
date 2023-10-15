#include <math.h>
#include "functions.hpp"

double func0_x(double x) {
        return 1;
}

double func1_x(double x) {
        return x;
}

double func2_x(double x) {
        return x * x;
}

double func3_x(double x) {
        return x * x * x;
}

double func5_x(double x) {
        return x * x * x * x * x;
}

double func9_x(double x) {
        return x * x * x * x * x * x * x * x * x;
}

double cosin(double x) {
	return cos(100 * x);
}

double exponen(double x) {
	return exp(-1000 * x);
}

double arcsin(double x) {
	return 1/sqrt(1 - x * x);
}

