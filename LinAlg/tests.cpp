#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>

#define eps 1e-5


typedef double (* function_pointer)(double x);

double Scalar(function_pointer f, int m, int N);
double Scalar_(int m, int N);
double Lambda(int n, double p, int N);
double C_Coeff(int m, int N, function_pointer f, double p);
double Y_Coeff(int m, int N, function_pointer f, double p);
double *TriangularSystem(int N, function_pointer f, double p);

double function1(double x);
double function2(double x);

using namespace std;

int main(void) {
         
	function_pointer f = function2;
        cout << "    Введите число узлов" << endl;
        
	int N = 0;
	cin >> N;

	double p = 0;
	cout << "    Введите число p" << endl;
        cin >> p;

	double *Y = TriangularSystem(N, f, p);
/*
	for(int i = 0; i < N; i++) {
		cout << Y[i] << endl;
	}
*/
	free(Y);

	return 0;

}	


double function1(double x) {

	return sin(M_PI * x);
}

double function2(double x) {
        if(fabs(x - 1) < eps) {
                return 0;
        }

        if(fabs(x) < eps) {
                return 0;
        }

        return exp(1/((2 * x - 1) * (2 * x - 1) - 1));
}




double Scalar(function_pointer f, int m, int N) {

	double sum = 0;

	for(int i = 1; i < N + 1; i++) {
		sum += sin((M_PI * m * i)/((double)N + 1)) * f(i/((double)N + 1));
	}
//        cout << m << "-----" << sum << endl;
	return sum;
}

double Scalar_(int m, int N) {
        
        double sum = 0;
        
        for(int i = 1; i < N + 1; i++) {
                sum += sin((M_PI * m * i)/((double)N + 1)) * sin((M_PI * m * i)/((double)N + 1));
        }       

        return sum;
}

double Lambda(int n, double p, int N) {
	return p - 2 * n * n * (cos((M_PI * n)/((double)N + 1)) - 1);
}

double C_Coeff(int m, int N, function_pointer f, double p) {
	
	return (Scalar(f, m, N))/(Lambda(m, p, N) * Scalar_(m, N));
}

double Y_Coeff(int m, int N, function_pointer f, double p) {
        double sum = 0;

	for(int i = 1; i < N + 1; i++) {
		sum += C_Coeff(i, N, f, p) * sin((M_PI * m * i)/((double)N + 1));
	}

	return sum;
}

double *TriangularSystem(int N, function_pointer f, double p) {

	double *Y = (double *)calloc(N, sizeof(double));

	for(int i = 1; i < N + 1; i++) {
		Y[i - 1] = Y_Coeff(i, N, f, p);
	}

	return Y;
}

