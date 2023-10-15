#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include "function_pointer.hpp"
#include "functions.hpp"
#include "IntegrationSimpson.hpp"
#include "IntegrationGauss.hpp"

#define eps 1e-2

using namespace std;

double derivative(function_pointer f, double x);
double derivative_N(function_pointer f, double x, int N);

double precision_gauss(function_pointer f, double a, double b);
double precision_simpson(function_pointer f, double a, double b);
double precision_gauss_N(function_pointer f, double a, double b, int N);
double precision_simpson_N(function_pointer f, double a, double b, int N);

int main(void) {
	cout << "  Введите номер функции  " << endl;
	int func_type = 0;
	cin >> func_type;

	function_pointer func = func0_x;

	switch(func_type) {
		case 0:
			func = func0_x;
			break;
		case 1:
			func = func1_x;
			break;
		case 2:
			func = func2_x;
			break;
		case 3:
			func = func3_x;
			break;
		case 4:
			func = func5_x;
			break;
		case 5:
			func = func9_x;
			break;
		case 6:
			func = cosin;
			break;
		case 7:
			func = exponen;
			break;
		case 8:
			func = arcsin;
			break;
	}

	double a = 0;
	double b = 0;

	cout << "   Введите концы отрезка   " << endl;
	cin >> a;
	cin >> b;

	double Gauss = 0;
	double Simpson = 0;
	double Gauss_N = 0;
	double Simpson_N = 0;

	cout << "   Введите разбиение отрезка   " << endl;
        
	int N = 0;
	cin >> N;

	Gauss = IntegralGauss(a, b, func);
	Gauss_N = IntegralGauss_N(a, b, func, N);

	Simpson = IntegralSimpson(a, b, func);
	Simpson_N = IntegralSimpson_N(a, b, func, N);

	cout << "   Гаусс простой:  " << Gauss << endl;
	cout << "   Гаусс по " << N << " разбиению:  " << Gauss_N << endl;
	cout << "   Симпсон простой:  " << Simpson << endl;
        cout << "   Симпсон по " << N << " разбиению:  " << Simpson_N << endl;
        
//	cout << derivative_N(func, (a + b)/2, 6) << endl;
	cout << " simpson = " << precision_simpson(func, a, b) << endl;
	cout << " gauss = " << precision_gauss(func, a, b) << endl;

//	cout << " simpson N = " << precision_simpson_N(func, a, b, N) << endl;
//        cout << " gauss N = " << precision_gauss_N(func, a, b, N) << endl;

	fstream out;
        out.open("ps.txt", std::ofstream::out | std::ofstream::trunc);

	for(int i = 0; i < N; i++) {
		out << i + 1 << " " << IntegralSimpson_N(a, b, func, i + 1) << endl;
        }

	out.close();

	fstream out1;
        out1.open("pg.txt", std::ofstream::out | std::ofstream::trunc);

        for(int i = 0; i < N; i++) {
                out1 << i + 1 << " " << IntegralGauss_N(a, b, func, i + 1) << endl;
        }

        out1.close();



	return 0;
}

double precision_gauss_N(function_pointer f, double a, double b, int N) {
	double precision = 0;
        
	for(int i = 1; i <= N; i++) {
	        precision += precision_gauss(f, a, a + i * (b - a)/((double)N));
	}

	return precision;
}

double precision_simpson_N(function_pointer f, double a, double b, int N) {
        double precision = 0;

        for(int i = 0; i < N; i++) {
                precision += precision_simpson(f, a, a + i * (b - a)/((double)N));
        }

        return precision;
}


double precision_gauss(function_pointer f, double a, double b) {
	double norm = 0;
        double value = 0;

	for(int i = 0; i <= 100; i++) {
		value = derivative_N(f, a + i * (b - a)/((double)100), 6);
		if(value > norm) {
			norm = value;
		}
	}

	return (norm * pow(b - a, 7))/2016000;


}	

double precision_simpson(function_pointer f, double a, double b) {
        double norm = 0;
        double value = 0;

        for(int i = 0; i <= 100; i++) {
                value = derivative_N(f, a + i * (b - a)/((double)100), 4);
                if(value > norm) {
                        norm = value;
                }
        }

        return (norm * pow(b - a, 5))/2880;


}


double derivative(function_pointer f, double x) {
	double right = f(x + eps);
	double left = f(x - eps);

	

	return (right - left)/(2 * eps);
}

double derivative_2(function_pointer f, double x) {
        double right = derivative(f, x + eps);
        double left = derivative(f, x - eps);

	cout << setprecision(16) << derivative(f, x - eps) << " " << derivative(f, x + eps) << endl;



        return (right - left)/(2 * eps);
}

double derivative_N(function_pointer f, double x, int N) {
        if(N == 1) {
		return derivative(f, x);
	}
        
//	cout << derivative_N(f, x + eps, N - 1) << " " << derivative_N(f, x - eps, N - 1) << endl;
	return (derivative_N(f, x + eps, N - 1) - derivative_N(f, x - eps, N - 1))/(2 * eps);

}


