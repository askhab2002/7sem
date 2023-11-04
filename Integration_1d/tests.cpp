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
#define eps1 1e-5


using namespace std;

double derivative(function_pointer f, double x);
double derivative_N(function_pointer f, double x, int N);

double precision_gauss(function_pointer f, double a, double b);
double precision_simpson(function_pointer f, double a, double b);
double precision_gauss_N(function_pointer f, double a, double b, int N);
double precision_simpson_N(function_pointer f, double a, double b, int N);

void test(function_pointer func, double a, double b, int N, double true_value, string file_txt, string plot, int kee);

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

	cout << "   Введите способ интегрирования. 0 - Симпсон, 1 - Гаусс" << endl;
	int kee = 0;
	cin >> kee;

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

	cout << " simpson N partition = " << precision_simpson_N(func, a, b, N) << endl;
        cout << " gauss N partition = " << precision_gauss_N(func, a, b, N) << endl;

	fstream out;
        out.open("ps.txt", std::ofstream::out | std::ofstream::trunc);

	for(int i = 0; i < N; i++) {
		out << i + 1 << " " << precision_simpson_N(func, a, b, i + 1) << endl;
        }

	out.close();

	fstream out1;
        out1.open("pg.txt", std::ofstream::out | std::ofstream::trunc);

        for(int i = 0; i < N; i++) {
                out1 << i + 1 << " " << precision_gauss_N(func, a, b, i + 1) << endl;
        }

        out1.close();

	func = func9_x;
	test(func, 0, 1, N, 0.1, "func_9x.txt", "plot1.gpi", kee);
	func = cosin;
	test(func, 0, M_PI, N, 0, "cosin.txt", "plot2.gpi", kee);
	func = exponen;
	test(func, 0, 1, N, 0.001, "exponen.txt", "plot3.gpi", kee);
	func = arcsin;
	test(func, -1 + eps1, 1 - eps1, N, M_PI, "arcsin.txt", "plot4.gpi", kee);
/*
        FILE *file = fopen("plot.gpi", "w+");

        fprintf(file, "set terminal png size 500, 500 enhanced font \"Helvetica Bold, 10\" \n");
        fprintf(file, "set output \"result.png\" \n");
        fprintf(file, "set xlabel \"Partition\" \n");
        fprintf(file, "set ylabel \"Error\" \n");
        fprintf(file, "set grid \n");
        fprintf(file, "set title \"Integration 1d, partition number = %d \" font \"Helvetica Bold, 10\" \n", N);
        fprintf(file, "plot \"func_9x.txt\" with line lc rgb \"blue\", \"cosin.txt\" with line lc rgb \"red\" ,  \"exponen.txt\" with line lc rgb \"black\", \"arcsin.txt\" with line lc rgb \"grey\" \n");

        fclose(file);
*/
	return 0;
}

void test(function_pointer func, double a, double b, int N, double true_value, string file_txt, string plot, int kee) {
	double Gauss = IntegralGauss(a, b, func);
        double Gauss_N = IntegralGauss_N(a, b, func, N);

        double Simpson = IntegralSimpson(a, b, func);
        double Simpson_N = IntegralSimpson_N(a, b, func, N);
        
	cout << file_txt << ": " << endl;
        cout << "   Гаусс простой:  " << Gauss << endl;
        cout << "   Гаусс по " << N << " разбиению:  " << Gauss_N << endl;
        cout << "   Симпсон простой:  " << Simpson << endl;
        cout << "   Симпсон по " << N << " разбиению:  " << Simpson_N << endl;

	fstream out;
        out.open(file_txt, std::ofstream::out | std::ofstream::trunc);
        
	if(kee == 0) {
                for(int i = 0; i < N; i++) {
                        out << i + 1 << " " << fabs(IntegralSimpson_N(a, b, func, i + 1) - true_value) << endl;
                }
	}

	else {
		for(int i = 0; i < N; i++) {
                        out << i + 1 << " " << fabs(IntegralGauss_N(a, b, func, i + 1) - true_value) << endl;
                }
	}

        out.close();

	FILE *file = fopen(plot.c_str(), "w+");

        fprintf(file, "set terminal png size 500, 500 enhanced font \"Helvetica Bold, 10\" \n");
        fprintf(file, "set output \"result.png\" \n");
        fprintf(file, "set xlabel \"Partition\" \n");
        fprintf(file, "set ylabel \"Error\" \n");
        fprintf(file, "set grid \n");
        fprintf(file, "set title \"Integration 1d, partition number = %d \" font \"Helvetica Bold, 10\" \n", N);
        fprintf(file, "plot \"%s\" with line lc rgb \"blue\" \n", file_txt.c_str());

        fclose(file);

//      cout << derivative_N(func, (a + b)/2, 6) << endl;
/*        cout << " simpson = " << precision_simpson(func, a, b) << endl;
        cout << " gauss = " << precision_gauss(func, a, b) << endl;

        cout << " simpson N partition = " << precision_simpson_N(func, a, b, N) << endl;
        cout << " gauss N partition = " << precision_gauss_N(func, a, b, N) << endl;
*/
}

double precision_gauss_N(function_pointer f, double a, double b, int N) {
	double norm = 0;
        double value = 0;

        for(int i = 0; i <= 100; i++) {
                value = derivative_N(f, a + i * (b - a)/((double)100), 6);
                if(value > norm) {
                        norm = value;
                }
        }

        return (norm * pow(b - a, 7))/(2016000 * pow((double) N, 6));
}

double precision_simpson_N(function_pointer f, double a, double b, int N) {
        double norm = 0;
        double value = 0;

        for(int i = 0; i <= 100; i++) {
                value = derivative_N(f, a + i * (b - a)/((double)100), 4);
                if(value > norm) {
                        norm = value;
                }
        }

        return (norm * pow(b - a, 5))/(2880 * pow( (double)N, 4));
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


