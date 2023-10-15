#include <stdio.h>
#include <iostream>
#include <cmath>
#include <math.h>
#include "function_pointer.hpp"
#include "functions.hpp"
#include "IntegrationSimpson.hpp"
#include "IntegrationGauss.hpp"

using namespace std;

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
        

	return 0;
}
