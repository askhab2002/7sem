#define _USE_MATH_DEFINES
#define eps 1e-6

#include <iostream>
#include <time.h>
#include <random>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "function_pointer.hpp"

using namespace std;

double *generate_equel(double left, double right, int num);


double func1(double x);
double func2(double x);
double func3(double x);
double func4(double x);
double cosinus(double x);

void function_out(function_pointer func, double *values, int values_number, string filename);
void fourier_out(function_pointer2 func, double *coeff, int num, double *values, int values_number, string filename);
double *coeff_out(function_pointer func, function_pointer cosin, double *values, int values_number);
double fourier(double *coeff, int num, double x);
double Scalar(function_pointer func, function_pointer cosin, double h, int m, int n);

double  *test(int node_number, function_pointer func);
double constant_p(double *coeff, int n, function_pointer func);

int main(int argc, char *argv[]) {

	if(argc != 3) {
                cout << "    Вы ввели некорректную командную строку" << endl;
                return 0;
        }


	
        string func_type_(argv[1]);
        int func_type = stoi(func_type_);
	string node_start_(argv[2]);
        int node_start = stoi(node_start_);

	function_pointer func;

        if(func_type == 1) {
                func = &func1;
        }

        if(func_type == 2) {
                func = func2;
        }

        if(func_type == 3) {
                func = func3;
        }

        if(func_type == 4) {
                func = func4;
        }

/*	string node_step_(argv[3]);
	int node_step = stoi(node_step_);
	string node_step_number_(argv[4]);
	int node_step_number = stoi(node_step_number_);

	int node_number = node_start;

	for(int i = 1; i < node_step_number + 1; i++) {
		node_number += node_step;
		test(node_number, func_type);
	}
*/      
	double *coeff = test(node_start, func);
        
        double p = constant_p(coeff, node_start, func);

	cout << " p = " << p << endl;
	
        free(coeff);

	return 0;
}


double *test(int node_number, function_pointer func) {

	function_pointer cosin = cosinus;

//	cout << " Введите границы отрезка" << endl;
        double left = 0;
        double right = 1;
//        cin >> left;
//        cin >> right;

	int part = 0;
        part = 15 * node_number;

	double *nodes = NULL;
	nodes = generate_equel(left, right, node_number);

	double *partition = generate_equel(left, right, part);
        
        double *coeff = coeff_out(func, cosin, nodes, node_number);


	function_pointer2 four;
	four = fourier;

	function_out(func, partition, part, "func.txt");
	fourier_out(four, coeff, node_number, partition, part, "fourier.txt");
	
	FILE *file = fopen("plot.gpi", "w+");

        fprintf(file, "set terminal png size 500, 500 enhanced font \"Helvetica Bold, 10\" \n");
        fprintf(file, "set output \"result.png\" \n");
        fprintf(file, "set xlabel \"X\" \n");
        fprintf(file, "set ylabel \"Y\" \n");
        fprintf(file, "set grid \n");
        fprintf(file, "set title \"Fourier, nodes number = %d \" font \"Helvetica Bold, 10\" \n", node_number);
        fprintf(file, "plot \"func.txt\" with line lc rgb \"blue\", \"fourier.txt\" with line lc rgb \"red\" \n");

        fclose(file);

//        free(coeff);
        free(partition);
        free(nodes);
       	

	return coeff;
}

double constant_p(double *coeff, int n, function_pointer func) {
        function_pointer2 four = fourier;

	double norm = 0;
	double h = 0;

	double *norms = (double *)calloc(10, sizeof(double));
	double *hs = (double *)calloc(10, sizeof(double));

	for(int i = 1; i <= 10; i++) {
	        h = 1 / ((double)n * (double)i - 0.5);

		for(double k = 0; k < 1; k += h) { 
	                if(fabs(func(h) - four(coeff, n, h)) > norm) {
		                norm = fabs(func(h) - four(coeff, n, h));
		        }
		}

		norms[i - 1] = log(1/norm);
		hs[i - 1] = log(1/h);
		norm = 0;
	}

	fstream out;
        out.open("p.txt", std::ofstream::out | std::ofstream::trunc);

        for(int i = 0; i < 10; i++) {
                out << hs[i] << " " << norms[i] << endl;

        }

        out.close();

	double p = 0;

	for(int i = 1; i < 10; i++) {
		p += (norms[i] - norms[i - 1])/(hs[i] - hs[i - 1]);
		
	}


	FILE *file = fopen("2plot.gpi", "w+");

        fprintf(file, "set terminal png size 500, 500 enhanced font \"Helvetica Bold, 10\" \n");
        fprintf(file, "set output \"p.png\" \n");
        fprintf(file, "set xlabel \"X\" \n");
        fprintf(file, "set ylabel \"Y\" \n");
        fprintf(file, "set grid \n");
        fprintf(file, "set title \"Fourier p constant = %lf,  nodes number = %d \" font \"Helvetica Bold, 10\" \n", p/9, n);
        fprintf(file, "plot \"p.txt\" with line lc rgb \"blue\"  \n");

        fclose(file);

        free(norms);
	free(hs);

	return p/9;


}

double *generate_equel(double left, double right, int num) {
        double step = 1/(num - 0.5);
        double *nodes = (double *)calloc( num - 1, sizeof(double) );

        if(num >= 2) {
                nodes[0] = left;
        }

        for(int i = 1; i < num - 1; i++) {
               nodes[i] = nodes[i - 1] + step;
        }

        return nodes;
}

double func1(double x) {
        return cos(x * 10 * M_PI);
}

double func2(double x) {
        return sin((x + M_PI/2) * M_PI);
}

double func3(double x) {
	if(fabs(x - 1) < eps) {
		return 0;
	}

	return exp(1/(x * x - 1));
}

double func4(double x) {
        if(fabs(x - 1) < eps) {
                return 0;
        }

	if(fabs(x) < eps) {
		return 0;
	}

        return exp(1/((2 * x - 1) * (2 * x - 1) - 1));
}

double cosinus(double x) {
	return cos(x);
}

void function_out(function_pointer func, double *values, int values_number, string filename) {
        fstream out;
        out.open(filename, std::ofstream::out | std::ofstream::trunc);

        for(int i = 0; i < values_number - 1; i++) {
                out << values[i] << " " << (*func)(values[i]) << endl;

        }

        out.close();
}

void fourier_out(function_pointer2 func, double *coeff, int num, double *values, int values_number, string filename) {
	fstream out;
        out.open(filename, std::ofstream::out | std::ofstream::trunc);

        for(int i = 0; i < values_number - 1; i++) {
                out << values[i] << " " << (*func)(coeff, num, values[i]) << endl;

        }

        out.close();
}

double Scalar(function_pointer func, function_pointer cosin, double h, int m, int n) {
        double scalar = 0;
        
	scalar += (*func)(0) * (*cosin)(0) / 2;
//        cout << " f(0) = " << (*func)(0) << endl;
//	cout << " f(1) = " << (*func)( h * (n - 1)) << endl;
	for(int k = 1; k < n; k++) {
		scalar += (*func)(k * h) * (*cosin)(M_PI * k * h * m);
	}

        cout << " m = " << m << " scalar = " << scalar << endl;
	return scalar * h;
}


double *coeff_out(function_pointer func, function_pointer cosin, double *values, int values_number) {

	double *coeff = (double *) calloc(values_number , sizeof(double));
        double scalar = 0;
	double h = 1/((double)values_number - 0.5);
	cout << " h = " << h << endl;
        
	

	for(int m = 0; m < values_number ; m++) {
		scalar = Scalar(func, cosin, h, m, values_number);
		
		coeff[m] = scalar * 2;
		cout << " coeff = " << coeff[m] << endl;
	} 
        
	double sum = 0;
        for(int i = 1; i < values_number; i++) {
	        sum += coeff[i];
	}
        coeff[0] = (*func)(0) - sum;
        cout << coeff[0] << endl;
	return coeff;
}

double fourier(double *coeff, int num, double x) {

	double sum = 0;

	for(int i = 0; i < num ; i++) {
		sum += coeff[i] * cos(M_PI * i * x);
	}

	return sum;
}
