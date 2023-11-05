#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include "FunctionPointer.hpp"
#include "Integration_2d.hpp"
#include "Triangulation.hpp"
#include "functions.hpp"

#define eps 1e-2


using namespace std;

double test(function_pointer f, int N);

int main(void) {
        int N = 0;
	cout << " Введите разбиение стороны " << endl;
	cin >> N;

	function_pointer f = function1;

	FILE *file1 = fopen("results.txt", "w+");

        for(int i = 1; i < N; i++) {
                fprintf(file1, "%d %lf\n", i + 1, test(f, i + 1));
        }

        fclose(file1);

        FILE *file = fopen("plot.gpi", "w+");

        fprintf(file, "set terminal png size 500, 500 enhanced font \"Helvetica Bold, 10\" \n");
        fprintf(file, "set output \"result.png\" \n");
        fprintf(file, "set xlabel \"Partition\" \n");
        fprintf(file, "set ylabel \"Error\" \n");
        fprintf(file, "set grid \n");
        fprintf(file, "set title \"Integration 2d, partition number = %d \" font \"Helvetica Bold, 10\" \n", N);
        fprintf(file, "plot \"results.txt\" with line lc rgb \"blue\" \n");

        fclose(file);

	return 0;
}

double test(function_pointer f, int N) {

	double ***triangles = (double ***) malloc(2 * N * N * sizeof(double **));
        double **nodes = (double **) malloc((N + 1) * (N + 1) * sizeof(double *));

        FILE *out;
	out = fopen("Triangulation.txt", "w+");

        Triangulation(1, 1, N, out, triangles, nodes);

	fclose(out);
/*
	for(int i = 0; i < 2 * N * N; i++) {
		cout << i + 1 << ": ";
		for(int j = 0; j < 3; j++) {
			cout << "(" << triangles[i][j][0] << ", " << triangles[i][j][1] << ")  ";
		}
		cout << endl;
	} 

  */     

//	cout << "    Интеграл по разбиению на " << 2 * N * N << " треугольников = " << Integral_2d_N(f, N, triangles) << endl;
 	
	double I = Integral_2d_N(f, N, triangles);

	for(int i = 0; i < (N + 1) * (N + 1); i++) {
		free(nodes[i]);
	}
	free(nodes);

        for(int i = 0; i < 2 * N * N; i++) {
/*	        for(int j = 0; j < 3; j++) {
		        free(triangles[i][j]);
		} */
	        free(triangles[i]);
	}
        free(triangles);

	return I;

 
}
