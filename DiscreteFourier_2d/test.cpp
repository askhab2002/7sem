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

double  *test(int node_number, double  **func, int second);

double *generate_equel(double left,  int num);
double Scalar(double  **func, function_pointer cosin, double h, int m, int n, int second);

double func1(double x, double y);
double func2(double x, double y);
double func3(double x, double y);
double func4(double x, double y);
double func5(double x, double y);
double cosinus(double x, double y);

double **do_matrix(function_pointer_2d func, double *nodes, int node_number);

double *coeff_out(double **func, function_pointer cosin, double *values, int values_number, int second);
double fourier(double **coeff, int num, double x, double y);
void function_out(function_pointer_2d func, double **coeff_new, int num, double *partition, int part, string filename);

double constant_p(double **coeff, int n, function_pointer_2d func);

int main(int argc, char *argv[]) {

        if(argc != 3) {
                cout << "    Вы ввели некорректную командную строку" << endl;
                return 0;
        }



        string func_type_(argv[1]);
        int func_type = stoi(func_type_);
        string node_number_(argv[2]);
        int node_number = stoi(node_number_);

        function_pointer_2d func;

//	cout << " ------------------------------------------- " << endl;

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

	if(func_type == 5) {
                func = func5;
        }

	double left = 0;

	double *nodes = generate_equel(left, node_number);

//	cout << " ------------------------------------------- " << endl;

	double **func_matrix = do_matrix(func, nodes, node_number);
/*
	for(int i = 0; i < node_number; i++) {
		for(int j = 0; j < node_number; j++) {
			cout << func_matrix[i][j] << " ";
		}
		cout << endl;
	}
*/
//	cout << " ------------------------------------------- " << endl;

        
	double **coeff = (double **) malloc(node_number * sizeof(double *));

	for(int i = 0; i < node_number; i++) {
		coeff[i] = test(node_number, func_matrix, i);
	}
        
//	cout << " ------------------------------------------- " << endl;

	double **coeff_new = (double **) malloc(node_number * sizeof(double *));

        for(int i = 0; i < node_number; i++) {
                coeff_new[i] = test(node_number, coeff, i);
        }

//	cout << " ------------------------------------------- " << endl;
        
	int part = 0;
        part = 10 * node_number;
	double *partition = generate_equel(left, part);
        
//	cout << "+++++++++++++++++++++++++++++++++++++" << endl;
	function_out(func, coeff_new, node_number, partition, part, "func.txt");


//	cout << "+++++++++++++++++++++++++++++++++++++" << endl;

        FILE *out = fopen("plot.gpi", "w+");


	fprintf(out, "#! /usr/bin/gnuplot -persist\n");
        fprintf(out, "set terminal png size 1000,1000 enhanced font \"Helvetica Bold, 20\"\n");
        fprintf(out, "set output \"%s.png\"\n\n", "result");

        fprintf(out, "set style line 1 lt 1 linecolor rgb \"red\" lw 1 pt 1\n");
        fprintf(out, "set style line 2 lt 1 linecolor rgb \"blue\" lw 1 pt 1\n");

    
        fprintf(out, "set xrange [0:1]\n");

        fprintf(out, "set title \"%s - %d knots \"\n", "result", node_number);

        fprintf(out, "set grid\n\n");

        fprintf(out, "splot  \"%s\" using 1:2:3 ls 1 title \"Interpolation Fourier Row\", ", "func.txt");
        fprintf(out, "\"%s\" using 1:2:4 ls 2 title \"Original function\"", "func.txt");

	fclose(out);
       
	double p = constant_p(coeff_new, node_number, func);

        cout << " p = " << p << endl;
	

	for(int i = 0; i < node_number; i++) {
		free(coeff[i]);
		free(func_matrix[i]);
		free(coeff_new[i]);
	}

	free(coeff);
	free(func_matrix);
	free(coeff_new);


        return 0;
}

double constant_p(double **coeff, int n, function_pointer_2d func) {
        function_pointer2_2d four = fourier;

        double norm = 0;
        double h = 0;

        double *norms = (double *)calloc(10, sizeof(double));
        double *hs = (double *)calloc(10, sizeof(double));

        for(int i = 1; i <= 10; i++) {
                h = 1 / ((double)n * (double)i);
  //              cout << " h = " << h << endl;
                for(double k = 0; k < 1; k += h) {
			for(double  l = 0; l < 1; l+= h) {
//				cout << four(coeff, n, h, l) << " ";
                                if(fabs(func(h, l) - four(coeff, n, h, l)) > norm) {
                                        norm = fabs(func(h, l) - four(coeff, n, h, l));
//					cout << fabs(func(h, l) - four(coeff, n, h, l)) << " ";
                                }
			}
                }

                norms[i - 1] = log(1/norm);
                hs[i - 1] = log(1/h);
                norm = 0;
//		cout << " norms = " << norms[i - 1] << " hs = " << hs[i - 1] << endl;
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

void function_out(function_pointer_2d func, double **coeff_new, int num, double *partition, int part, string filename) {

	fstream out;
        out.open(filename, std::ofstream::out | std::ofstream::trunc);

	for(int i = 0; i < part - 1; i++) {
                
		for(int j = 0; j < part - 1; j++) {
			out << partition[i] << " " << partition[j] << " " << fourier(coeff_new, num, partition[i], partition[j]) << " " << (*func)(partition[i], partition[j]) << endl;

		}
	}

	out.close();
}

double fourier(double **coeff, int num, double x, double y) {
	double sum = 0;

        for(int i = 0; i < num ; i++) {
		for(int j = 0; j < num; j++) { 
			sum += coeff[i][j] * cos(M_PI * i * x) * cos(M_PI * j * y);
		}
        }

        return sum;
}

double **do_matrix(function_pointer_2d func, double *nodes, int node_number) {
	double **matrix = (double **) malloc( node_number * sizeof(double *));

	for(int i = 0; i < node_number; i++) {
		matrix[i] = (double *) malloc(node_number * sizeof(double));
		for(int j = 0; j < node_number; j++) {
			matrix[i][j] = (*func)(nodes[i], nodes[j]);
		}
	}

	return matrix;
}



double *test(int node_number, double  **func, int second) {

        function_pointer cosin = cosin;

//      cout << " Введите границы отрезка" << endl;
        double left = 0;
//        cin >> left;
//        cin >> right;

        
        double *nodes = NULL;
        nodes = generate_equel(left, node_number);

        

        double *coeff = coeff_out(func, cosin, nodes, node_number, second);


        
        free(nodes);


        return coeff;
}
        

double func1(double x , double y) {
        return cos(x * 10 * M_PI) +  cos(y * 10 * M_PI);
}

double func2(double x , double y) {
	if(fabs(x - 1) < eps) {
                return 0;
        }

        if(fabs(x) < eps) {
                return 0;
        }

	if(fabs(y - 1) < eps) {
                return 0;
        }

        if(fabs(y) < eps) {
                return 0;
        }
        return exp(1/((2 * y - 1) * (2 * y - 1) - 1)) * exp(1/((2 * x - 1) * (2 * x - 1) - 1))  ;
}

double func3(double x, double y) {
        if(fabs(x - 1) < eps) {
                return 0;
        }

        return exp(1/(x * x - 1)) * exp(1/(y * y - 1));
}

double func4(double x, double y) {
        if(fabs(x - 1) < eps) {
                return 0;
        }

        if(fabs(x) < eps) {
                return 0;
        }

        return exp(1/((2 * x - 1) * (2 * x - 1) - 1)) ;
}

double func5(double x , double y) {
        return cos(x * 10 * M_PI) *  cos(y * 10 * M_PI);
}


double cosinus(double x) {
        return cos(x);
}

double *generate_equel(double left, int num) {
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

double Scalar(double  **func, function_pointer cosin, double h, int m, int n, int second) {
/*
	for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                        cout << func[i][j] << " ";
                }
                cout << endl;
        }
	cout << cos(0) << endl; */
//	cout << " ------------------------------------------- +++++---" << endl;
        double scalar = 0;

        scalar += func[0][second] * cos(0) / 2;

//	cout << " ------------------------------------------- " << endl;

        for(int k = 1; k < n; k++) {
                scalar += func[k][second] * cos(M_PI * k * h * m);
        }

//        cout << " m = " << m << " scalar = " << scalar << endl;
        return scalar * h;
}

double *coeff_out(double **func, function_pointer cosin, double *values, int values_number, int second) {

        double *coeff = (double *) calloc(values_number , sizeof(double));
        double scalar = 0;
        double h = 1/((double)values_number - 0.5);
//        cout << " h = " << h << endl;



        for(int m = 0; m < values_number ; m++) {
                scalar = Scalar(func, cosin, h, m, values_number, second);

                coeff[m] = scalar * 2;
//                cout << " coeff = " << coeff[m] << endl;
        }

//	cout << " ------------------------------------------- ++++" << endl;
       
        double sum = 0;
        for(int i = 1; i < values_number; i++) {
                sum += coeff[i];
        }
        coeff[0] = func[0][ second] - sum;
        cout << coeff[0] << endl;
//	cout << " +++++++++++++++++++++++++++++++++" << endl;
        return coeff;
}
