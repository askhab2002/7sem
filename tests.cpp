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

//в командную строку вводиться код вида узлов( 0 - самому вводить узлы, 1 - равные отрезки, 2 - чебышевские),
// название файла, откуда считываются узлы(если вводим самому). 

double **interpolation_table(function_pointer func, double *nodes, int node_number);
void interpolation_table_out(double **table, int node_number);

double *jordan(double **table, int node_number);
double *deep_jordan(double **matrix,  double *b, int n);

double random_double_range(double left, double right);

void function_out(function_pointer func, double *values, int values_number);
void polynomial_out(double *jordan, int node_number, double *values, int values_number);
void lagrange_out(double **table, int node_number, double *values,int values_number);
void comparison_table_out(double **table, double *jordan, int node_number, double *values, int values_number);

double *generate_equel(double left, double right, int num);
double *generate_ch(double left, double right, int num);
double *generate_rand(double left, double right, int num);
double func(double x);

double polynomial(double *jordan, int node_number, double value);
double lagrange(double **table, int node_number, double value);


double func1(double x);
double func2(double x);
double func3(double x);

int comparison(const void *a, const void *b);

int main(int argc, char *argv[]) {
	if(argc != 3 && argc != 4) {
                cout << "    Вы ввели некорректную командную строку" << endl;
                return 0;
        }
        

        string nodes_type_(argv[1]);
	int nodes_type = stoi(nodes_type_);
	string func_type_(argv[2]);
	int func_type = stoi(func_type_);
	int node_number = 0;

	cout << " Введите левый и правый конец интерполяции" << endl;
	double left, right;
	cin >> left;
	cin >> right;
        
	double *nodes = NULL;

	if(argc == 4) {
                string nodes_input_(argv[3]);
		string file_name = nodes_input_;
		ifstream in(file_name);

		if(in.is_open()) {
			in >> node_number;
		
			for(int i = 0; !in.eof(); i++) {
				in >> nodes[i];
                        }

                }
                in.close();		
	}
        
	else {
		cout << " Введите число узлов" << endl;
		cin >> node_number;
                

		if(nodes_type == 0) {
		        for(int i = 0; i < node_number; i++) {
			       cin >> nodes[i];
		        }
	        }
        }

	int part = 0;
	cout << " Введите число разбиения промежутка для построения графика" << endl;
	cin >> part;

        
        if(nodes_type == 1) {
	        nodes = generate_equel(left, right, node_number);
        }

	if(nodes_type == 2) {
		nodes = generate_ch(left, right, node_number);
        }

	if(nodes_type == 3) {
		nodes = generate_rand(left, right, node_number);
	}
        
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

	double **table = interpolation_table(func, nodes, node_number);
	interpolation_table_out(table, node_number);

	double *polynomial = jordan(table, node_number);
	
        double *partition = generate_equel(left, right, part);


        function_out(func, partition, part);
	polynomial_out(polynomial, node_number, partition, part);
	lagrange_out(table, node_number, partition, part);
	comparison_table_out(table, polynomial, node_number, partition, part);

	free(nodes);
	free(table[0]);
	free(table[1]);
	free(table);
	free(partition);
	free(polynomial);

        return 0;
               	       

}

int comparison(const void *a, const void *b) {
	if(fabs(*(const double*)a - *(const double*)b) < eps) {
		return 0;
	}

	if(*(const double*)a - *(const double*)b > eps) {
		return 1;
	}

	return -1;
}

double func1(double x) {
	return cos(x);
}

double func3(double x) {
        return fabs(x);
}

double func2(double x) {
	return (1/(1 + 25 * x * x));
}

double **interpolation_table(function_pointer func, double *nodes, int node_number) {
	double **table = (double **)calloc( 2, sizeof(double *));
	table[0] = (double *)calloc( node_number, sizeof(double));
	table[1] = (double *)calloc( node_number, sizeof(double));

	for(int i = 0; i < node_number; i++) {
		table[0][i] = nodes[i];
		table[1][i] = (*func)(nodes[i]);
	}

	return table;
}

void interpolation_table_out(double **table, int node_number) {
	ofstream out;
        out.open("interpolation_nodes.txt", std::ofstream::out | std::ofstream::trunc);

//        cout << "num " << node_number << endl;	
        for(int i = 0; i < node_number; i++) {
		out << table[0][i] << " " << table[1][i] << endl;
	}

	out.close();
}

double *jordan(double **table, int node_number) {
	double **matrix = (double **)calloc( node_number, sizeof(double *));
	for(int i = 0; i < node_number; i++) {
		matrix[i] = (double *)calloc( node_number, sizeof(double));
		matrix[i][0] = 1;

		for(int j = 1; j < node_number; j++) {
			matrix[i][j] = matrix[i][j - 1] * table[0][i];
		}

		
	}

	double *b = (double *)calloc( node_number, sizeof(double));
	for(int i = 0; i < node_number; i++) {
		b[i] = table[1][i];
	}

	

	return deep_jordan(matrix,  b, node_number);
}

double *deep_jordan(double **a,  double *b, int n) {

	double d, s;
	double *x = (double *)calloc( n, sizeof(double));

	 for (int k = 0; k < n; k++) { 
		 for (int j = k + 1 ; j < n; j++) { 
			 d = a[j][k] / a[k][k]; 
			 for (int i = k ; i < n; i++) {
				 a[j][i] = a[j][i] - d * a[k][i]; 
			 }
			 b[j] = b[j] - d * b[k]; 
		 }
	 }
	 
	 for (int k = n - 1; k >= 0; k--) {
		 d = 0;
		 for (int j = k; j < n; j++) {
			 s = a[k][j] * x[j]; 
			 d = d + s; 
		 }
		 
		 x[k] = (b[k] - d) / a[k][k]; 
	 }
         
	 for(int i = 0; i < n; i++) {
		 free(a[i]);
         }
	 free(a);
	 free(b);

	 return x;

}

double random_double_range(double left, double right) {
//	double range = right - left;
//	double div = (double)RAND_MAX/range;

//	std::uniform_real_distribution<double> unif(left, right);
//        std::default_random_engine re;
//
//        srandom(time(NULL));
        
//	long ltime;


//        ltime = time (NULL);
//        int stime = (unsigned int) ltime/2;
//        srand(stime);
	double f = (double)rand() / RAND_MAX;
        return left + f * (right - left);
}

double *generate_equel(double left, double right, int num) {
	double step = (right - left)/(num - 1);
        double *nodes = (double *)calloc( num, sizeof(double) );
        
        if(num >= 2) {
		nodes[0] = left;
	}

        for(int i = 1; i < num; i++) {
               nodes[i] = nodes[i - 1] + step;
        }

        return nodes;
}

double *generate_ch(double left, double right, int num) {
//	double *nodes = (double *)calloc( num,  sizeof(* nodes)) ;
        double *nodes = new double[num];
	if(num >= 2) {
	        nodes[0] = left;
                nodes[num - 1] = right;
	}
        
//	cout << left << " " << right << " ";

        for(int i = 1; i < num ; i++) {
                nodes[i] = (right + left)/2 + cos((2 * (double)i - 1)/(2 * (double)num) * M_PI) * (right - left)/2;
//		cout << nodes[i]  << " ";
        }

	cout << endl;

	qsort(nodes, num, sizeof(double), comparison);

	return nodes;
}

double *generate_rand(double left, double right, int num) {
	double *nodes = (double *)calloc( num, sizeof(double));
        
	if(num >= 2) {
	        nodes[0] = left;
                nodes[num - 1] = right;
	}

        for(int i = 1; i < num - 1; i++) {
                nodes[i] = random_double_range(left, right);
        }
        
	qsort(nodes, num, sizeof(double), comparison);

	return nodes;
}

double polynomial(double *jordan, int node_number, double value) {
	double sum = 0;
	for(int i = 0; i < node_number; i++) {
		sum += (pow( value, i) ) * jordan[i];
	}

	return sum;
}

double lagrange(double **table, int size, double x) {

	double *x_values = table[0];
	double *y_values = table[1];

	double lagrange = 0;
	double pol;

	for (int i = 0; i < size; i++)
	{
		pol = 1;
		for (int j = 0; j < size; j++)
		{
			if (j == i) continue;
			pol *= (x - x_values[j])/(x_values[i] - x_values[j]);
		}
		lagrange += pol * y_values[i];
	}
	return lagrange;
}
	        	
void polynomial_out(double *jordan, int node_number, double *values, int  values_number) {
	fstream out;
	out.open("polynomial_out.txt", std::ofstream::out | std::ofstream::trunc);

	for(int i = 0; i < values_number; i++) {
		out << values[i] << " " << polynomial(jordan, node_number, values[i]) << endl;
	}

	out.close();
}

void lagrange_out(double **table, int node_number, double *values, int values_number) {
        fstream out;
        out.open("lagrange_out.txt", std::ofstream::out | std::ofstream::trunc);

        for(int i = 0; i < values_number; i++) {
                out << values[i] << " " << lagrange(table, node_number, values[i]) << endl;
        }

        out.close();
}

void comparison_table_out(double **table, double *jordan, int node_number, double *values, int values_number) {
	fstream out;
	out.open("comparison_table.txt", std::ofstream::out | std::ofstream::trunc);

	for(int i = 0; i < values_number; i++) {
		out << values[i] << " " << polynomial(jordan, node_number, values[i]) << " " << lagrange(table, node_number, values[i]) << endl;

	}

	out.close();
}
		
				
void function_out(function_pointer func, double *values, int values_number) {
	fstream out;
	out.open("function_out.txt", std::ofstream::out | std::ofstream::trunc);

	for(int i = 0; i < values_number; i++) {
		out << values[i] << " " << (*func)(values[i]) << endl;

	}

	out.close();
}



