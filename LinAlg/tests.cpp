#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>

#define eps1 1e-5
#define eps2 1e-5


typedef double (* function_pointer)(double x);

double Scalar(double *f, int m, int N);
double Lambda(int n, double p, int N);
double C_Coeff(int m, int N, double *Coeff);
double Y_Coeff(int m, int N, double *f, double p);
double *TriangularSystem(int N, double *f, double p);
double *TriangularSystemP(int N, double *f, double *p);

double Tau(double **A, int N, double *q);
double *Step(double **A, int N, double *x, double *b, double Tau);
double RealError(double **A, int N, double *x, double *b);
double Richardson(double *x, double **A, double *b, double Tau, int N, double eps, int mIter);

double **DoMatrix(int N, double p);
double *DoB(function_pointer f, int N);
double *DoP(int N);
double **DoMatrixP(int N, double *p);
double *DoF(function_pointer f, int N);

double *TestR(function_pointer f, int N, double eps, double *error, int mIter, double p);

double *Step1(double **A, double *b, double *x_k, int N, double *p, double Tau);
double BSolver(double *x, double **A, double *b, double *p, double Tau, int N, double eps, int mIter);
double *TestB(double *f, int N, double eps, double *error, int mIter);

double function1(double x);
double function2(double x);

double cosinus(double x);
double cosi(int m, int k, int N);
double *coeff_out(function_pointer func, function_pointer cosin, int values_number);
double Scalar_(function_pointer func, function_pointer cosin, double h, int m, int n);

using namespace std;

int main(void) {
         
	function_pointer f = function2;
        function_pointer cosin = cosinus;

        cout << "    Введите число узлов" << endl;
        
	int N = 0;
	cin >> N;

/*	double *F = DoF(f, N);

	cout << " f = " << endl;
	for(int i = 0; i < N; i++) {
		cout << F[i] << endl;
	}
*/
	double p = 0;
//	cout << "    Введите число p" << endl;
//        cin >> p;

	cout << "---------------" << endl;

	double *Y = coeff_out(f, cosin, N + 1); // = TriangularSystem(N, F, p);

	cout << " Y : " << endl;
	for(int i = 0; i < N + 1; i++) {
		cout << Y[i] << endl;
	}

	cout << "    C: " << endl;
	double *C = (double *)calloc(N, sizeof(double));
	for(int i = 0; i < N; i++) {
                C[i] = Y[i]/Lambda(i + 1, p, N + 1);
		cout << C[i] << endl;
        }

	double *ans = (double *)calloc(N, sizeof(double));
        double sum = 0;

	cout << "   ans: " << endl;
	for(int k = 0; k < N; k++) {
		sum = 0;
		for(int m = 1; m < N + 1; m++) {
		       sum += cosi(m, k + 1, N + 1) * C[m - 1];
	        }
		ans[k] = sum;
		cout << ans[k] << endl;
        }		

	free(Y);
/*
	cout << "  Введите число итераций для Ричардсона" << endl;
	int mIter = 0;
	cin >> mIter;

	double error = 0;

	double eps = eps2;

	double *x = TestR(f, N, eps, &error, mIter, p);

	

        cout << " error = " << error << endl;

	cout << " x: " << endl;

        for(int i = 0; i < N; i++) {
		cout << x[i] << endl;
	}

        free(x);	

	double *z = TestB(F, N, eps, &error, mIter);

	cout << " error = " << error << endl;

        cout << " z: " << endl;

        for(int i = 0; i < N; i++) {
                cout << z[i] << endl;
        }

        free(z);
	free(F);
*/
	return 0;

}

double Lambda(int n, double p, int N) {
//      cout << "p  = " << p << endl;
        return p + 2 * ((double)N - 0.5) * ((double)N - 0.5) * (1 - cos((M_PI * n)/((double)N - 0.5)));
}

double cosinus(double x) {
	return cos(x);
}

double cosi(int m, int k, int N) {
	return cos((M_PI * m * k)/(N - 0.5));
}

double Scalar_(function_pointer func, function_pointer cosin, double h, int m, int n) {
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

double *coeff_out(function_pointer func, function_pointer cosin, int values_number) {
        
	cout << " ==================== " << endl;
	double *coeff = (double *) calloc(values_number , sizeof(double));
        double scalar = 0;
	double h = 1/((double)values_number - 0.5);
	cout << " h = " << h << endl;



	for(int m = 0; m < values_number ; m++) {
		scalar = Scalar_(func, cosin, h, m, values_number);

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

double *TestR(function_pointer f, int N, double eps, double *error, int mIter, double p) {

	double **A = DoMatrix(N, p);
	double *b = DoB(f, N);
 
	cout << " A: " << endl;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			cout << A[i][j] << " ";
		}
		cout << endl;
	}

	cout << " b: " << endl;
	for(int i = 0; i < N; i++) {
		cout << b[i] << endl;
	}

	double *x = (double *)calloc(N , sizeof(double));

	
	for(int i = 0; i < N; i++) {
		x[i] = 0;
	}

	double q = 0;
	double tau = Tau(A, N, &q);

	cout << " tau = " << tau << endl;
        
	*error = Richardson(x, A, b, tau, N, eps, mIter);

	for(int i = 0; i < N; i++) {
		free(A[i]);
	}
	free(A);
	free(b);

	return x;
}

double Tau(double **A, int N, double *q) {
	double m = 0;
	double M = 0;
	double sum = 0;

	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			sum += fabs(A[i][j]);
		}
                 
		if(sum > M) {
			M = sum;
		}

		if(m > A[i][i] + fabs(A[i][i]) - sum) {
			m = A[i][i] + fabs(A[i][i]) - sum;
		}

		sum = 0;
	}

	*q = (M - m)/(M + m);
        cout << " m = " << m << " M = " << M << endl;
	return 2/(m + M);
}

double *Step(double **A, int N, double *x, double *b, double Tau) {
	double *new_x = (double *)calloc(N, sizeof(double));
        double sum = 0;

	for(int i = 0; i < N; i++) {
		new_x[i] = x[i] + Tau * b[i];

		for(int j = 0; j < N; j++) {
			sum -= Tau * A[i][j] * x[j];
		}

		new_x[i] += sum;
		sum = 0;
	}

        free(x);
	return new_x;
}

double RealError(double **A, int N, double *x, double *b) {
       double *new_b = (double *)calloc(N, sizeof(double));
       double *bb = (double *)calloc(N, sizeof(double));
       for(int i = 0; i < N; i++) {
	       bb[i] = b[i];
       }

       double sum = 0;

       for(int i = 0; i < N; i++) {
	       for(int j = 0; j < N; j++) {
		       new_b[i] += A[i][j] * x[j];
	       }
       }

       for(int i = 0; i < N; i++) {
	       bb[i] -= new_b[i];
       }

       for(int i = 0; i < N; i++) {
	       sum += bb[i] * bb[i];
       }
       
       free(new_b);
       free(bb);

       return sqrt(sum);
}

double Richardson(double *x, double **A, double *b, double Tau, int N, double eps, int mIter) {

	double **x_k = (double **)malloc(mIter * sizeof(double *));
        double precision = 0;

//	cout << " ---------------------- " << endl;

	int end = mIter - 1;

	x_k[0] = (double *)calloc(N, sizeof(double));
	for(int i = 0; i < N; i++) {
		x_k[0][i] = x[i];
	}

	for(int i = 1; i < mIter; i++) {
//		cout << " ---------------------- " << endl;
		x_k[i] = Step(A, N, x_k[i - 1], b, Tau);
		precision = RealError(A, N, x_k[i], b);

//		free(x_k[i - 1]);
//		cout << "---------" << precision << endl;
		if(precision < eps) {
			cout <<  "++++ " << endl;
			
			end = i;
			cout << " i = " << i << endl;
			break;
		}
                

//		cout << i << " " ;

	}

	cout << " iterations = " << end << endl;

	for(int i = 0; i < N; i++) {
		x[i] = x_k[end][i];
	}

	free(x_k[end]);
	free(x_k);
         
	return precision;
}

double **DoMatrix(int N, double p) {
       double **matrix = (double **)malloc(N *  sizeof(double *));

       matrix[0] = (double *)calloc(N, sizeof(double));
       matrix[N - 1] = (double *)calloc(N, sizeof(double));

       matrix[0][0] = N * N * 2. + p;
       matrix[0][1] = -1. * N * N;
       matrix[N - 1][N - 2] = -1. * N * N;
       matrix[N - 1][N - 1] = 2. * N * N + p ;

       for(int i = 1; i < N - 1; i++) {
	       matrix[i] = (double *)calloc(N, sizeof(double));
	       matrix[i][i - 1]  = -1. * N * N;
	       matrix[i][i] = 2. * N * N + p;
	       matrix[i][i + 1] = -1. * N * N;
       }

       return matrix;

}

double *DoB(function_pointer f, int N) {

	double *b = (double *)calloc(N, sizeof(double));

//	cout << "b:" << endl;
	for(int i = 0; i < N; i++) {
		b[i] = (*f)((i + 1)/((double)N + 1));
//		cout << b[i] << endl;
	}

	return b;
}

double *DoP(int N) {
	double *p = (double *)calloc(N, sizeof(double));

	for(int i = 0; i < N; i++) {
		p[i] = 1 + sin(M_PI * i/((double)N + 1));
	}

	return p;
}

double **DoMatrixP(int N, double *p) {
       double **matrix = (double **)malloc(N *  sizeof(double *));

       matrix[0] = (double *)calloc(N, sizeof(double));
       matrix[N - 1] = (double *)calloc(N, sizeof(double));

       matrix[0][0] = 2. + p[0];
       matrix[0][1] = -1.;
       matrix[N - 1][N - 2] = -1.;
       matrix[N - 1][N - 1] = 2. + p[N - 1];

       for(int i = 1; i < N - 1; i++) {
               matrix[i] = (double *)calloc(N, sizeof(double));
               matrix[i][i - 1]  = -1.;
               matrix[i][i] = 2. + p[i];
               matrix[i][i + 1] = -1.;
       }

       return matrix;

}

double *DoF(function_pointer f, int N) {

	double *F = (double *)calloc(N, sizeof(double));

	for(int i = 0; i < N + 1; i++) {
		cout << i + 1 << " " << (i)/((double)N + 0.5) << endl;
		F[i] = (*f)((i)/((double)N + 0.5));
	}

	return F;
}


double function1(double x) {

	return cos(M_PI * x );
}

double function2(double x) {
        if(fabs(x - 1) < eps1) {
                return 0;
        }

        if(fabs(x) < eps1) {
                return 0;
        }

        return exp(1/((2 * x - 1) * (2 * x - 1) - 1));
}




double Scalar(double *f, int m, int N) {


	double sum = f[0] * (0.5/((double)N + 0.5));
   //     double sum1 = 0;

	for(int i = 1; i < N + 1; i++) {
		sum += (cos((M_PI * m * i)/((double)N + 0.5)) * f[i])/((double)N + 0.5);
//		sum1 += sin((M_PI * m * i)/((double)N + 1)) * sin((M_PI * m * i)/((double)N + 1))/((double)N + 1);
	}
        cout << m << "-----" << sum << endl;
	return sum;
}



double C_Coeff(int m, int N, double *f, double p) {
	
	return 2 * Scalar(f, m, N)/Lambda(m, p, N);
}

double Y_Coeff(int m, int N, double *Coeff) {
        double sum = 0;

	for(int i = 1; i < N + 1; i++) {
		sum += Coeff[i - 1] * cos((M_PI * m * i)/((double)N + 0.5));
	}

	return sum;
}

double *TriangularSystem(int N, double *f, double p) {

	double *Y = (double *)calloc(N, sizeof(double));
        double *Coeff = (double *)calloc(N, sizeof(double));
	for(int i = 0; i < N; i++) {
		Coeff[i] = C_Coeff(i + 1, N, f, p);
//		cout << C[i] << " ";
	}
//	cout << endl;

	for(int i = 1; i < N + 1; i++) {
		Y[i - 1] = Y_Coeff(i, N, Coeff);
	}

	free(Coeff);
	return Y;
}

double *TriangularSystemP(int N, double *f, double *p) {

        double *Y = (double *)calloc(N, sizeof(double));
	double *Coeff = (double *)calloc(N, sizeof(double));
        for(int i = 0; i < N; i++) {
                Coeff[i] = C_Coeff(i + 1, N, f, p[i]);
        }

        for(int i = 1; i < N + 1; i++) {
                Y[i - 1] = Y_Coeff(i, N, Coeff);
        }

	free(Coeff);
        return Y;
}

double *Step1(double **A, double *b, double *x_k, int N, double *p, double Tau) {
	double *values = (double *)calloc(N, sizeof(double));

	double sum = 0;

	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			sum += A[i][j] * x_k[j];
		}

		values[i] = b[i] - sum;
		sum = 0;
	}

	double *y_k = TriangularSystemP(N, values, p);

	for(int i = 0; i < N; i++) {
		values[i] = x_k[i] + Tau * y_k[i];
	}

	free(y_k);

	return values;
}


double BSolver(double *x, double **A, double *b, double *p, double Tau, int N, double eps, int mIter) {

        double **x_k = (double **)malloc(mIter * sizeof(double *));
        double precision = 0;

//      cout << " ---------------------- " << endl;

        int end = mIter - 1;

        x_k[0] = (double *)calloc(N, sizeof(double));
        for(int i = 0; i < N; i++) {
                x_k[0][i] = x[i];
        }

        for(int i = 1; i < mIter; i++) {
//              cout << " ---------------------- " << endl;
                x_k[i] = Step1(A, b, x_k[i - 1], N, p, Tau);
                precision = RealError(A, N, x_k[i], b);

                free(x_k[i - 1]);
//              cout << "---------" << precision << endl;
                if(precision < eps) {
                        cout <<  "++++ " << endl;

                        end = i;
                        cout << " i = " << i << endl;
                        break;
                }

//              cout << i << " " ;

        }

        cout << " iterations = " << end << endl;

        for(int i = 0; i < N; i++) {
                x[i] = x_k[end][i];
        }

        free(x_k[end]);
        free(x_k);

        return precision;
}

double *TestB(double *f, int N, double eps, double *error, int mIter) {

	double *p = DoP(N);
	double **A = DoMatrixP(N, p);

	double *x = (double *)calloc(N, sizeof(double));
//        cout << "=============" << endl;

        for(int i = 0; i < N; i++) {
                x[i] = 0;
        }

        double q = 0;
        double tau = Tau(A, N, &q);

        cout << " tau = " << tau << endl;

        *error = BSolver(x, A, f, p, tau, N, eps, mIter);

        for(int i = 0; i < N; i++) {
                free(A[i]);
        }
        free(A);
        

	free(p);

        return x;
}
