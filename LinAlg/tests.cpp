#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>

#define eps1 1e-5
#define eps2 1e-5


typedef double (* function_pointer)(double x);

double Lambda(int n, double p, int N);
double Lambda_k(int n, int k, int N, double h);
double *TriangularSystemP(int N, double *b, function_pointer cosin);

double Tau(double **A, int N, double *q);
double *Step(double **A, int N, double *x, double *b, double Tau);
double RealError(double **A, int N, double *x, double *b);
double Richardson(double *x, double **A, double *b, double Tau, int N, double eps, int mIter);

double **DoMatrix(int N, double p);
double *DoB(function_pointer f, int N);
double *DoP(int N);
double **DoMatrixP(int N, double *p);
double *DoF(function_pointer f, int N);

double *TestR(double  *b, int N, double eps, double *error, int mIter, double p);

double *StepB(double **A, double **inv_A, double *b, double *x_k, int N, double Tau);
double BSolver(double *x, double **A, double *b, double Tau, int N, double eps, int mIter);
double *TestB(double  *f, int N, double eps, double *error, int mIter);

double function1(double x);
double function2(double x);

double cosinus(double x);
double cosi(int m, int k, int N);
double *coeff_out(double *func, function_pointer cosin, int values_number);
double Scalar_(double *func, function_pointer cosin, double h, int m, int n);
double fourier(double *coeff, int num, double x);
double *error(double *Y, double *f, double p, int N);
double *error_B(double *Y, double *f, double h, int N);
double error_B_(double *x, double *b, double **A, int N);

void test_1(int N, function_pointer f, function_pointer cosin, double p);
void test_2(int N, function_pointer f, double p, int mIter);
void test_3(int N, function_pointer f, int mIter);

double **MatrixMult(double **A, double **B, int N);
double *MatrixVector(double **A, double *x, int N);
double **Inverse(int N);

using namespace std;

int main(void) {
         
	function_pointer f = function2;
        function_pointer cosin = cosinus;

        cout << "    Введите число корней" << endl;
        
	int N = 0;
	cin >> N;

	double p = 0;
	cout << "    Введите число p" << endl;
        cin >> p;

	cout << "---------------" << endl;
	cout << "  Точный метод  " << endl;

	test_1(N, f, cosin, p);


	cout << "  Введите максимальное число итераций" << endl;
	int mIter = 0;

	cin >> mIter;

//	cout << "   Метод Ричардсона   " << endl;

//	test_2(N, f, p, mIter);

	cout << endl;
	cout << "   Метод с обуславливателем  " << endl;
	test_3(N, f, mIter);

	return 0;

}

void test_3(int N, function_pointer f, int mIter) {

	double error = mIter;
	error = 0;
//	double eps = eps2;

	double h = 1/((double)N - 0.5);

	double *b = (double *)calloc(N, sizeof(double));

	for(int i = 0; i < N; i++) {
		b[i] = (*f)(i * h);
	}

	cout << "   Функция или массив?(0, 1) " << endl;
        int index = 0;
        cin >> index;
        if(index == 1) {
                for(int i = 0; i < N/2; i++) {
                        b[i] = i;
                }

                index = 0;
                for(int i = N - 1; i >= N/2; i--) {
                        b[i] = -index;
                        index++;
                }
        }

	double **inv_A = Inverse(N);
	double **A = DoMatrix(N, 1);
	double *p = DoP(N);
	double **A_ = DoMatrixP(N, p);
	free(p);

//	double *Y = TestB(b, N, eps, &error, mIter);

	double *x = (double *)calloc(N, sizeof(double));

	double q = 0;
        double tau = Tau(A, N, &q);
	double *Y = StepB(A, inv_A, b, x, N, tau);
        
	free(x);

	error = error_B_(Y, b, A, N);

	cout << "  error = " << error << endl;
/*
	double *errors = error_B(Y, b, h, N);

        double error_norm = 0;

        double b_norm = 0;
*/
  //      cout << "  errors: " << endl;
        for(int i = 0; i < N; i++) {
//                error_norm += errors[i] * errors[i];
//                b_norm += b[i] * b[i];

		free(inv_A[i]);
		free(A[i]);
		free(A_[i]);

        }

	free(inv_A);
	free(A);
	free(A_);

//        cout << " error_norm = " << sqrt(error_norm) << "  relative error = " << sqrt(error_norm)/sqrt(b_norm) << endl;

	free(Y);
	free(b);
//	free(errors);

	return;
}

double *StepB(double **A, double **inv_A, double *b, double *x, int N, double Tau) {

	double *values = MatrixVector(A, x, N);

	for(int i = 0; i < N; i++) {
		values[i] = b[i] - values[i];
	}

	double *y = MatrixVector(inv_A, values, N);


	for(int i = 0; i < N; i++) {
		values[i] = x[i] + Tau * y[i];
	}

	free(y);

	return values;
}

double **MatrixMult(double **A, double **B, int N) {

	double **C = (double **)malloc(N * sizeof(double *));

	for(int i = 0; i < N; i++) {
		C[i] = (double *)calloc(N, sizeof(double));
	}

	double sum = 0;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			for(int k = 0; k < N; k++) {
				sum += A[i][k] * B[k][j];
			}
			C[i][j] = sum;
			sum = 0;
		}
	}

	return C;
}

double *MatrixVector(double **A, double *x, int N) {

	double *y = (double *)calloc(N, sizeof(double));

	double sum = 0;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {

			sum += A[i][j] * x[j];
		}

		y[i] = sum;

		sum = 0;
	}

	return y;
}



double **Inverse(int N) {

	double **inv_A = (double **)malloc(N * sizeof(double *));
        double *f = (double *)calloc(N, sizeof(double));
        double *non = (double *)calloc(N, sizeof(double));

	function_pointer cosin = cosinus;

        for(int i = 0; i < N; i++) {
                inv_A[i] = (double *)calloc(N, sizeof(double));
                f[i] = 0;
        }



        for(int j = 0; j < N; j++) {
                f[j] = 1;
                if(j != 0) {
                        f[j - 1] = 0;
                }

                non = TriangularSystemP(N, f, cosin);
                for(int i = 0; i < N; i++) {
                        inv_A[i][j] = non[i];
//                      cout << inv_A[i][j] << " ";
                }
                free(non);
//                cout << endl;

        }

	free(f);


	return inv_A;
}
		


double BSolver(double *x, double **A, double *b, double Tau, int N, double eps, int mIter) {

        double **x_k = (double **)malloc(mIter * sizeof(double *));
        double precision = 0;

	double *err = (double *)calloc(N, sizeof(double));
	
        int end = mIter - 1;

	double **inv_A = Inverse(N);

        x_k[0] = (double *)calloc(N, sizeof(double));
        for(int i = 0; i < N; i++) {
                x_k[0][i] = x[i];
        }

        for(int i = 1; i < mIter; i++) {
		precision = 0;

                x_k[i] = StepB(A, inv_A, b, x_k[i - 1], N, Tau);
		
		err = error_B(x_k[i], b, 1/((double)N - 0.5), N);

		for(int j = 0; j < N; j++) {
			precision += err[j] * err[j];
			
		}

		precision = sqrt(precision);

                free(x_k[i - 1]);

                if(precision < eps) {
                        cout <<  "++++ " << endl;

                        end = i;
                        cout << " i = " << i << endl;
			
                        break;
                }

                free(err);

        }

	

	cout << " errorBB = " << precision;

        cout << " iterations = " << end << endl;

        for(int i = 0; i < N; i++) {
                x[i] = x_k[end][i];
        }

        free(x_k[end]);
        free(x_k);

	

	for(int i = 0; i < N; i++) {
		free(inv_A[i]);
	}
	free(inv_A);

        return precision;
}

double *TestB(double *b, int N, double eps, double *error, int mIter) {

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

        *error = BSolver(x, A, b, tau, N, eps, mIter);

        for(int i = 0; i < N; i++) {
                free(A[i]);
        }
        free(A);
        

	free(p);

        return x;
}

void test_2(int N, function_pointer f, double p, int mIter) {
	double error_ = 0;

  
  	double eps = eps2;

	double *b = DoB(f, N);

	cout << "   Функция или массив?(0, 1) " << endl;
        int index = 0;
        cin >> index;
        if(index == 1) {
                for(int i = 0; i < N/2; i++) {
                        b[i] = i;
                }

                index = 0;
                for(int i = N - 1; i >= N/2; i--) {
                        b[i] = -index;
                        index++;
                }
        }

        double *Y = TestR(b, N, eps, &error_, mIter, p);



        cout << " error = " << error_ << endl;


        double h = 1/((double)N - 0.5);


        double *errors = error(Y, b, p,  N);

        double error_norm = 0;

	double b_norm = 0;

        cout << "  errors: " << endl;
        for(int i = 0; i < N; i++) {
                error_norm += errors[i] * errors[i];
                b_norm += (*f)(i * h) * (*f)(i * h);

        }

        cout << " error_norm = " << sqrt(error_norm) << "  relative error = " << sqrt(error_norm)/sqrt(b_norm) << endl;

        free(Y);
	free(b);
        free(errors);

	return;
}

void test_1(int N, function_pointer f, function_pointer cosin, double p) {

	double h = 1/((double)N - 0.5);

        double *F = (double *)calloc(N, sizeof(double));

        for(int i = 0; i < N; i++) {
                F[i] = (*f)(i * h);
        }

	cout << "   Функция или массив?(0, 1) " << endl;
	int index = 0;
	cin >> index;
	if(index == 1) {
		for(int i = 0; i < N/2; i++) {
			F[i] = i;
		}

		index = 0;
		for(int i = N - 1; i >= N/2; i--) {
			F[i] = -index;
			index++;
		}
	}

	double *D = coeff_out(F, cosin, N); 

        double *C = (double *)calloc(N, sizeof(double));

        index = 0;
	for(int i = 0; i < N; i++) {
                if(i == 0 && fabs(Lambda(i, p, N)) < eps2) {
                       C[i] = 0;
		       index = 1;
                       continue;
                }
                C[i] = D[i]/Lambda(i, p, N);
        }


        double *Y = (double *)calloc(N, sizeof(double));
        
        for(int k = 0; k < N; k++) {
              Y[k] = fourier(C, N, k * h);
  
        }

	double *errors = error(Y, F, p, N);
        
	double error_norm = 0;
	double b_norm = 0;

        for(int i = 0; i < N; i++) {
                  error_norm += errors[i] * errors[i];
		  b_norm += F[i] * F[i];

        }

	cout << "error_norm = " << sqrt(error_norm) << "   relative error = " << sqrt(error_norm)/sqrt(b_norm) << endl;

        free(Y);
        free(D);
        free(C);
	free(errors);

	free(F);

	return ;

}



double *error(double *Y, double *f, double p, int N) {
	double *errors = (double *)calloc(N, sizeof(double));

	errors[0] = f[0] - p * Y[0] - 2 * (N - 0.5) * (N - 0.5) * (Y[0] - Y[1]);
	errors[N - 1] = f[N - 1] - p * Y[N - 1] - (N - 0.5) * (N - 0.5) * (Y[N - 1] - Y[N - 2]);

	for(int i = 1; i < N - 1; i++) {
		errors[i] = f[i] - p * Y[i] + (N - 0.5) * (N - 0.5) * (Y[i + 1] - 2 * Y[i] + Y[i - 1]);
	}

	return errors;
}

double *error_B(double *Y, double *f, double h, int N) {
	double *errors = (double *)calloc(N, sizeof(double));
        double *p_B = (double *)calloc(N, sizeof(double));

	for(int i = 0; i < N; i++) {
		p_B[i] = 1 + sin(M_PI * i * h) * sin(M_PI * i * h);
	}

        errors[0] = f[0] - p_B[0] * Y[0] - 2 * (N - 0.5) * (N - 0.5) * (Y[0] - Y[1]);
        errors[N - 1] = f[N - 1]  - p_B[N - 1] * Y[N - 1] - (N - 0.5) * (N - 0.5) * (Y[N - 1] - Y[N - 2]);

        for(int i = 1; i < N - 1; i++) {
                errors[i] = f[i] - p_B[i] * Y[i] + (N - 0.5) * (N - 0.5) * (Y[i + 1] - 2 * Y[i] + Y[i - 1]);
        }

	free(p_B);

        return errors;
}

double error_B_(double *x, double *b, double **A, int N) {
        double *y = MatrixVector(A, x, N);
        double norm = 0;
        
        for(int i = 0; i < N; i++) {
                 if(fabs(y[i] - b[i]) > norm) {
                         norm = fabs(y[i] - b[i]);
                 }
        }
        
        return norm;
}

double fourier(double *coeff, int num, double x) {

        double sum = 0;

        for(int i = 0; i < num ; i++) {
                sum += coeff[i] * cos(M_PI * i * x);
        }

        return sum;
}

double Lambda(int n, double p, int N) {
        return p + 2 * ((double)N - 0.5) * ((double)N - 0.5) * (1 - cos((M_PI * n)/((double)N - 0.5)));
}

double Lambda_k(int n, int k, int N, double h) {
	double p_k = 1 + sin(M_PI * k * h) * sin(M_PI * k * h);

        return p_k + 2 * ((double)N - 0.5) * ((double)N - 0.5) * (1 - cos((M_PI * n)/((double)N - 0.5)));
}

double cosinus(double x) {
	return cos(x);
}

double cosi(int m, int k, int N) {
	return cos((M_PI * m * k)/((double)N - 0.5));
}

double Scalar_(double *func, function_pointer cosin, double h, int m, int n) {
        double scalar = 0;

	scalar += func[0] * (*cosin)(0) / 2;

	for(int k = 1; k < n; k++) {
		scalar += func[k] * (*cosin)(M_PI * k * h * m);
	}

	return scalar * h;
}

double *coeff_out(double *func, function_pointer cosin, int values_number) {
        
	double *coeff = (double *) calloc(values_number , sizeof(double));
        double scalar = 0;
	double h = 1/((double)values_number - 0.5);

	for(int m = 0; m < values_number ; m++) {
		scalar = Scalar_(func, cosin, h, m, values_number);

		coeff[m] = scalar * 2;

	}

	double sum = 0;
        for(int i = 1; i < values_number; i++) {
	        sum += coeff[i];
	}
        coeff[0] = func[0] - sum;

        

	return coeff;
}

double *TestR(double  *b, int N, double eps, double *error, int mIter, double p) {

	double **A = DoMatrix(N, p);

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

		for(int j = i - 1; j < i + 2; j++) {
			if(j == -1 || j == N) {
				continue;
			}
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
	       for(int j = i - 1; j < i + 2; j++) {
                       if(j == -1 || j == N) {
                                continue;
                        }
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



	int end = mIter - 1;

	x_k[0] = (double *)calloc(N, sizeof(double));
	for(int i = 0; i < N; i++) {
		x_k[0][i] = x[i];
	}

	for(int i = 1; i < mIter; i++) {

		x_k[i] = Step(A, N, x_k[i - 1], b, Tau);
		precision = RealError(A, N, x_k[i], b);

		if(precision < eps) {
			cout <<  "++++ " << endl;
			
			end = i;
			cout << " i = " << i << endl;
			break;
		}

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

       matrix[0][0] = (N - 0.5) * (N - 0.5) * 2. + p;
       matrix[0][1] = -2. * (N - 0.5) * (N - 0.5);
       matrix[N - 1][N - 2] = -1. * (N - 0.5) * (N - 0.5);
       matrix[N - 1][N - 1] = 1. * (N - 0.5) * (N - 0.5) + p ;

       for(int i = 1; i < N - 1; i++) {
	       matrix[i] = (double *)calloc(N, sizeof(double));
	       matrix[i][i - 1]  = -1. * (N - 0.5) * (N - 0.5);
	       matrix[i][i] = 2. * (N - 0.5) * (N - 0.5) + p;
	       matrix[i][i + 1] = -1. * (N - 0.5) * (N - 0.5) ;
       }

       return matrix;

}

double *DoB(function_pointer f, int N) {

	double *b = (double *)calloc(N, sizeof(double));

	for(int i = 0; i < N; i++) {
		b[i] = (*f)((i)/((double)N - 0.5));
	}

	return b;
}

double *DoP(int N) {
	double *p = (double *)calloc(N, sizeof(double));

	for(int i = 0; i < N; i++) {
		p[i] = 1 + sin(M_PI * i/((double)N - 0.5)) * sin(M_PI * i/((double)N - 0.5));
	}

	return p;
}

double **DoMatrixP(int N, double *p) {
       double **matrix = (double **)malloc(N *  sizeof(double *));

       matrix[0] = (double *)calloc(N, sizeof(double));
       matrix[N - 1] = (double *)calloc(N, sizeof(double));

       matrix[0][0] = (N - 0.5) * (N - 0.5) * 2. + p[0];
       matrix[0][1] = -2. * (N - 0.5) * (N - 0.5);
       matrix[N - 1][N - 2] = -1. * (N - 0.5) * (N - 0.5);
       matrix[N - 1][N - 1] = (N - 0.5) * (N - 0.5) * 1. + p[N - 1];

       for(int i = 1; i < N - 1; i++) {
               matrix[i] = (double *)calloc(N, sizeof(double));
               matrix[i][i - 1]  = -1. * (N - 0.5) * (N - 0.5);
               matrix[i][i] = (N - 0.5) * (N - 0.5) * 2. + p[i];
               matrix[i][i + 1] = -1. * (N - 0.5) * (N - 0.5);
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




double *TriangularSystemP(int N, double *b, function_pointer cosin) {

	double h = 1/((double)N - 0.5);

        

        double *D = coeff_out(b, cosin, N); 

        double *C = (double *)calloc(N, sizeof(double));

        for(int i = 0; i < N; i++) {
                if(i == 0 && fabs(Lambda(i, 1, N)) < eps2) {
                       C[i] = 0;
                       continue;
                }
                C[i] = D[i]/Lambda(i, 1, N);
        }


        double *Y = (double *)calloc(N, sizeof(double));
        
        for(int k = 0; k < N; k++) {
                Y[k] = fourier(C, N, k * h);
        }

	free(D);
	free(C);

        return Y;
}


