#include <math.h>
#include <iostream>
#include "FunctionPointer.hpp"
#include "Triangulation.hpp"

void Triangulation(double Lx, double Ly, int N, FILE *file, double ***triangles, double **nodes) {
	double x = Lx/N;
	double y = Ly/N;

	fprintf(file, "%d\n", (N + 1) * (N + 1));
	fprintf(file, "%d\n", 2 * N * N);
	fprintf(file, "%d\n", 3 * N * N - 2 * N);
	fprintf(file, "%d\n", 4 * N);

	

	for(int i = 0; i < (N + 1) * (N + 1); i++) {
		nodes[i] = (double *)calloc(2, sizeof(double));
		fprintf(file, "%d: %lf %lf\n", i + 1, x * (i  % (N + 1)), y * (i/(N + 1)));
		nodes[i][0] = x * (i  % (N + 1));
		nodes[i][1] = y * (i/(N + 1));
	}

	int k = 0;
	int corner = 0;




	for(int i = 1; i < (N + 1) ; i++) {
		for(int j = 1; j < N + 1; j++) {
//		        triangles[2 * ((i - 1) * N + (j - 1))] = (double **)calloc(3, sizeof(double));
//			triangles[2 * ((i - 1) * N + (j - 1)) + 1] = (double **)calloc(3, sizeof(double));
		        k++;
			corner = j + (i - 1) * (N + 1);
		        fprintf(file, "%d: %d %d %d\n", k, corner, corner + 1, corner + (N + 1));
		        k++;
		        fprintf(file, "%d: %d %d %d\n", k, corner + 1, corner + (N + 1), corner + (N + 1) + 1);
/*			triangles[2 * ((i - 1) * N + (j - 1))][0] = nodes[corner];
			triangles[2 * ((i - 1) * N + (j - 1))][1] = nodes[corner + 1];
			triangles[2 * ((i - 1) * N + (j - 1))][2] = nodes[corner + (N + 1)];
			triangles[2 * ((i - 1) * N + (j - 1)) + 1][0] = nodes[corner + 1];
			triangles[2 * ((i - 1) * N + (j - 1)) + 1][1] = nodes[corner + (N + 1)];
			triangles[2 * ((i - 1) * N + (j - 1)) + 1][2] = nodes[corner + (N + 1) + 1]; */
		}
	}

        for(int i = 0; i < N * N; i++) {
		triangles[2 * i] = (double **) calloc(3, sizeof(double));
		triangles[2 * i][0] = nodes[i + 1];
		triangles[2 * i][1] = nodes[i + 2];
		triangles[2 * i][2] = nodes[i + 2 + (N + 1)];

		triangles[2 * i + 1] = (double **) calloc(3, sizeof(double));
                triangles[2 * i + 1][0] = nodes[i + 1];
                triangles[2 * i + 1][1] = nodes[i + 1 + (N + 1)];
                triangles[2 * i + 1][2] = nodes[i + 2 + (N + 1)];
	}	

	k = 1;
	corner = 1; 

	fprintf(file, "%d: %d %d\n", k, corner, corner + 1 + (N + 1));

	k++;
	fprintf(file, "%d: %d %d\n", k, corner + 1, corner + 1 + (N + 1));

	corner = N;
	k++;
	fprintf(file, "%d: %d %d\n", k, corner, corner + 1 + (N + 1));

	for(int i = 2, j = 1; i < N + 1; i++) {
		k++;
		corner = j + (i - 1) * (N + 1);
		fprintf(file, "%d: %d %d\n", k, corner, corner + 1);
		k++;
		fprintf(file, "%d: %d %d\n", k, corner, corner + 1 + (N + 1));
		k++;
		fprintf(file, "%d: %d %d\n", k, corner + 1, corner + 1 + (N + 1));
	}

	for(int i = 2, j = N; i < N + 1; i++) {
                k++;
                corner = j + (i - 1) * (N + 1);
                fprintf(file, "%d: %d %d\n", k, corner, corner + 1);
                k++;
                fprintf(file, "%d: %d %d\n", k, corner, corner + 1 + (N + 1));
        }

	for(int i = 1, j = 2; j < N; j++) {
                k++;
                corner = j + (i - 1) * (N + 1);
                fprintf(file, "%d: %d %d\n", k, corner + 1, corner + 1 + (N + 1));
                k++;
                fprintf(file, "%d: %d %d\n", k, corner, corner + 1 + (N + 1));
        }


	for(int i = 2; i < N + 1; i++) {
		for(int j = 2; j < N; j++) {
			k++;
			fprintf(file, "%d: %d %d\n", k, corner, corner + 1 + (N + 1));
			k++;
			fprintf(file, "%d: %d %d\n", k, corner, corner + 1 + (N + 1));
			k++;
			fprintf(file, "%d: %d %d\n", k, corner, corner + 1 + (N + 1));
		}
	}

	k = 1;

	for(int i = 1, j = 1; j < N + 1; j++) {
	       corner = j + (i - 1) * (N + 1);
               fprintf(file, "%d: %d %d\n", k, corner, corner + 1);
	       k++;      
        }

	for(int i = N + 1, j = 1; j < N + 1; j++) {
		corner = j + (i - 1) * (N + 1);
               fprintf(file, "%d: %d %d\n", k, corner, corner + 1);
               k++;
        }

	for(int i = 1, j = 1; i < N + 1; i++) {
               corner = j + (i - 1) * (N + 1);
               fprintf(file, "%d: %d %d\n", k, corner, corner + (N + 1));
               k++;
        }

	for(int i = 1, j = N + 1; i < N + 1; i++) {
               corner = j + (i - 1) * (N + 1);
               fprintf(file, "%d: %d %d\n", k, corner, corner + (N + 1));
               k++;
        }
/*
	for(int i = 0; i < (N + 1) * (N + 1); i++) {
		free(nodes[i]);
	}

	free(nodes);

*/
}	

			








	
