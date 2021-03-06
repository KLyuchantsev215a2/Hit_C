#include"pch.h"
//#include<math.h>
#include"main.h"
//#include<iostream>
#include<stdio.h>
#include<conio.h>
#include <iostream>
#include <fstream>
#include <omp.h>
//#include "particle.h"
#include <time.h> 

int main()
{
	//input:
	//rho = initial density
		//H = height of a plate
		//L = plate width
		//N = number of particles approximating the body for now 9 + 9 = 18
		// v_0 = initial velocity
		//outpt:
	//dynamics of the body
		//mu, k the constants of the material for the generalized Hooke Law
		//E - elastic modulus for materials
	double rho_0 = 1;
	double v_0 = 1;
	double Time = 1;
	double mu = 20;
	double k = 50;
	double E = 9 * k*mu / (3 * k + mu);

	int sqn = 5;
	int S = sqn * sqn;
	int N = sqn * sqn;
	//S = H * L;

	double m = rho_0 * S / N;
	double k_h = 1.;
	double dt = 1;

	int fr = 10;


	particle *P;
	P = new particle[N];


	std::vector<double> nabla_W(2);


	for (int i = 0; i < N; i++) {
		//P[i].SIG = ComputeStress(P[i].F,i, mu, k);
	}

	double *rho;
	rho = new double[N];
	for (int i = 0; i < N; i++) {
		rho[i] = rho_0;
	}

	double *h;
	h = new double[N];
	for (int i = 0; i < N; i++) {
		h[i] = sqrt(m*k_h/ rho[i]);
	}

	std::vector<std::vector<double>> x(N);
	for (int i = 0; i < N; i++) {
		x[i].resize(2);
		x[i][0] = 0;
		x[i][1] = 1;
	}

	std::vector<std::vector<double>> v(N);
	for (int i = 0; i < N; i++) {
		v[i].resize(2);
		v[i][0] = 0;
		v[i][1] = 1;
	}

	for (int n = 0; n< int(Time / dt); n++) {
		
		//for (int i = 0; i < N; i++)
			//P[i].L.zeros();

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				for (int beta = 0; i < 2; i++) {
					nabla_W[beta] = Compute_nabla_W(i, j, x, h, beta);
					double PI = 0; //ComputeViscocity(x, v, rho, i, j, h, E);
					// v = ComputeVelocity(i, j, beta, dt, m, v, rho, SIG, PI, nabla_W);
				}
			}
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				for (int beta = 0; i < 2; i++) {
					nabla_W [beta] = Compute_nabla_W(i, j, x, h, beta);
					//P[i].L = ComputeL(P[i].L, i, j, beta, m, v, rho, nabla_W);
				}
			}
		}



		for (int i = 0; i < N; i++) {
			// F(1:2, 1 : 2, i) = F(1:2, 1 : 2, i) + dt * L(1:2, 1 : 2, i)*F(1:2, 1 : 2, i);
			//mat dtLL = dt *P[i].L;
			//P[i].F = expmat(dtLL)*P[i].F;
		}

		for (int i = 0; i < N; i++) {
			//P[i].SIG = ComputeStress(P[i].F,i, mu, k);
		}

		for (int i = 0; i < N; i++) {
			x[i][0] = x[i][0] + dt * v[i][0];
			x[i][1] = x[i][1]+ dt * v[i][1];
		}

		for (int i = 0; i < N; i++) {
			rho[i] = ComputeRho(i, x, m, N, h);
		}


		/*if (int(n / fr) == n / fr) {
			n = n;
			x_coord(1:N) = x(1, 1, 1:N);
			y_coord(1:N) = x(1, 2, 1:N);
			subplot(2, 2, 1);
			scatter(x_coord, y_coord);
			detF = ones(1, N);

			for (int i = 0; i < N; i++)
				detF(1, i) = det(F(1:2, 1 : 2, i));


			tri = delaunay(x_coord, y_coord);
			subplot(2, 2, 2);
			trisurf(tri, x_coord, y_coord, detF);

			subplot(2, 2, 3);
			trisurf(tri, x_coord, y_coord, rho(1:N));

			errSIG = ones(1, N);

			for (int i = 0; i < N; i++) {
				SIG15 = SIG(1:2, 1 : 2, 15);
				errSIG(1, i) = norm((SIG(1:2, 1 : 2, 15) - SIG(1:2, 1 : 2, i)));
			}
			subplot(2, 2, 4);
			trisurf(tri, x_coord, y_coord, errSIG);
		}*/
	}

}


