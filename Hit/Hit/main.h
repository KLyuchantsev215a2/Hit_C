#include <armadillo>
#include <vector>
#define pi 3.14
using namespace arma;
class particle {
public:
	mat F{ {1, 0}, {0,1} };
	mat L{ {0, 0}, {0,0} };
	mat SIG{ {0, 0}, {0,0} };
};

double _pow(double base, int exponent) {
	double result = 1.;
	for (int i = 0; i < exponent; ++i) {
		result *= base;
	}

	return result;
}

double norm_vector(std::vector<double> x) {
	return sqrt(x[1] * x[1] + x[2] * x[2]);
}

std::vector<double> sub_coordinat(std::vector<double> x, std::vector<double> y) {
	x[1] -= y[1];
	x[2] -= y[2];
	return x;
}
mat ComputeL(mat L, int i, int j, int beta, double m, std::vector<std::vector<double>>  v, double *rho, std::vector<double> nabla_W) {
	L(1, beta) -= (m / rho[j] * (v[j][0] - v[i][0])*nabla_W[beta]);
	L(2, beta) -= (m / rho[j] * (v[j][1] - v[i][1])*nabla_W[beta]);
	return L;
}
double ComputeW(int i, int j, std::vector<std::vector<double>> x, double m, int N, double *h) {
	std::vector<double> r = sub_coordinat(x[i], x[j]);
	double W_i = 0;
	double W_j = 0;
	double h_i = h[i];
	double h_j = h[j];

	double q_i = norm_vector(r) / h_i;
	double q_j = norm_vector(r) / h_j;

	double C_i = 1 / (pi*h_i*h_i);
	double C_j = 1 / (pi*h_j*h_j);

	if ((q_i >= 0) && (q_i <= 1))
		W_i = C_i * (15. / 7.)*(2. / 3. - q_i * q_i + 1. / 2. * q_i*q_i*q_i);

	if ((q_i > 1) && (q_i <= 2))
		W_i = C_i * (5. / 14.)*_pow((2. - q_i), 3);

	if ((q_j >= 0) && (q_j <= 1))
		W_j = C_j * (15. / 7.)*(2. / 3. - q_j * q_j + 1. / 2. * q_j*q_j*q_j);
	if ((q_j > 1) && (q_j <= 2))
		W_j = C_j * (5. / 14.)*_pow((2. - q_j), 3);

	return 0.5*(W_i + W_j);
}

double ComputeRho(int i, std::vector<std::vector<double>> x, double m, int N, double *h) {

	double ro = 0;

	for (int j = 1; j <= N; j++)
		ro += m * ComputeW(i, j, x, m, N, h);
	return  ro;
}
mat ComputeStress(mat F, int i, double mu, double k) {
	mat F_3 = F;
	F_3(3, 3) = 1;
	double J = det(F);
	mat E = { { 1,0,0 }, {0, 1, 0}, {0, 0, 1} };
	mat B = F_3 * trans(F_3);  //left Cauchy-Green tensor
	mat Biso = pow(J, -2. / 3.)*B;    // isochoric part of the left Cauchy - Green tensor
	mat devBiso = Biso - 1 / 3 * trace(Biso)*E;   // deviatoric part Biso
	mat stress = (2. * mu*devBiso + k / 10. * (pow(J, 5) - pow(J, -5))*E) / J;

	mat stress_2(2, 2);
	stress_2 = stress;
	return stress_2;
}
double Compute_nabla_W(int i, int j, std::vector<std::vector<double>> x, double *h, int beta) {
	std::vector<double> r = sub_coordinat(x[i], x[j]);
	std::vector<double>nabla_W_i(2);
	std::vector<double>nabla_W_j(2);

	nabla_W_i[0] = 0;
	nabla_W_i[1] = 0;
	nabla_W_j[0] = 0;
	nabla_W_j[1] = 0;

	double h_i = h[i];
	double h_j = h[j];


	double q_i = norm_vector(r) / h_i;
	double q_j = norm_vector(r) / h_j;
	double C_i = 1 / (pi*h_i*h_i);
	double C_j = 1 / (pi*h_j*h_j);

	if ((q_i > 0) && (q_i < 1))
		nabla_W_i[beta] = C_i * (15. / 7.)*(-2. * q_i + 3. / 2. * q_i*q_i)*(r[beta] / (h_i*norm_vector(r)));
	if ((q_i >= 1) && (q_i <= 2))
		nabla_W_i[beta] = C_i * (15. / 7.)*(-1. / 2.)*(2. - q_i) *(2. - q_i) * (r[beta] / (h_i*norm_vector(r)));

	if ((q_j > 0) && (q_j < 1))
		nabla_W_j[beta] = C_j * (15. / 7.)*(-2. * q_j + 3. / 2. * q_j*q_j)*(r[beta] / (h_j*norm_vector(r)));
	if ((q_j >= 1) && (q_j <= 2))
		nabla_W_j[beta] = C_j * (15. / 7.)*(-1. / 2.)*(2. - q_j) *(2. - q_j) * (r[beta] / (h_j*norm_vector(r)));

	return(0.5*(nabla_W_j[beta] + nabla_W_i[beta]));

}