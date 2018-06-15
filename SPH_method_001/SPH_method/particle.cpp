#include"particle.h"
#include<iostream>
#include <omp.h>
double _pow(double base, int exponent) {
	double result = 1.;
	for (int i = 0; i < exponent; ++i) {
		result *= base;
	}

	return result;
}

 double norm_vector(std::vector<double> x){
	double norm_vector_x = 0;
	
	for (int i = 0; i < fabs(x.size()); i++)
		norm_vector_x += x[i] * x[i];

	return sqrt(norm_vector_x);
}

std::vector<double> sub_coordinat(std::vector<double> x, std::vector<double> y){
	int size = x.size();

	if (x.size() > y.size())
		size = y.size();

	for (int i = 0; i < size; i++)
		x[i] = x[i] - y[i];

	return x;
}


double particle::W(particle P){
	std::vector<double> x_i = sub_coordinat(this->x_coordinat, P.x_coordinat);
	
	double h_a =                                        (this->h_radius_kernel);
	double h_b =											 P.h_radius_kernel;

	double q_a =                                        norm_vector(x_i) / h_a;
	double q_b =                                        norm_vector(x_i) / h_b;
	
	double C_a =                                                 2. / (3.*h_a);
	double C_b =                                                 2. / (3.*h_b);
	
	double W_a=0;
	double W_b=0;
	
	if ((q_a <= 1) && (q_a>= 0))
		W_a= C_a *       (1 - (3. / 2.)*   _pow(q_a, 2) + (3. / 4.) * _pow(q_a, 3));

	if ((q_a >1) && (q_a <= 2))
		W_a= C_a*             (1. / 4.)*                      (_pow(2 - q_a, 3));

	if ((q_b <= 1) && (q_b >= 0))
		W_b = C_b *       (1 - (3. / 2.)*   _pow(q_b, 2) + (3. / 4.) * _pow(q_b, 3));

	if ((q_b >1) && (q_b <= 2))
		W_b = C_b*             (1. / 4.)*                      (_pow(2 - q_b, 3));


	return 0.5*(W_a+W_b);
}



double particle::W_gradient(particle P){
	std::vector<double> x_i = sub_coordinat(this->x_coordinat, P.x_coordinat);

	double h_a =                                        (this->h_radius_kernel);
	double h_b =                                              P.h_radius_kernel;

	double pi =                                              3.14159265358979;

	double q_a =                                           norm_vector(x_i) / h_a;
	double q_b =                                           norm_vector(x_i) / h_b;
	double sign =                                 (x_i[0] > 0) - (x_i[0] < 0);

	double const_C_a =                                 (2. / (3.*h_a))*(sign / h_a);
	double const_C_b =                                 (2. / (3.*h_b))*(sign / h_b);

	double W_grad_a = 0;
	double W_grad_b = 0;

	if ((q_a <= 1) && (q_a >= 0))
		W_grad_a= const_C_a*(    -3. * q_a + (9. / 4.)*_pow(q_a, 2)   );

	if ((q_a > 1) && (q_a <= 2))
		W_grad_a= const_C_a*(   (-3. / 4.)*    _pow(2 - q_a, 2)     );
	
	if ((q_b <= 1) && (q_b >= 0))
		W_grad_b = const_C_b*(-3. * q_b + (9. / 4.)*_pow(q_b, 2));

	if ((q_b > 1) && (q_b <= 2))
		W_grad_b = const_C_b*((-3. / 4.)*    _pow(2 - q_b, 2));

	return 0.5*(W_grad_a+W_grad_b);
}

double gas_particle::ro_particle(gas_particle *P, int B, int i) {
	double ro = 0;	
	#pragma omp parallel for reduction(+:ro)
	for (int j = i - B; j <= i + B; j++)
			ro += this->m* (this->W(P[j]));
	return  ro;	
}
double dust_particle::ro_particle(dust_particle *P, int B, int i) {
	double ro = 0;
	#pragma omp parallel for reduction(+:ro)
	for (int j = i - B; j <= i + B; j++) 
		ro += this->m* (this->W(P[j]));
	return  ro;
}


double gas_particle::Viscosity(gas_particle X, double gamma) {
	double r_ab = ((this->x_coordinat[0]) - (X.x_coordinat[0]));
	/*узнать точные условия на вязкость*/  /*близкие частицы*/               /*только проблемные частицы*/
	if ((this->V - X.V)* r_ab>0)
		return 0;
	/*
	double a = 1.;   //никак в статье одномерный случай
	double b = 2.;

	double h_ab =         (this->h_radius_kernel + X.h_radius_kernel) / 2.;
	double V_ab =                                          (this->V)-(X.V);
	double ro_ab =                                      (X.ro+this->ro)/2.;
	double c_ab =                                      (X.cs+ this->cs)/2.;

	double nu =                                                   h_ab*0.1;
	double mu = (X.h_radius_kernel*V_ab*r_ab) / (_pow(r_ab, 2) + _pow(nu, 2));

	double viscocity=                  (-a*c_ab*mu + b*_pow(mu, 2)) / ro_ab ;//знак
	*/
	
	double a = 1.;
	double ro_ab = (X.ro + this->ro) / 2.;
	double V_ab = (this->V) - (X.V);

	double w_ab = V_ab* ((r_ab > 0) - (r_ab < 0));
	double v_sig = this->cs + X.cs - 3 * w_ab;

	double viscocity = (-a / 2.)*((v_sig*w_ab) / ro_ab);
	
	return viscocity;

}




