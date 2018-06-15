#include<vector>

class particle {
public:
	double ro;
	//double t_stop_s;

	std::vector<double>     x_coordinat;
	double              h_radius_kernel;
	double                            m;
	double                            V;

	//double v_s;

	//int					       sector;

	double                           W(particle x);
	//double  W_s(particle x);
	double	                W_gradient(particle x);
	
};

class gas_particle : public particle {
public:

	double                   P;
	double                   U;
	double                  cs;
	double                 v_s;
	double                   E;

	double                 e_s;
	double t_stop_s;
	double  ro_particle(gas_particle *X, int B, int i);
	//int							 sector;
	double  Viscosity(gas_particle X, double gamma);

	gas_particle* particle_initialization(std::vector<std::vector<double> > x, double V, double E, double m, int N, int N_gran, double gamma, double h_radius_kernel, int N_L, int N_R, int N_gran_L, int N_gran_R,double k) {
		
		gas_particle *X;
		X = new gas_particle[N + N_gran];

		double ro_L = 1.;
		double ro_R = 0.125;

		double P_L = 1.;
		double P_R = 0.1;

		double U_L = 2.5;
		double U_R = 2.;

		double cs_L = 1.18321595662;
		double cs_R = 1.5830052443;
		

		/*переделать*/
		for (int i = N_gran_L; i < N + N_gran_L; i++) {  //инициализация данных о частице
			X[i].x_coordinat = (x[i]);
			X[i].m = m;
			X[i].V = V;
			X[i].E = E;
			

			if (i<N_L + N_gran_L) {
				X[i].h_radius_kernel = h_radius_kernel;
				X[i].ro = ro_L;
				X[i].P = P_L;
				X[i].U = U_L;
				X[i].cs = cs_L;
			}

			else {

				X[i].h_radius_kernel = h_radius_kernel;
				X[i].ro = ro_R;
				X[i].P = P_R;
				X[i].U = U_R;
				X[i].cs = cs_R;
			}
			X[i].h_radius_kernel =m*k / X[i].ro;
		//	X[i].h_radius_kernel = 0.005;
		}

		for (int i = 0; i < N_gran_L; i++) {  //инициализация данных о частице
			X[i].x_coordinat = (x[N_gran_L-i-1]);
			X[i].m = m;
			X[i].V = V;
			
			X[i].E = E;


			X[i].h_radius_kernel = h_radius_kernel;
			X[i].ro = ro_L;
			X[i].P = P_L;
			X[i].U = U_L;
			X[i].cs = cs_L;
			X[i].h_radius_kernel = m* k / X[i].ro;
		//	X[i].h_radius_kernel = 0.005;
		}

		for (int i = N + N_gran_L; i < N_gran + N; i++) {
			X[i].x_coordinat = (x[i]);
			X[i].m = m;
			X[i].V = V;
			
			X[i].E = E;

			X[i].h_radius_kernel = h_radius_kernel;
			X[i].ro = ro_R;
			X[i].P = P_R;
			X[i].U = U_R;
			X[i].cs = cs_R;
			
			X[i].h_radius_kernel =m* k / X[i].ro;
		//	X[i].h_radius_kernel = 0.005;
		}
		return(X);
	}
};

class dust_particle : public particle {
public:
	double ro;
	//double t_stop_s;
	double v_s_dust;
	double ro_particle(dust_particle *P, int B, int i);
	dust_particle* particle_dust_initialization(std::vector<std::vector<double> > x, double V, double m, int N, int N_gran, double gamma, int N_L, int N_R, int N_gran_L, int N_gran_R, double k, double d_2_g,double h_radius_kernel) {
		dust_particle *X;
		X = new dust_particle[N + N_gran];

		double ro_L = 1.*d_2_g;
		double ro_R = 0.125*d_2_g;
		/*переделать*/
		for (int i = N_gran_L; i < N + N_gran_L; i++) {  //инициализация данных о частице

			X[i].x_coordinat = (x[i]);
			X[i].m = m;
			X[i].V = V;

			if (i<N_L + N_gran_L) {
				X[i].h_radius_kernel = h_radius_kernel;
				X[i].ro = ro_L;

			}
			else {
				X[i].h_radius_kernel = h_radius_kernel;
				X[i].ro = ro_R;
			}
			X[i].h_radius_kernel = m*k / X[i].ro;
			//X[i].h_radius_kernel = 0.005;
		}

		for (int i = 0; i < N_gran_L; i++) {  //инициализация данных о частице
			X[i].x_coordinat = (x[N_gran_L - i - 1]);
			X[i].m = m;
			X[i].V = V;
			X[i].h_radius_kernel = h_radius_kernel;
			X[i].ro = ro_L;
			X[i].h_radius_kernel = m* k / X[i].ro;
			//X[i].h_radius_kernel = 0.005;
		}

		for (int i = N + N_gran_L; i < N_gran + N; i++) {
			X[i].x_coordinat = (x[i]);
			X[i].m = m;
			X[i].V = V;
			X[i].h_radius_kernel = h_radius_kernel;
			X[i].ro = ro_R;
			X[i].h_radius_kernel = m* k / X[i].ro;
			//X[i].h_radius_kernel = 0.005;
		}
		return(X);
	}

};



