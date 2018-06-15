/*SPH метод, задача разрыва (одномерная) */
#include<stdio.h>
#include<conio.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include "particle.h"
#include <time.h> 
/*
std::vector<std::vector<double> > start_coordinat(int N,int N_R,int N_L, int N_gran, int N_gran_L) {
	std::vector<std::vector<double> > x(N + N_gran);
	/*инициализация коррдинат газа*/ /*переделать
	for (int i = N_gran_L; i < N + N_gran_L; i++) {
		x[i].resize(3);

		if (i<N_L + N_gran_L)
			x[i][0] = (0.5 / double(N_L))*(i - N_gran_L);
		else
			x[i][0] = 0.5 + (0.5 / double(N_R))*(i - N_L - N_gran_L + 1);

		x[i][1] = 0;
		x[i][2] = 0;
	}

	for (int i = 0; i <N_gran_L; i++) {
		x[i].resize(3);

		x[i][0] = -(0.5 / double(N_L))*(i + 1);

		x[i][1] = 0;
		x[i][2] = 0;
	}
	for (int i = N + N_gran_L; i < N_gran + N; i++) {
		x[i].resize(3);

		x[i][0] = 1 + 0.5 / double(N_R)*(i - (N + N_gran_L) + 1);

		x[i][1] = 0;
		x[i][2] = 0;
	}
	return x;
}
*/

void main() {                 //сделать адекватно если будет время
	clock_t elapsedTime;
	elapsedTime = clock();
	/*аппроксимация найти оптимальные параметры*/		
	/*найти аналитически!!!*/
	double CFL =                                0.1;//курант фридрихс леви
	double k =                                  1.4;//константа адаптиного радиуса сглаживания 
	double eps =                                 1.;
	
	double K =    2000.;
	double tau = 0.00001;//h_radius_kernel/3.34664 * CFL; // условие куранта
	double d_2_g = 1.;

	/*Неварьируемые константы*/
	int N_sector =                           500;
	int N =                                 999*2;
	int N_gran =                             99*4;
	int N_dust =                            999*2;
	int N_gran_dust =                        99*4;

	double gamma =                      7.0 / 5.0;
	double V =                                  0;
	double E =                                  0;
	double m =              1.125*0.5 / double(N);
	double m_dust =         d_2_g*1.125*0.5 / double(N_dust);
	double h_radius_kernel=   m*k;

	double t_max = 0.2;
	int T = t_max / tau;

	//начальное распределение газа
	int N_L =                           N / 9 * 8;
	int N_R =                                 N/9 ;
	int N_gran_L =                 N_gran / 9 * 8;
	int N_gran_R =                     N_gran / 9;

	//начальное распределение пыли
	int N_L_dust = N_dust / 9 * 8;
	int N_R_dust = N_dust / 9;
	int N_gran_L_dust = N_gran_dust / 9 * 8;
	int N_gran_R_dust = N_gran_dust / 9;

	//int N_sector = (1 + 4 * h_radius_kernel) / (2 * h_radius_kernel);
	std::vector<std::vector<double> > x(N + N_gran);   //координаты каждой частицы газа
	std::vector<std::vector<double> > x_dust(N_dust + N_gran_dust);   //координаты каждой частицы пыли
	
																	  /*инициализация коррдинат газа*/ /*переделать*/
	for (int i = N_gran_L; i < N + N_gran_L; i++) {
		x[i].resize(3);

		if (i<N_L + N_gran_L)
			x[i][0] = (0.5 / double(N_L))*(i - N_gran_L);
		else
			x[i][0] = 0.5 + (0.5 / double(N_R))*(i - N_L - N_gran_L + 1);

		x[i][1] = 0;
		x[i][2] = 0;
	}

	for (int i = 0; i <N_gran_L; i++) {
		x[i].resize(3);

		x[i][0] = -(0.5 / double(N_L))*(i + 1);

		x[i][1] = 0;
		x[i][2] = 0;
	}
	for (int i = N + N_gran_L; i < N_gran + N; i++) {
		x[i].resize(3);

		x[i][0] = 1 + 0.5 / double(N_R)*(i - (N + N_gran_L) + 1);

		x[i][1] = 0;
		x[i][2] = 0;
	}

	/*начальное распределение пыли*/
	for (int i = N_gran_L_dust; i < N_dust + N_gran_L_dust; i++) {
		x_dust[i].resize(3);

		if (i<N_L_dust + N_gran_L_dust)
			x_dust[i][0] = (0.5 / double(N_L_dust))*(i - N_gran_L_dust);
		else
			x_dust[i][0] = 0.5 + (0.5 / double(N_R_dust))*(i - N_L_dust - N_gran_L_dust + 1);

		x_dust[i][1] = 0;
		x_dust[i][2] = 0;
	}

	for (int i = 0; i <N_gran_L_dust; i++) {
		x_dust[i].resize(3);

		x_dust[i][0] = -(0.5 / double(N_L_dust))*(i + 1);

		x_dust[i][1] = 0;
		x_dust[i][2] = 0;
	}
	for (int i = N_dust + N_gran_L_dust; i < N_gran_dust + N_dust; i++) {
		x_dust[i].resize(3);

		x_dust[i][0] = 1 + 0.5 / double(N_R_dust)*(i - (N_dust + N_gran_L_dust) + 1);

		x_dust[i][1] = 0;
		x_dust[i][2] = 0;
	}

	gas_particle tmp_particle_initialization;
	dust_particle tmp_particle_dust_inicialization;
	
	gas_particle *X; //вектор газа
	dust_particle *X_dust;//вектор пыли

	double *V_tmp;
	double *V_tmp_dust;

	/*
	double *v_s;
	double *v_s_dust;

	double *t_stop_s;

	double *v_s_cout;
	double *v_s_dust_cout;

	
	v_s = new double[N_sector];
	v_s_dust = new double[N_sector];

	t_stop_s = new double[N_sector];

	v_s_cout = new double[N_sector];
	v_s_dust_cout = new double[N_sector];
	*/

	V_tmp = new double[N + N_gran_L];
	V_tmp_dust = new double[N_dust + N_gran_L_dust];
	
	for (int i = N_gran_L; i < N+N_gran_L; i++) {
		V_tmp[i] = 0;
		V_tmp_dust[i] = 0;
	}

	std::ofstream output_ro;
	output_ro.open("output_ro.txt");

	std::ofstream output_ro1;
	output_ro1.open("output_ro1.txt");

	std::ofstream output_t_max1;
	std::ofstream output_t_max2;
	std::ofstream output_t_0_1;
	std::ofstream t;
	output_t_0_1.open("output_t_0_1.txt");
	output_t_max1.open("output_t_max1.txt");
	output_t_max2.open("output_t_max2.txt");
	t.open("t.txt");

	/*инициализация частиц*/
	X = tmp_particle_initialization.particle_initialization(x, V, E, m, N,N_gran,gamma,h_radius_kernel,N_L,N_R,N_gran_L,N_gran_R,k);
	X_dust = tmp_particle_dust_inicialization.particle_dust_initialization(x_dust, V, m_dust, N_dust, N_gran_dust, gamma, N_L_dust, N_R_dust, N_gran_L_dust, N_gran_R_dust, k,d_2_g,h_radius_kernel);

	int B = 2* X[N_L+2+N_gran_L].h_radius_kernel*N_L / 0.5+1; //максимальное количество нужных частиц
	int B_dust =  2*X_dust[N_L+2+N_gran_L].h_radius_kernel*N_L_dust / 0.5 + 1;
	//B = 1;
	//B_dust = 1;
	std::cout << B <<" "<< B_dust<< "\n";
	
	int B_u_v_s = 2;  //усрелняем кол-во соседей
	int B_u_v_s_d = 0;
	
	/*метод SPH*/
	double Viscosity = 0;
	double gradient=0;
	int cout_swap = 0;

for (int j = 0; j < T; j++) {

		if (j % 20 == 1) {
			std::cout << cout_swap << " " << X_dust[2 + N_gran_L].ro << " " << X_dust[1 + N_gran_L].x_coordinat[0] << " " << X[1 + N_gran_L].h_radius_kernel << " " << X_dust[1 + N_gran_L].V << " " << 100.0*j*tau / t_max << "%" << "\n";
			std::cout << X[N_L + N_gran_L + 1].ro << " " << X[N_L + N_gran_L + 1].x_coordinat[0] << " " << X[N_L + 1 + N_gran_L].h_radius_kernel << " " << X[N_L + 1 + N_gran_L].V << " " << 100.0*j*tau / t_max << "%" << "\n";
		}

		if (j == 1) {
			for (int i = 0; i < N + N_gran; i++)
				output_ro << X[i].x_coordinat[0] << " " << X[i].ro << "\n";
			for (int i = 0; i < N_dust + N_gran_dust; i++)
				output_ro1 << X_dust[i].x_coordinat[0] << " " << X_dust[i].ro << "\n";
		}

		/*
		double min = 1;
		int b_tmp = 0;

		for (int a = N_gran_L; a < N + N_gran_L; a++) {
			for (int b = N_gran_L_dust; b < N_dust + N_gran_L_dust; b++) {

				if (min >= fabs((X[a].x_coordinat[0] - X_dust[b].x_coordinat[0])))
				{

					min = fabs((X[a].x_coordinat[0] - X_dust[b].x_coordinat[0]));
					b_tmp = b;

				}
			}
			min = 1;
			X[a].v_s = 0;

			for (int b = b_tmp - B_u_v_s; b <= b_tmp + B_u_v_s; b++) {
				//X[a].v_s+=m* V_tmp_dust[b]* X_dust[b].W_s(X[a])/X[b].ro;
				X[a].v_s += V_tmp_dust[b];
			}
			X[a].v_s = X[a].v_s / (2 * (B_u_v_s)+1);
			//X[a].e_s = 1;

		}


		double min_1 = 1;
		int a_tmp = 0;
		for (int b = N_gran_L_dust; b < N_dust + N_gran_L_dust; b++) {

			for (int a = N_gran_L; a < N + N_gran_L; a++) {
				if (min_1 >= fabs((X[a].x_coordinat[0] - X_dust[b].x_coordinat[0])))
				{

					min_1 = fabs((X[a].x_coordinat[0] - X_dust[b].x_coordinat[0]));
					a_tmp = a;

				}

			}
			min_1 = 1;
			X_dust[b].v_s_dust = 0;
			for (int a = a_tmp - B_u_v_s; a <= a_tmp + B_u_v_s; a++) {
				//X_dust[b].v_s_dust +=X_dust[b].m* V_tmp[a]*X_dust[b].W_s(X[a])/X_dust[a].ro;
				X_dust[b].v_s_dust += V_tmp[a];
			}
			X_dust[b].v_s_dust = X_dust[b].v_s_dust / (2 * (B_u_v_s)+1);
		}
		/*
		double min = 1;
		int b_tmp = 0;

		for (int a = N_gran_L; a < N + N_gran_L; a++) {
			for (int b = N_gran_L_dust; b <  N_dust + N_gran_L_dust; b++) {

				if (min > fabs((X[a].x_coordinat[0] - X_dust[b].x_coordinat[0])))
				{

					min = fabs((X[a].x_coordinat[0] - X_dust[b].x_coordinat[0]));
					b_tmp = b;

				}
			}

			X[a].v_s = V_tmp_dust[b_tmp] ;


			min = 1;

		}*/

		/*Расчеты средних скоростей пыли и газа*/
		for (int a = N_gran_L; a < N + N_gran_L; a++) {
			X[a].v_s = 0;
			for (int b = a - B_u_v_s_d; b <= a + B_u_v_s_d; b++) {
				//X[a].v_s+=m* V_tmp_dust[b]* X_dust[b].W_s(X[a])/X[b].ro;
				X[a].v_s += V_tmp_dust[b];
			}
			X[a].v_s = X[a].v_s / (2 * (B_u_v_s_d)+1);
		}
		//X[a].e_s = 1;

		for (int b = N_gran_L_dust; b < N_dust + N_gran_L_dust; b++) {
			X_dust[b].v_s_dust = 0;
			for (int a = b - B_u_v_s; a <= b + B_u_v_s; a++) {
				//X_dust[b].v_s_dust +=X_dust[b].m* V_tmp[a]*X_dust[b].W_s(X[a])/X_dust[a].ro;
				X_dust[b].v_s_dust += V_tmp[a];
				//	X_dust[b].t_stop_s += X_dust[a].ro;
			}
			//X_dust[b].t_stop_s = X_dust[b].t_stop_s / (2 * (B_u_v_s));
			//X[b].t_stop_s = X_dust[b].t_stop_s;
			X_dust[b].v_s_dust = X_dust[b].v_s_dust / (2 * (B_u_v_s)+1);
		}

		/*
		for (int a = 0; a < N_sector; a++) {
			v_s[a] = 0;
			v_s_cout[a] = 0;
			v_s_dust[a] = 0;
			v_s_dust_cout[a]=0;
			t_stop_s[a] = 0;
		}
		for (int a = N_gran_L; a < N + N_gran_L; a++) {
			X[a].sector = X[a].x_coordinat[0] * N_sector ;
			if (X[a].sector >=(N_sector))
				X[a].sector = N_sector;
			if (X[a].sector <= 0)
				X[a].sector = 0;
			v_s[X[a].sector] += X[a].V;
			v_s_cout[X[a].sector]++;

		}

		for (int a = N_gran_L_dust; a < N_dust + N_gran_L_dust; a++) {
			X_dust[a].sector = X_dust[a].x_coordinat[0] * N_sector;
			if (X_dust[a].sector >= (N_sector))
				X_dust[a].sector = N_sector;
			if (X_dust[a].sector <= 0)
				X_dust[a].sector = 0;
			v_s_dust[X_dust[a].sector] += X_dust[a].V;
			v_s_dust_cout[X_dust[a].sector]++;

			t_stop_s[X_dust[a].sector] += X_dust[a].ro;

		}
		for (int a = N_gran_L; a < N + N_gran_L; a++) {
			X[a].v_s = v_s_dust[X[a].sector] / v_s_dust_cout[X[a].sector];
			X[a].e_s = (m_dust*v_s_dust_cout[X[a].sector]) / (m* v_s_cout[X[a].sector]);
			X[a].t_stop_s = (t_stop_s[X[a].sector] / v_s_dust_cout[X[a].sector])/(X[a].ro*K);

		}
		for (int a = N_gran_L_dust; a < N_dust + N_gran_L_dust; a++) {
			X_dust[a].t_stop_s = (t_stop_s[X_dust[a].sector] / v_s_dust_cout[X_dust[a].sector])/(X_dust[a].ro*K);
			X_dust[a].v_s_dust = v_s[X_dust[a].sector] / v_s_cout[X_dust[a].sector];
		}

		*/

		
		/*Расчеты для газа скорости и энергии*/
		
		//#pragma omp parallel for private(Viscosity,gradient,b) 
		for (int a = N_gran_L; a < N + N_gran_L; a++) {	
				//#pragma omp parallel for private(Viscosity,gradient) reduction(+:V_tmp[:10])
				for (int b = a - B; b < a + B; b++) {
					if (((fabs((X[a].x_coordinat[0] - X[b].x_coordinat[0]))) >= 2 * X[a].h_radius_kernel) && (fabs((X[a].x_coordinat[0] - X[b].x_coordinat[0])) >= 2 * X[b].h_radius_kernel))
						continue;
					Viscosity = X[a].Viscosity(X[b], gamma);
					gradient = X[a].W_gradient(X[b]);
					V_tmp[a] = V_tmp[a] - tau*(X[b].m*gradient*(X[a].P / (X[a].ro*X[a].ro) + X[b].P / (X[b].ro*X[b].ro) + Viscosity));
					X[a].U = X[a].U + tau*((X[a].P / (X[a].ro*X[a].ro))  *   X[b].m  *  (X[a].V - X[b].V)  *  gradient + 0.5*X[b].m*Viscosity*(X[a].V - X[b].V)*gradient);
					//	std::cout << K1*(V_tmp[a] - X[a].v_s )<< "\n";
				}
			//	V_tmp[a] = V_tmp[a] - tau*((X[a].e_s/ X[a].t_stop_s)*(V_tmp[a] - X[a].v_s));
			V_tmp[a] = V_tmp[a] - tau*((d_2_g*K*X[a].ro)*(V_tmp[a] - X[a].v_s));

		}
		

		/*расчет скорости пыли*/
	#pragma omp parallel for
	for (int a = N_gran_L_dust; a < N_dust + N_gran_L_dust; a++) {
		//V_tmp_dust[a] =  V_tmp_dust[a] +( (tau*(X_dust[a].ro /X_dust[a].t_stop_s))  )*(X_dust[a].v_s_dust - V_tmp_dust[a]);
		V_tmp_dust[a] = V_tmp_dust[a] + ((tau*(K*X_dust[a].ro)))*(X_dust[a].v_s_dust - V_tmp_dust[a]);

	}
	
	#pragma omp parallel for
	for (int i = N_gran_L; i < N + N_gran_L; i++) {
			X[i].x_coordinat[0] += tau*V_tmp[i];
	}
	
	#pragma omp parallel for
	for (int i = N_gran_L_dust; i < N_dust + N_gran_L_dust; i++) {
			X_dust[i].x_coordinat[0] += tau*V_tmp_dust[i];
		}
	
	
	for (int i = N_gran_L_dust; i < N_dust + N_gran_L_dust - 1; i++) {
		for (int j = 0; j < N_dust + N_gran_L_dust - i - 1; j++) {
			if (X_dust[j].x_coordinat[0] > X_dust[j + 1].x_coordinat[0]) {
				std::swap(X_dust[j], X_dust[j + 1]);
				cout_swap++;
			}

		}
	}


	for (int i = N_gran_L_dust; i < N_dust + N_gran_L_dust; i++) {
		X_dust[i].h_radius_kernel = X_dust[i].m*k / X_dust[i].ro;	
		X_dust[i].ro = X_dust[i].ro_particle(X_dust, B_dust, i);
		X_dust[i].V = V_tmp_dust[i];
	}
	

	for (int i = N_gran_L; i < N + N_gran_L; i++) {
		X[i].h_radius_kernel = X[i].m*k / X[i].ro;
		X[i].ro = X[i].ro_particle(X, B, i);
		X[i].P = X[i].ro*X[i].U*(gamma - 1.0);
		X[i].cs = sqrt((gamma*X[i].P) / X[i].ro);
		X[i].V = V_tmp[i];
	}

	if (j == T / 2) {
		for (int i = 0; i < N + N_gran; i++)
			output_t_0_1 << 1.0 / double(N)*i << " " << X[i].x_coordinat[0]-0.5  << " " << X[i].ro << " " << X[i].V << " " << X[i].P << " " << X[i].U << "\n";
	}
}
	
	
	for (int i = N_gran_L; i < N+N_gran_L; i++){
		X[i].E = X[i].m*(X[i].U + (1.0 / 2.0)*X[i].V*X[i].V);
		output_t_max1 << 1.0 / double(N)*i <<" "<< X[i].x_coordinat[0]-0.5<< " " << X[i].ro << " " << X[i].V<<" "<< X[i].P <<" "<<X[i].U<<" "<<X[i].E<<" " << X[i].h_radius_kernel << "\n";
	}
	for (int i = N_gran_L_dust; i < N_dust + N_gran_L_dust; i++) {
		output_t_max2 << 1.0 / double(N_dust)*i << " " << X_dust[i].x_coordinat[0]-0.5 << " " << X_dust[i].ro << " " << X_dust[i].V <<" "<< X_dust[i].h_radius_kernel<<"\n";
	}
	for (int i = 0; i < N_gran + N; i++){
		X[i].E = X[i].m*(X[i].U + (1.0 / 2.0)*X[i].V*X[i].V);
		t << 1.0 / double(N)*i << " " << X[i].x_coordinat[0]-0.5 << " " << X[i].ro << " " << X[i].V << " " << X[i].P << " " << X[i].U << " " << X[i].E << " " << X[i].h_radius_kernel << "\n";
	}
	std::cout << "Finish";

	output_t_0_1.close();
	output_t_max1.close();
	output_t_max2.close();

	elapsedTime = clock() - elapsedTime;
	printf("It took me %d clicks (%f seconds).\n", elapsedTime, ((float)elapsedTime) / CLOCKS_PER_SEC);

	_getch();
}