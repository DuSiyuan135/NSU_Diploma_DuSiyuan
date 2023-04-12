#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string> 

#define PI 3.141592653589793
#define KC 25
#define LOWERLIMIT 4


using namespace std;

typedef struct Cell {
	double cen_x;
	double cen_y;
	double u;
	double v;
	double con;
	int num_par = KC;
	double sum_par = 0;
}Cell;

typedef struct Particle {
	double x;
	double y;
	double con;
}Particle;

void InitializationCell(int M, int N, int T_, double a, double b, double c, double d, double hx, double hy, double Q, double w, 
	Cell** cell_pp);
void InitializationParticle(int M, int N, int kc, double a, double b, double c, double d, double hx, double hy, 
	Cell** cell_pp, Particle** particle_pp);
void InitializationAddParticles(int* num_empty_p, Particle*** add_particle_ppp, int k);
void ParticleMotion(int M, int N, int kc, int T, int T_, double t, double a, double b, double c, double d, double hx, double hy, double Q, double w,
	Cell** cell_pp, Particle** particle_pp, Particle*** add_particle_ppp, int* num_empty_p);
void CalculationCellConcentration(int M, int N, int kc, int T_, double a, double b, double c, double d, double hx, double hy, double tau,
	Cell** cell_pp, Particle** particle_pp, Particle*** add_particle_ppp, int* num_empty_p, int k);
void AddParticles(int index_i, int index_j, int index_k, int kc, int T_, double a, double c, double hx, double hy, double t, 
	Particle*** add_particle_ppp, Cell** cell_pp, int k);
void PrintCellCenterCoordinatesX(int M, int N, Cell** cell_pp);
void PrintCellCenterCoordinatesY(int M, int N, Cell** cell_pp);
void PrintCellCenterCoordinates(int M, int N, Cell** cell_pp);
void PrintCellConcentration(int M, int N, Cell** cell_pp, int k);
void PrintCellVelocity(int M, int N, Cell** cell_pp, int k);
void PrintCellVelocity_u(int M, int N, Cell** cell_pp, int k);
void PrintCellVelocity_v(int M, int N, Cell** cell_pp, int k);
void PrintParticleCenterCoordinatesX(int M, int N, int kc, Particle** particle_pp, int k);
void PrintParticleCenterCoordinatesY(int M, int N, int kc, Particle** particle_pp, int k);
void PrintParticleConcentration(int M, int N, int kc, Particle** particle_pp, int k);
double ConcentrationFunction(double x, double y, double t, int T_);
double VelocityFunction_u(double x, double y, double Q, double w);
double VelocityFunction_v(double x, double y, double Q, double w);

int main()
{
	double t, a, b, c, d, hx, hy, Q, w;
	int N, M, T, T_, kc = 5;

	Q = 0.1; w = 0.1;

	a = -100; b = 100;
	c = -100; d = 100;
	N = 200; M = 200;
	T = 3000; T_ = 300;
	t = 1500;
	
	hx = (b - a) / N;
	hy = (d - c) / M;

	Cell** cell_pp = new Cell*[M];
	for (int j = 0; j < M; j++)
	{
		cell_pp[j] = new Cell[N];
	}

	Particle** particle_pp = new Particle*[M*kc];
	for (int j = 0; j < M*kc; j++)
	{
		particle_pp[j] = new Particle[N*kc];
	}

	Particle*** add_particle_ppp = new Particle**[T + 1];

	int* num_empty_p = new int[T + 1];
	for (int k = 0; k <= T; k++)
	{
		num_empty_p[k] = 0;
	}

	InitializationCell(M, N, T_, a, b, c, d, hx, hy, Q, w, cell_pp);
	InitializationParticle(M, N, kc, a, b, c, d, hx, hy, cell_pp, particle_pp);

	// PrintCellCenterCoordinatesX(M, N, cell_pp);
	// PrintCellCenterCoordinatesY(M, N, cell_pp);
	// PrintCellCenterCoordinates(M, N, cell_pp);

	// PrintCellConcentration(M, N, cell_pp, 0);

	// PrintCellVelocity_u(M, N, cell_pp, 0);
	// PrintCellVelocity_v(M, N, cell_pp, 0);

	// PrintParticleCenterCoordinatesX(M, N, kc, particle_pp, 0);
	// PrintParticleCenterCoordinatesY(M, N, kc, particle_pp, 0);
	// PrintParticleConcentration(M, N, kc, particle_pp, 0);

	ParticleMotion(M, N, kc, T, T_, t, a, b, c, d, hx, hy, Q, w, cell_pp, particle_pp, add_particle_ppp, num_empty_p);

	// PrintCellConcentration(M, N, cell_pp, T);

	// PrintCellVelocity_u(M, N, cell_pp, T);
	// PrintCellVelocity_v(M, N, cell_pp, T);

	// PrintParticleCenterCoordinatesX(M, N, kc, particle_pp, T);
	// PrintParticleCenterCoordinatesY(M, N, kc, particle_pp, T);
	// PrintParticleConcentration(M, N, kc, particle_pp, T);



	for (int j = 0; j < M; j++)
	{
		delete[] cell_pp[j];
	}
	delete[] cell_pp;

	for (int j = 0; j < M*kc; j++)
	{
		delete[] particle_pp[j];
	}
	delete[] particle_pp;

	for (int k = 1; k <= T; k++)
	{
		for (int i = 0; i < num_empty_p[k]; i++)
		{
			delete[] add_particle_ppp[k][i];
		}
	}
	for (int k = 1; k <= T; k++)
	{
		delete[] add_particle_ppp[k];
	}
	delete[] add_particle_ppp;

	delete[] num_empty_p;


	system("pause");

	return 0;
}

void InitializationCell(int M, int N, int T_, double a, double b, double c, double d, double hx, double hy, double Q, double w, Cell** cell_pp)
{
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			cell_pp[j][i].cen_x = a + i*hx + 0.5*hx;
			cell_pp[j][i].cen_y = c + j*hy + 0.5*hy;

			cell_pp[j][i].con = ConcentrationFunction(cell_pp[j][i].cen_x, cell_pp[j][i].cen_y, 0, T_);

			cell_pp[j][i].u = VelocityFunction_u(cell_pp[j][i].cen_x, cell_pp[j][i].cen_y, Q, w);
			cell_pp[j][i].v = VelocityFunction_v(cell_pp[j][i].cen_x, cell_pp[j][i].cen_y, Q, w);
		}
	}
}

void InitializationParticle(int M, int N, int kc, double a, double b, double c, double d, double hx, double hy, Cell** cell_pp, Particle** particle_pp)
{
	for (int j = 0; j < M*kc; j++)
	{
		for (int i = 0; i < N*kc; i++)
		{
			particle_pp[j][i].x = a + i*hx / kc + 0.5*hx / kc;
			particle_pp[j][i].y = c + j*hy / kc + 0.5*hy / kc;
			// cout << (int)floor((particle_pp[j][i].y - c) / hy) << "," << (int)floor((particle_pp[j][i].x - a) / hx) << " ";
			particle_pp[j][i].con = cell_pp[(int)floor((particle_pp[j][i].y - c) / hy)][(int)floor((particle_pp[j][i].x - a) / hx)].con;
		}
		// cout << endl;
	}
}

void ParticleMotion(int M, int N, int kc, int T, int T_, double t, double a, double b, double c, double d, double hx, double hy, double Q, double w,
	Cell** cell_pp, Particle** particle_pp, Particle*** add_particle_ppp, int* num_empty_p)
{
	double u = 0, v = 0;
	double index_x, index_y;
	double tau = t / T;

	PrintCellConcentration(M, N, cell_pp, 0);

	for (int k = 1; k <= T; k++)
	{
		// Compute the motion of the initial particle.
		for (int j = 0; j < M*kc; j++)
		{
			for (int i = 0; i < N*kc; i++)
			{
				// cout << (int)floor((particle_pp[j][i].y - c) / hy) << "," << (int)floor((particle_pp[j][i].x - a) / hx) << " ";
				index_x = (int)floor((particle_pp[j][i].x - a) / hx);
				index_y = (int)floor((particle_pp[j][i].y - c) / hy);
				if (index_x < N && index_x >= 0 && index_y < M && index_y >= 0)
				{
					// cell_pp[(int)index_y][(int)index_x].u = VelocityFunction_u(cell_pp[(int)index_y][(int)index_x].cen_x, cell_pp[(int)index_y][(int)index_x].cen_y, k*tau, omega);
					// cell_pp[(int)index_y][(int)index_x].v = VelocityFunction_v(cell_pp[(int)index_y][(int)index_x].cen_x, cell_pp[(int)index_y][(int)index_x].cen_y, k*tau, omega);
					// u = cell_pp[(int)index_y][(int)index_x].u;
					// v = cell_pp[(int)index_y][(int)index_x].v;
					u = VelocityFunction_u(cell_pp[(int)index_y][(int)index_x].cen_x, cell_pp[(int)index_y][(int)index_x].cen_y, Q, w);
					v = VelocityFunction_v(cell_pp[(int)index_y][(int)index_x].cen_x, cell_pp[(int)index_y][(int)index_x].cen_y, Q, w);
					particle_pp[j][i].x += u*tau;
					particle_pp[j][i].y += v*tau;
				}
			}
			//cout << endl;
		}

		// Calculate the motion of the newly added particles.
		for (int kk = 0; kk < k; kk++)
		{
			for (int i = 0; i < num_empty_p[kk]; i++)
			{
				for (int j = 0; j < KC; j++)
				{
					index_x = (int)floor((add_particle_ppp[kk][i][j].x - a) / hx);
					index_y = (int)floor((add_particle_ppp[kk][i][j].y - c) / hy);
					if (index_x < N && index_x >= 0 && index_y < M && index_y >= 0)
					{
						// cell_pp[(int)index_y][(int)index_x].u = VelocityFunction_u(cell_pp[(int)index_y][(int)index_x].cen_x, cell_pp[(int)index_y][(int)index_x].cen_y, k*tau, omega);
						// cell_pp[(int)index_y][(int)index_x].v = VelocityFunction_v(cell_pp[(int)index_y][(int)index_x].cen_x, cell_pp[(int)index_y][(int)index_x].cen_y, k*tau, omega);
						// u = cell_pp[(int)index_y][(int)index_x].u;
						// v = cell_pp[(int)index_y][(int)index_x].v;
						u = VelocityFunction_u(cell_pp[(int)index_y][(int)index_x].cen_x, cell_pp[(int)index_y][(int)index_x].cen_y, Q, w);
						v = VelocityFunction_v(cell_pp[(int)index_y][(int)index_x].cen_x, cell_pp[(int)index_y][(int)index_x].cen_y, Q, w);
						add_particle_ppp[kk][i][j].x += u*tau;
						add_particle_ppp[kk][i][j].y += v*tau;
					}
				}
			}
		}

		CalculationCellConcentration(M, N, kc, T_, a, b, c, d, hx, hy, tau, cell_pp, particle_pp, add_particle_ppp, num_empty_p, k);
	}
}

void CalculationCellConcentration(int M, int N, int kc, int T_, double a, double b, double c, double d, double hx, double hy, double tau,
	Cell** cell_pp, Particle** particle_pp, Particle*** add_particle_ppp, int* num_empty_p, int k)
{
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			cell_pp[j][i].num_par = 0;
			cell_pp[j][i].sum_par = 0;
		}
	}

	// Calculate the initial particle number and the sum of initial particle concentrations in each cell.
	double index_x, index_y;
	for (int j = 0; j < M*kc; j++)
	{
		for (int i = 0; i < N*kc; i++)
		{
			index_x = (int)floor((particle_pp[j][i].x - a) / hx);
			index_y = (int)floor((particle_pp[j][i].y - c) / hy);
			if (index_x < N && index_x >= 0 && index_y < M && index_y >= 0)
			{
				cell_pp[(int)index_y][(int)index_x].num_par++;
				cell_pp[(int)index_y][(int)index_x].sum_par += particle_pp[j][i].con;
			}
		}
	}

	// Calculate the number of added particles and the sum of added particle concentrations in each cell.
	for (int kk = 0; kk < k; kk++)
	{
		for (int i = 0; i < num_empty_p[kk]; i++)
		{
			for (int j = 0; j < (kc*kc); j++)
			{
				index_x = (int)floor((add_particle_ppp[kk][i][j].x - a) / hx);
				index_y = (int)floor((add_particle_ppp[kk][i][j].y - c) / hy);
				if (index_x < N && index_x >= 0 && index_y < M && index_y >= 0)
				{
					cell_pp[(int)index_y][(int)index_x].num_par++;
					cell_pp[(int)index_y][(int)index_x].sum_par += add_particle_ppp[kk][i][j].con;
				}
			}
		}
	}

	// If the number of particles in the cell is not zero, calculate the concentration. 
	// Otherwise count the number of cells with zero particles.
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			// if (cell_pp[j][i].num_par == 0)
			if (cell_pp[j][i].num_par < LOWERLIMIT)
			{
				num_empty_p[k]++;
				// cell_pp[j][i].con = -1;
				if (cell_pp[j][i].num_par != 0)
				{
					// cell_pp[j][i].con = ConcentrationFunction(cell_pp[j][i].cen_x, cell_pp[j][i].cen_y, k*tau, T_);
					cell_pp[j][i].con = cell_pp[j][i].sum_par / cell_pp[j][i].num_par;
				}
			}
			else // (cell_pp[j][i].num_par != 0)
			{
				cell_pp[j][i].con = cell_pp[j][i].sum_par / cell_pp[j][i].num_par;
			}
		}
	}

	InitializationAddParticles(num_empty_p, add_particle_ppp, k);

	// If there are cells with zero particles, add new particles to them.
	if (num_empty_p[k] > 0)
	{
		int index_k = 0;
		int flag = 0;
		// cout << "k= " << k << "\t num_empty= " << num_empty_p[k] << endl;
		for (int j = 0; j < M; j++)
		{
			for (int i = 0; i < N; i++)
			{
				//if (cell_pp[j][i].num_par == 0)
				if (cell_pp[j][i].num_par < LOWERLIMIT)
				{
					AddParticles(i, j, index_k, kc, T_, a, c, hx, hy, k*tau, add_particle_ppp, cell_pp, k);
					index_k++;
					if (index_k == num_empty_p[k])
					{
						flag = 1;
						break;
					}
				}
			}
			if (flag == 1)
			{
				break;
			}
		}
	}

	PrintCellConcentration(M, N, cell_pp, k);
}

void InitializationAddParticles(int* num_empty_p, Particle*** add_particle_ppp, int k)
{
	add_particle_ppp[k] = new Particle*[num_empty_p[k]];
	for (int kk = 0; kk < num_empty_p[k]; kk++)
	{
		add_particle_ppp[k][kk] = new Particle[KC];
	}
	cout << "In the step" << k << ", " << num_empty_p[k] << "*" << KC << " particles were added." << endl;
}

void AddParticles(int index_i, int index_j, int index_k, int kc, int T_, double a, double c, double hx, double hy, double t, Particle*** add_particle_ppp, Cell** cell_pp, int k)
{
	for (int i = 0; i < KC; i++)
	{
		add_particle_ppp[k][index_k][i].x = a + index_i*hx + hx / kc / 2 + i%kc*hx / kc;
		add_particle_ppp[k][index_k][i].y = c + index_j*hy + hy / kc / 2 + i/kc*hy / kc;
		//add_particle_ppp[k][index_k][i].con = ConcentrationFunction(add_particle_ppp[k][index_k][i].x, add_particle_ppp[k][index_k][i].y, t, T_);

		if (add_particle_ppp[k][index_k][i].x*add_particle_ppp[k][index_k][i].x + add_particle_ppp[k][index_k][i].y*add_particle_ppp[k][index_k][i].y <= 9)
		{
			add_particle_ppp[k][index_k][i].con = ConcentrationFunction(add_particle_ppp[k][index_k][i].x, add_particle_ppp[k][index_k][i].y, t, T_);
		}
		else
		{
			add_particle_ppp[k][index_k][i].con = cell_pp[index_i][index_j].con;
		}
	}
}

void PrintCellCenterCoordinatesX(int M, int N, Cell** cell_pp)
{
	ofstream OutFile("Data\\Cell Center Coordinates X.txt");
	// OutFile << "Cell Center Coordinates X" << endl;
	cout << "Cell Center Coordinates X" << endl;
	for (int j = M - 1; j >= 0; j--)
	{
		for (int i = 0; i < N; i++)
		{
			cout << cell_pp[j][i].cen_x << " ";
			OutFile << cell_pp[j][i].cen_x << " ";
		}
		cout << endl;
		OutFile << endl;
	}

	OutFile.close();
}

void PrintCellCenterCoordinatesY(int M, int N, Cell** cell_pp)
{
	ofstream OutFile("Data\\Cell Center Coordinates Y.txt");
	// OutFile << "Cell Center Coordinates Y" << endl;
	cout << "Cell Center Coordinates Y" << endl;
	for (int j = M - 1; j >= 0; j--)
	{
		for (int i = 0; i < N; i++)
		{
			cout << cell_pp[j][i].cen_y << " ";
			OutFile << cell_pp[j][i].cen_y << " ";
		}
		cout << endl;
		OutFile << endl;
	}

	OutFile.close();
}

void PrintCellCenterCoordinates(int M, int N, Cell** cell_pp)
{
	ofstream OutFile("Data\\Cell Center Coordinates.txt");
	// OutFile << "Cell Center Coordinates" << endl;
	cout << "Cell Center Coordinates" << endl;
	for (int j = M - 1; j >= 0; j--)
	{
		for (int i = 0; i < N; i++)
		{
			cout << cell_pp[j][i].cen_x << "," << cell_pp[j][i].cen_y << " ";
			OutFile << cell_pp[j][i].cen_x << "," << cell_pp[j][i].cen_y << " ";
		}
		cout << endl;
		OutFile << endl;
	}

	OutFile.close();
}

void PrintCellConcentration(int M, int N, Cell** cell_pp, int k)
{
	string name = "k=" + to_string(k) + ".txt";
	ofstream OutFile("Data\\Cell Concentration\\" + name);
	// OutFile << "Concentration" << endl;
	// cout << "Concentration" << endl;
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			// cout << cell_pp[j][i].con << " ";
			OutFile << cell_pp[j][i].con << " ";
		}
		// cout << endl;
		OutFile << endl;
	}
	OutFile.close();
}

void PrintCellVelocity(int M, int N, Cell** cell_pp, int k)
{
	string name = "k=" + to_string(k) + ".txt";
	ofstream OutFile("Data\\Cell Velocity\\" + name);

	OutFile << "Velocity_u" << endl;
	cout << "Velocity_u" << endl;
	for (int j = M - 1; j >= 0; j--)
	{
		for (int i = 0; i < N; i++)
		{
			cout << cell_pp[j][i].u << " ";
			OutFile << cell_pp[j][i].u << " ";
		}
		cout << endl;
		OutFile << endl;
	}

	OutFile << "Velocity_v" << endl;
	cout << "Velocity_v" << endl;
	for (int j = M - 1; j >= 0; j--)
	{
		for (int i = 0; i < N; i++)
		{
			cout << cell_pp[j][i].v << " ";
			OutFile << cell_pp[j][i].v << " ";
		}
		cout << endl;
		OutFile << endl;
	}
	OutFile.close();
}

void PrintCellVelocity_u(int M, int N, Cell** cell_pp, int k)
{
	string name = "k=" + to_string(k) + " u.txt";
	ofstream OutFile("Data\\Cell Velocity\\" + name);

	// OutFile << "Velocity_u" << endl;
	cout << "Velocity_u" << endl;
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			cout << cell_pp[j][i].u << " ";
			OutFile << cell_pp[j][i].u << " ";
		}
		cout << endl;
		OutFile << endl;
	}
	OutFile.close();
}

void PrintCellVelocity_v(int M, int N, Cell** cell_pp, int k)
{
	string name = "k=" + to_string(k) + " v.txt";
	ofstream OutFile("Data\\Cell Velocity\\" + name);

	// OutFile << "Velocity_v" << endl;
	cout << "Velocity_v" << endl;
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			cout << cell_pp[j][i].v << " ";
			OutFile << cell_pp[j][i].v << " ";
		}
		cout << endl;
		OutFile << endl;
	}
	OutFile.close();
}


void PrintParticleCenterCoordinatesX(int M, int N, int kc, Particle** particle_pp, int k)
{
	string name = "k=" + to_string(k) + " x.txt";
	ofstream OutFile("Data\\Particle Coordinates\\" + name);
	//cout << "Particle Coordinates X" << endl;
	for (int j = 0; j < M*kc; j++)
	{
		for (int i = 0; i < N*kc; i++)
		{
			//cout << particle_pp[j][i].x << " ";
			OutFile << particle_pp[j][i].x << " ";
		}
		//cout << endl;
		OutFile << endl;
	}
	OutFile.close();
}

void PrintParticleCenterCoordinatesY(int M, int N, int kc, Particle** particle_pp, int k)
{
	string name = "k=" + to_string(k) + " y.txt";
	ofstream OutFile("Data\\Particle Coordinates\\" + name);
	//cout << "Particle Coordinates Y" << endl;
	for (int j = 0; j < M*kc; j++)
	{
		for (int i = 0; i < N*kc; i++)
		{
			//cout << particle_pp[j][i].y << " ";
			OutFile << particle_pp[j][i].y << " ";
		}
		//cout << endl;
		OutFile << endl;
	}
	OutFile.close();
}

void PrintParticleConcentration(int M, int N, int kc, Particle** particle_pp, int k)
{
	string name = "k=" + to_string(k) + ".txt";
	ofstream OutFile("Data\\Particle Concentration\\" + name);
	cout << "Particle Concentration" << endl;
	for (int j = 0; j < M*kc; j++)
	{
		for (int i = 0; i < N*kc; i++)
		{
			//cout << particle_pp[j][i].con << " ";
			OutFile << particle_pp[j][i].con << " ";
		}
		//cout << endl;
		OutFile << endl;
	}
	OutFile.close();
}

double ConcentrationFunction(double x, double y, double t, int T_)
{
	if (x*x + y*y <= 9)
	{
		// fmod: division with remainder of type double
		// 2kT_ <= t <= (2k+1)T_
		if ((int)floor(t / T_)%2==0 || fmod(t, T_)==0) 
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		return 0;
	}
}

double VelocityFunction_u(double x, double y, double Q, double w)
{
	double u = Q * x / (2 * PI*w*sqrt(x*x + y*y));
	return u;
}

double VelocityFunction_v(double x, double y, double Q, double w)
{
	double v = Q * y / (2 * PI*w*sqrt(x*x + y*y));
	return v;
}