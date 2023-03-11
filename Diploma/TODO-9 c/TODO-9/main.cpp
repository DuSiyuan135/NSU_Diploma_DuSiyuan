#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string> 


using namespace std;

typedef struct Cell {
	double cen_x;
	double cen_y;
	double u;
	double v;
	double con;
	int num_par = 0;
	double sum_par = 0;
}Cell;

typedef struct Particle {
	double x;
	double y;
	double con;
}Particle;

void InitializationCell(int M, int N, double a, double b, double c, double d, double hx, double hy, Cell** cell_pp);
void InitializationParticle(int M, int N, int kc, double a, double b, double c, double d, double hx, double hy, Cell** cell_pp, Particle** particle_pp);
void ParticleMotion(int M, int N, int kc, int T, double t, double a, double b, double c, double d, double hx, double hy, Cell** cell_pp, Particle** particle_pp);
void CalculationCellConcentration(int M, int N, int kc, double a, double b, double c, double d, double hx, double hy, Cell** cell_pp, Particle** particle_pp);
void PrintCellCenterCoordinatesX(int M, int N, Cell** cell_pp);
void PrintCellCenterCoordinatesY(int M, int N, Cell** cell_pp);
void PrintCellCenterCoordinates(int M, int N, Cell** cell_pp);
void PrintCellConcentration(int M, int N, Cell** cell_pp);
void PrintCellVelocity(int M, int N, Cell** cell_pp);
double ConcentrationFunction(double x, double y);
double VelocityFunction_u(double x, double y);
double VelocityFunction_v(double x, double y);

int main()
{
	double t, a, b, c, d, hx, hy;
	int N, M, T, kc = 5;

	a = -100; b = 100;
	c = -100; d = 100;
	N = 100; M = 100;
	T = 100; t = 100;

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

	InitializationCell(M, N, a, b, c, d, hx, hy, cell_pp);
	InitializationParticle(M, N, kc, a, b, c, d, hx, hy, cell_pp, particle_pp);

	ParticleMotion(M, N, kc, T, t, a, b, c, d, hx, hy, cell_pp, particle_pp);
	CalculationCellConcentration(M, N, kc, a, b, c, d, hx, hy, cell_pp, particle_pp);

	//PrintCellCenterCoordinatesX(M, N, cell_pp);
	//PrintCellCenterCoordinatesY(M, N, cell_pp);
	//PrintCellCenterCoordinates(M, N, cell_pp);
	PrintCellConcentration(M, N, cell_pp);
	//PrintCellVelocity(M, N, cell_pp);

	for (int j = 0; j < M; j++)
	{
		delete[] cell_pp[j];
	}
	delete[] cell_pp;


	system("pause");

	return 0;
}

void InitializationCell(int M, int N, double a, double b, double c, double d, double hx, double hy, Cell** cell_pp)
{
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			cell_pp[j][i].cen_x = a + i*hx + 0.5*hx;
			cell_pp[j][i].cen_y = c + j*hy + 0.5*hy;

			cell_pp[j][i].con = ConcentrationFunction(cell_pp[j][i].cen_x, cell_pp[j][i].cen_y);

			cell_pp[j][i].u = VelocityFunction_u(cell_pp[j][i].cen_x, cell_pp[j][i].cen_y);
			cell_pp[j][i].v = VelocityFunction_v(cell_pp[j][i].cen_x, cell_pp[j][i].cen_y);
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
			//cout << (int)floor((particle_pp[j][i].y - c) / hy) << "," << (int)floor((particle_pp[j][i].x - a) / hx) << " ";
			particle_pp[j][i].con = cell_pp[(int)floor((particle_pp[j][i].y - c) / hy)][(int)floor((particle_pp[j][i].x - a) / hx)].con;
		}
		//cout << endl;
	}
}

void ParticleMotion(int M, int N, int kc, int T, double t, double a, double b, double c, double d, double hx, double hy, Cell** cell_pp, Particle** particle_pp)
{
	double u = 0, v = 0;
	double index_x, index_y;
	double tau = t / T;
	for (int k = 0; k <= T; k++)
	{
		for (int j = 0; j < M*kc; j++)
		{
			for (int i = 0; i < N*kc; i++)
			{
				//cout << (int)floor((particle_pp[j][i].y - c) / hy) << "," << (int)floor((particle_pp[j][i].x - a) / hx) << " ";
				index_x = (int)floor((particle_pp[j][i].x - a) / hx);
				index_y = (int)floor((particle_pp[j][i].y - c) / hy);
				if (index_x < N && index_y < M)
				{
					u = cell_pp[(int)index_y][(int)index_x].u;
					v = cell_pp[(int)index_y][(int)index_x].v;
					particle_pp[j][i].x += u*tau;
					particle_pp[j][i].y += v*tau;
				}
			}
			//cout << endl;
		}
	}
}

void CalculationCellConcentration(int M, int N, int kc, double a, double b, double c, double d, double hx, double hy, Cell** cell_pp, Particle** particle_pp)
{
	double index_x, index_y;
	for (int j = 0; j < M*kc; j++)
	{
		for (int i = 0; i < N*kc; i++)
		{
			index_x = (int)floor((particle_pp[j][i].x - a) / hx);
			index_y = (int)floor((particle_pp[j][i].y - c) / hy);
			if (index_x < N && index_y < M)
			{
				cell_pp[(int)index_y][(int)index_x].num_par++;
				cell_pp[(int)index_y][(int)index_x].sum_par += particle_pp[j][i].con;
			}
		}
	}

	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			if (cell_pp[j][i].num_par == 0)
			{
				cell_pp[j][i].con = 0;
			}
			else
			{
				cell_pp[j][i].con = cell_pp[j][i].sum_par / cell_pp[j][i].num_par;
			}
		}
	}
}

void PrintCellCenterCoordinatesX(int M, int N, Cell** cell_pp)
{
	ofstream OutFile("Data\\Cell Center Coordinates X.txt");
	OutFile << "Cell Center Coordinates X" << endl;
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
	OutFile << "Cell Center Coordinates Y" << endl;
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
	OutFile << "Cell Center Coordinates" << endl;
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

void PrintCellConcentration(int M, int N, Cell** cell_pp)
{
	ofstream OutFile("Data\\Cell Concentration.txt");
	OutFile << "Concentration" << endl;
	cout << "Concentration" << endl;
	for (int j = M - 1; j >= 0; j--)
	{
		for (int i = 0; i < N; i++)
		{
			cout << cell_pp[j][i].con << " ";
			OutFile << cell_pp[j][i].con << " ";
		}
		cout << endl;
		OutFile << endl;
	}
	OutFile.close();
}

void PrintCellVelocity(int M, int N, Cell** cell_pp)
{
	ofstream OutFile("Data\\Cell Velocity.txt");
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

double ConcentrationFunction(double x, double y)
{
	if (x <= 10 && x >= -10)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

double VelocityFunction_u(double x, double y)
{
	return 0.1;
}

double VelocityFunction_v(double x, double y)
{
	return 0.1;
}