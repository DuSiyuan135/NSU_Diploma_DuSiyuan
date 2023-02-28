#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string> 


using namespace std;

void CellCenterCoordinates(double a, double b, double c, double d, double hx, double hy);
void Concentration(int M, int N, double a, double b, double c, double d, double hx, double hy);
void Velocity(int M, int N, double a, double b, double c, double d, double hx, double hy);
double ConcentrationFunction(double x, double y);
double VelocityFunction_u(double x, double y);
double VelocityFunction_v(double x, double y);

int main()
{
	double a, b, c, d, hx, hy;
	int N, M, T;

	a = -100; b = 100;
	c = -100; d = 100;
	N = 100; M = 100;
	T = 100;

	hx = (b - a) / N;
	hy = (d - c) / M;

	
	CellCenterCoordinates(a, b, c, d, hx, hy);
	Concentration(M, N, a, b, c, d, hx, hy);
	Velocity(M, N, a, b, c, d, hx, hy);
	
	//ofstream outfile("data\\cell boundary coordinates.txt");
	//outfile << "cell boundary coordinates" << endl;
	//cout << "cell boundary coordinates" << endl;
	//for (double y = d; y  >= c; y -= hy)
	//{
	//	for (double x = a; x <= d; x += hx)
	//	{
	//		cout << x << "," << y << " ";
	//		outfile << x << "," << y << " ";
	//	}
	//	cout << endl;
	//	outfile << endl;
	//}
	//outfile.close();

	system("pause");

	return 0;
}

//void File1(int M, int N, double a, double b, double c, double d, double hx, double hy)
//{
//	struct CenterCoordinates
//	{
//		double x;
//		double y;
//	};
//
//	struct Concentration
//	{
//		double x;
//		double y;
//	};
//
//	struct Velocity
//	{
//		double x;
//		double y;
//	};
//
//	struct Cell
//	{
//		struct CenterCoordinates cc;
//		struct Concentration con;
//		struct Velocity c_v;
//	};
//
//}

void CellCenterCoordinates(double a, double b, double c, double d, double hx, double hy)
{
	ofstream OutFile("Data\\Cell Center Coordinates.txt");
	OutFile << "Cell Center Coordinates" << endl;
	cout << "Cell Center Coordinates" << endl;
	for (double y = d-0.5*hy; y >= c; y -= hy)
	{
		for (double x = a+0.5*hx; x <= d; x += hx)
		{
			cout << x << "," << y << " ";
			OutFile << x << "," << y << " ";
		}
		cout << endl;
		OutFile << endl;
	}
	OutFile.close();
}

void Concentration(int M, int N, double a, double b, double c, double d, double hx, double hy)
{
	double** con = new double*[N];
	for (int j = 0; j < N; j++)
	{
		con[j] = new double[M];
	}

	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			con[j][i] = ConcentrationFunction(a + 0.5*hx + i*hx, c + hy + 0.5*hy + j*hy);
		}
	}
	ofstream OutFile("Data\\Cell Center Coordinates.txt", ios::app);
	OutFile << "Concentration" << endl;
	cout << "Concentration" << endl;
	for (int i = M-1; i >= 0; i--)
	{
		for (int j = 0; j < N; j++)
		{
			cout << con[i][j] << " ";
			OutFile << con[i][j] << " ";
		}
		cout << endl;
		OutFile << endl;
	}
	OutFile.close();

	for (int i = 0; i < M; i++)
	{
		delete[] con[i];
	}
	delete[] con;
}

void Velocity(int M, int N, double a, double b, double c, double d, double hx, double hy)
{
	double** u= new double*[N];
	double** v = new double*[N];
	for (int j = 0; j < N; j++)
	{
		u[j] = new double[M];
		v[j] = new double[M];
	}

	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			u[j][i] = VelocityFunction_u(a + 0.5*hx + i*hx, c + hy + 0.5*hy + j*hy);
			v[j][i] = VelocityFunction_v(a + 0.5*hx + i*hx, c + hy + 0.5*hy + j*hy);
		}
	}
	ofstream OutFile("Data\\Cell Center Coordinates.txt", ios::app);
	OutFile << "Velocity_u" << endl;
	cout << "Velocity_u" << endl;
	for (int i = M - 1; i >= 0; i--)
	{
		for (int j = 0; j < N; j++)
		{
			cout << u[i][j] << " ";
			OutFile << u[i][j] << " ";
		}
		cout << endl;
		OutFile << endl;
	}

	OutFile << "Velocity_v" << endl;
	cout << "Velocity_v" << endl;
	for (int i = M - 1; i >= 0; i--)
	{
		for (int j = 0; j < N; j++)
		{
			cout << v[i][j] << " ";
			OutFile << v[i][j] << " ";
		}
		cout << endl;
		OutFile << endl;
	}
	OutFile.close();

	for (int i = 0; i < M; i++)
	{
		delete[] u[i];
		delete[] v[i];
	}
	delete[] u;
	delete[] v;
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