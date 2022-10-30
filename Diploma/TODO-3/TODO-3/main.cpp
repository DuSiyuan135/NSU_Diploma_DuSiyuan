#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string> 

using namespace std;

double Wide_w(double x, double w0, double cL, double L);
double Velocity_u(double x, double w, double Q);
double AccurateLocation(double t, double w0, double cL, double L, double Q);
double Euler(double tau, double x, double w0, double cL, double L, double Q);
double RungeKuttaFourthOrder(double tau, double x, double w0, double cL, double L, double Q);
double IterativeMethod(double A, double B, double x0);


int main()
{
	double u0 = 0.1, a = 0, b = 100, Q = 0.1, w0 = 0.1, cL = 0.2, L = 100;

	cout << "Euler:" << endl;
	for (double tau = 2; tau >= 0.5; tau /= 2)
	{
		cout << "tau= " << tau << "\t";
		double t = 0, err = 0, x_E = 0;
		x_E = AccurateLocation(t, w0, cL, L, Q);

		string site = "tau=" + to_string(tau);
		ofstream OutFile("Result\\Euler\\" + site + ".txt");
		OutFile << "tau=" << tau << endl;
		OutFile << "t=" << t << " \t" << "x=" << a << " \t" << "x_E=" << x_E << " \t" << "Err=" << err << endl;
		// printf("t=%.5lf\t x_E=%.5lf\t x=%.5lf\t err=%.5lf\t\n", t, a, AccurateLocation(t, w0, cL, L, Q), 0);
		for (double x = a; x_E < b;)
		{
			x = Euler(tau, x, w0, cL, L, Q);
			t += tau;
			x_E = AccurateLocation(t, w0, cL, L, Q);
			if (err < abs(x_E - x))
				err = abs(x_E - x);
			// printf("t=%.5lf\t x_E=%.5lf\t x=%.5lf\t err=%.5lf\t", t, x, AccurateLocation(t, w0, cL, L, Q), abs(AccurateLocation(t, w0, cL, L, Q) - x));
			// cout << endl;
			OutFile << "t=" << t << " \t" << "x=" << x << " \t" << "x_E=" << x_E << " \t" << "Err=" << err << endl;
		}
		printf("err=%.7lf\t t=%.5lf\t", err, t);
		cout << endl;

		OutFile.close();
	}

	cout << "RungeKutta:" << endl;
	for (double tau = 2; tau >= 0.5; tau /= 2)
	{
		cout << "tau= " << tau << "\t";
		double t = 0, err = 0, x_E = 0;
		x_E = AccurateLocation(t, w0, cL, L, Q);

		string site = "tau=" + to_string(tau);
		ofstream OutFile("Result\\RungeKutta\\" + site + ".txt");
		OutFile << "tau=" << tau << endl;
		OutFile << "t=" << t << " \t" << "x=" << a << " \t" << "x_E=" << x_E << " \t" << "Err=" << err << endl;
		// printf("t=%.5lf\t x_E=%.5lf\t x=%.5lf\t err=%.5lf\t\n", t, a, AccurateLocation(t, w0, cL, L, Q), 0);
		for (double x = a; x_E < b;)
		{
			x = RungeKuttaFourthOrder(tau, x, w0, cL, L, Q);
			t += tau;
			x_E = AccurateLocation(t, w0, cL, L, Q);
			if (err < abs(x_E - x))
				err = abs(x_E - x);
			// printf("t=%.5lf\t x_E=%.5lf\t x=%.5lf\t err=%.5lf\t", t, x, AccurateLocation(t, w0, cL, L, Q), abs(AccurateLocation(t, w0, cL, L, Q) - x));
			// cout << endl;
			OutFile << "t=" << t << " \t" << "x=" << x << " \t" << "x_E=" << x_E << " \t" << "Err=" << err << endl;
		}
		printf("err=%.7lf\t t=%.5lf\t", err, t);
		cout << endl;

		OutFile.close();
	}

	system("pause");

	return 0;
}

double Wide_w(double x, double w0, double cL, double L)
{
	return w0*(1 - ((1 - cL)*x / L)*((1 - cL)*x / L));
}

double Velocity_u(double x, double w, double Q)
{
	return Q / w;
}

double AccurateLocation(double t, double w0, double cL, double L, double Q)
{
	double A, B, x0;
	A = 1/(3*L*L)*(1-cL)*(1-cL);
	B = Q*t/w0;
	x0 = 50;

	return IterativeMethod(A, B, x0);
}

double IterativeMethod(double A, double B, double x0) // x = Ax^3+B
{
	double x_k, x_k1; 
	x_k = x0;

	while (1)
	{
		x_k1 = A*x_k*x_k*x_k + B;
		if (abs(x_k1 - x_k) < 1e-7)
		{
			return x_k1;
		}
		x_k = x_k1;
	}
}

double Euler(double tau, double x, double w0, double cL, double L, double Q)
{
	double u = Velocity_u(x, Wide_w(x, w0, cL, L), Q);

	return (x + tau*u);
}

double RungeKuttaFourthOrder(double tau, double x, double w0, double cL, double L, double Q)
{
	double k1, k2, k3, k4;

	k1 = Velocity_u(x, Wide_w(x, w0, cL, L), Q);
	k2 = Velocity_u(x + tau / 2 * k1, Wide_w(x + tau / 2 * k1, w0, cL, L), Q);
	k3 = Velocity_u(x + tau / 2 * k2, Wide_w(x + tau / 2 * k2, w0, cL, L), Q);
	k4 = Velocity_u(x + tau * k3, Wide_w(x + tau * k3, w0, cL, L), Q);

	return (x + tau / 6 * (k1 + 2 * k2 + 2 * k3 + k4));
}