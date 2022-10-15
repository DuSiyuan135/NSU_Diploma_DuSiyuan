#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;

double Wide_w(double x, double w0, double cL, double L);
double Velocity_u(double x, double w, double Q);
double AccurateLocation(double t, double w0, double cL, double L, double Q);
double Euler(double tau, double x, double w0, double cL, double L, double Q);
double RungeKuttaFourthOrder(double tau, double x, double w0, double cL, double L, double Q);

int main()
{
	double u0 = 0.1, a = 0, b = 100, Q = 0.1, w0 = 0.1, cL = 0.2, L = 100;

	cout << "Euler:" << endl;
	for (double tau = 2; tau >= 0.5; tau /= 2)
	{
		cout << "tau= " << tau << endl;
		double t = 0, err = 0;
		// printf("t=%.5lf\t x_E=%.5lf\t x=%.5lf\t err=%.5lf\t\n", t, a, AccurateLocation(t, w0, cL, L, Q), 0);
		for (double x = a; x < b;)
		{
			x = Euler(tau, x, w0, cL, L, Q);
			t += tau;
			if (err < abs(AccurateLocation(t, w0, cL, L, Q) - x))
				err = abs(AccurateLocation(t, w0, cL, L, Q) - x);
			 printf("t=%.5lf\t x_E=%.5lf\t x=%.5lf\t err=%.5lf\t", t, x, AccurateLocation(t, w0, cL, L, Q), abs(AccurateLocation(t, w0, cL, L, Q) - x));
			 cout << endl;
		}
		printf("t=%.5lf\t err=%.7lf\t", t, err);
		cout << endl;
	}

	cout << "RungeKutta:" << endl;
	for (double tau = 2; tau >= 0.5; tau /= 2)
	{
		cout << "tau= " << tau << endl;
		double t = 0, err = 0;
		// printf("t=%.5lf\t x_R=%.5lf\t x=%.5lf\t err=%.5lf\t\n", t, a, AccurateLocation(t, w0, cL, L, Q), 0);
		for (double x = a; x < b;)
		{
			x = RungeKuttaFourthOrder(tau, x, w0, cL, L, Q);
			t += tau;
			if (err < abs(AccurateLocation(t, w0, cL, L, Q) - x))
				err = abs(AccurateLocation(t, w0, cL, L, Q) - x);
			 printf("t=%.5lf\t x_R=%.5lf\t x=%.5lf\t err=%.5lf\t", t, x, AccurateLocation(t, w0, cL, L, Q), abs(AccurateLocation(t, w0, cL, L, Q) - x));
			 cout << endl;
		}
		printf("t=%.5lf\t err=%.7lf\t", t, err);
		cout << endl;
	}

	system("pause");

	return 0;
}

double Wide_w(double x, double w0, double cL, double L)
{
	return w0*(1 - (1 - cL)*x / L);
}

double Velocity_u(double x, double w, double Q)
{
	return Q/w;
}

double AccurateLocation(double t, double w0, double cL, double L, double Q)
{
	return (L/(1 - cL) + sqrt(w0*L*(2*cL*Q*t - 2*Q*t + L*w0)) / ((cL - 1)*w0));
	// return (125 - 5 * sqrt(625 - 10 * t));
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