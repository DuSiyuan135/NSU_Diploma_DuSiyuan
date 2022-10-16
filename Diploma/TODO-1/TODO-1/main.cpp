#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string> 


using namespace std;

double Velocity_u(double x, double u0);
double AccurateLocation(double t, double u0, double a);
double Euler(double tau, double u0, double x);
double RungeKuttaFourthOrder(double tau, double u0, double x);


int main()
{
	double u0 = 0.1, a = 0, b = 100;
	
	cout << "Euler:" << endl;
	for (double tau = 2; tau >= 0.5; tau /= 2)
	{
		cout << "tau= " << tau << "\t";
		double t = 0, err = 0, x_E = 0;
		x_E = AccurateLocation(t, u0, a);

		string site = "tau=" + to_string(tau);
		ofstream OutFile("Result\\Euler\\" + site + ".txt");
		OutFile << "tau=" << tau << endl;
		OutFile << "t=" << t << " \t" << "x=" << a << " \t" << "x_E=" << x_E << " \t" << "Err=" << err << endl;
		// printf("t=%.5lf\t x=%.5lf\t x_E=%.5lf\t err=%.5lf\t\n", t, a, AccurateLocation(t, u0, a), 0);

		for (double x = a; x_E < b;)
		{
			x = Euler(tau, u0, x);
			t += tau;
			x_E = AccurateLocation(t, u0, a);
			if (abs(x - x_E) > err)
			{
				err = abs(x - x_E);
				if (err<1e-7)
				{
					err = 0;
				}
			}
			// printf("t=%.5lf\t x=%.5lf\t x_E=%.5lf\t err=%.5lf\t", t, x, AccurateLocation(t, u0, a), abs(x-AccurateLocation(t, u0, a)));
			// cout << endl;
			OutFile << "t=" << t << " \t" << "x=" << x << " \t" << "x_E=" << x_E << " \t" << "Err=" << err << endl;
		}
		printf("err=%.5lf\t t=%.5lf\t", err, t);
		cout << endl;

		OutFile.close();
	}

	cout << "RungeKutta:" << endl;
	for (double tau = 2; tau >= 0.5; tau /= 2)
	{
		cout << "tau= " << tau << "\t";
		double t = 0, err = 0, x_E = 0;
		x_E = AccurateLocation(t, u0, a);

		string site = "tau=" + to_string(tau);
		ofstream OutFile("Result\\RungeKutta\\" + site + ".txt");
		OutFile << "tau=" << tau << endl;
		OutFile << "t=" << t << " \t" << "x=" << a << " \t" << "x_E=" << x_E << " \t" << "Err=" << err << endl;
		// printf("t=%.5lf\t x=%.5lf\t x_E=%.5lf\t err=%.5lf\t\n", t, a, AccurateLocation(t, u0, a), 0);
		
		for (double x = a; x_E < b;)
		{
			x = RungeKuttaFourthOrder(tau, u0, x);
			t += tau;
			x_E = AccurateLocation(t, u0, a);
			if (abs(x - x_E) > err)
			{
				err = abs(x - x_E);
				if (err<1e-7)
				{
					err = 0;
				}
			}
			// printf("t=%.5lf\t x=%.5lf\t x_E=%.5lf\t err=%.5lf\t", t, x, AccurateLocation(t, u0, a), abs(x - AccurateLocation(t, u0, a)));
			// cout << endl;
			OutFile << "t=" << t << " \t" << "x=" << x << " \t" << "x_E=" << x_E << " \t" << "Err=" << err << endl;
		}
		printf("err=%.5lf\t t=%.5lf\t", err, t);
		cout << endl;

		OutFile.close();
	}

	system("pause");

	return 0;
}

double Velocity_u(double x, double u0)
{
	return u0;
}

double AccurateLocation(double t, double u0, double a)
{
	return (u0*t + a);
}

double Euler(double tau, double u0, double x)
{
	double u = Velocity_u(x, u0);

	return (x + tau*u);
}

double RungeKuttaFourthOrder(double tau, double u0, double x)
{
	double k1, k2, k3, k4;

	k1 = Velocity_u(x, u0);
	k2 = Velocity_u(x + tau / 2 * k1, u0);
	k3 = Velocity_u(x + tau / 2 * k2, u0);
	k4 = Velocity_u(x + tau * k3, u0);
	
	return (x+tau/6*(k1+2*k2+2*k3+k4));
}
