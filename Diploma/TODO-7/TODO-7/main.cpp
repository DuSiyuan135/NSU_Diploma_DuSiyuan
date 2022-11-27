#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string> 

#define PI 3.1415926536


using namespace std;

double Velocity_u(double r, double w, double Q);
double AccurateLocation(double t, double R_in, double w, double Q);
double Euler(double tau, double r, double w, double Q);
double RungeKuttaFourthOrder(double tau, double r, double w, double Q);


int main()
{
	double Q = 0.1, w = 0.1, R_in = 1, R_out = 10;

	cout << "Euler:" << endl;
	for (double tau = 2; tau >= 0.5; tau /= 2)
	{
		cout << "tau= " << tau << "\t";
		double t = 0, err = 0, r_E = 0;
		r_E = AccurateLocation(t, R_in, w, Q);

		string site = "tau=" + to_string(tau);
		ofstream OutFile("Result\\Euler\\" + site + ".txt");
		OutFile << "tau=" << tau << endl;
		OutFile << "t=" << t << " \t" << "r=" << R_in << " \t" << "x_E=" << r_E << " \t" << "Err=" << err << endl;
		// printf("t=%.5lf\t x_E=%.5lf\t x=%.5lf\t err=%.5lf\t\n", t, a, AccurateLocation(t, w0, cL, L, Q), 0);
		for (double r = R_in; r_E < R_out;)
		{
			r = Euler(tau, r, w, Q);
			t += tau;
			r_E = AccurateLocation(t, R_in, w, Q);
			if (err < abs(r_E - r))
				err = abs(r_E - r);
			// printf("t=%.5lf\t x_E=%.5lf\t x=%.5lf\t err=%.5lf\t", t, x, AccurateLocation(t, w0, cL, L, Q), abs(AccurateLocation(t, w0, cL, L, Q) - x));
			// cout << endl;
			OutFile << "t=" << t << " \t" << "r=" << r << " \t" << "r_E=" << r_E << " \t" << "Err=" << err << endl;
		}
		printf("err=%.10lf\t t=%.5lf\t", err, t);
		cout << endl;

		OutFile.close();
	}

	cout << "RungeKutta:" << endl;
	for (double tau = 2; tau >= 0.5; tau /= 2)
	{
		cout << "tau= " << tau << "\t";
		double t = 0, err = 0, r_E = 0;
		r_E = AccurateLocation(t, R_in, w, Q);

		string site = "tau=" + to_string(tau);
		ofstream OutFile("Result\\RungeKutta\\" + site + ".txt");
		OutFile << "tau=" << tau << endl;
		OutFile << "t=" << t << " \t" << "x=" << R_in << " \t" << "r_E=" << r_E << " \t" << "Err=" << err << endl;
		// printf("t=%.5lf\t x_E=%.5lf\t x=%.5lf\t err=%.5lf\t\n", t, a, AccurateLocation(t, w0, cL, L, Q), 0);
		for (double r = R_in; r_E < R_out;)
		{
			r = RungeKuttaFourthOrder(tau, r, w, Q);
			t += tau;
			r_E = AccurateLocation(t, R_in, w, Q);
			if (err < abs(r_E - r))
				err = abs(r_E - r);
			// printf("t=%.5lf\t x_E=%.5lf\t x=%.5lf\t err=%.5lf\t", t, x, AccurateLocation(t, w0, cL, L, Q), abs(AccurateLocation(t, w0, cL, L, Q) - x));
			// cout << endl;
			OutFile << "t=" << t << " \t" << "r=" << r << " \t" << "r_E=" << r_E << " \t" << "Err=" << err << endl;
		}
		printf("err=%.10lf\t t=%.5lf\t", err, t);
		cout << endl;

		OutFile.close();
	}

	system("pause");

	return 0;
}

double Velocity_u(double r, double w, double Q)
{
	return (Q / (2 * PI*w*r));
}

double AccurateLocation(double t, double R_in, double w, double Q)
{
	return (sqrt(Q*t / (PI*w) + R_in*R_in));
}

double Euler(double tau, double r, double w, double Q)
{
	double u = Velocity_u(r, w, Q);

	return (r + tau*u);
}

double RungeKuttaFourthOrder(double tau, double r, double w, double Q)
{
	double k1, k2, k3, k4;

	k1 = Velocity_u(r, w, Q);
	k2 = Velocity_u(r + tau / 2 * k1, w, Q);
	k3 = Velocity_u(r + tau / 2 * k2, w, Q);
	k4 = Velocity_u(r + tau * k3, w, Q);

	return (r + tau / 6 * (k1 + 2 * k2 + 2 * k3 + k4));
}