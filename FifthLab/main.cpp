#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>

using namespace std;

double f(double x)
{
	return (x * x) - 1 - log(x);
}

double fder(double x)
{
	return 2 * x - 1 / x;
}

double Hermit(double x, double* gridX, double* gridY, double* gridYDer, int n)
{
	double hermit = 0.0;
	double sum = 0.0;
	double mult = 1.0;

	for (int i = 0; i < n; ++i)
	{
		mult = 1.0;
		sum = 0.0;

		for (int j = 0; j < n; ++j)
		{
			if (j != i)
			{
				sum += (x - gridX[i]) / (gridX[i] - gridX[j]);
			}
		}

		for (int k = 0; k < n; ++k)
		{
			if (k != i)
			{
				mult *= ((x - gridX[k]) / (gridX[i] - gridX[k])) * ((x - gridX[k]) / (gridX[i] - gridX[k]));
			}
		}

		hermit += ((x - gridX[i]) * gridYDer[i] + (1 - 2 * sum) * gridY[i]) * mult;
	}

	return hermit;
}

double* ChebGrid(int n, double a, double b)
{
	double* grid = new double[n];

	double first = (b + a) / 2;
	double second = (b - a) / 2;

	double pi = M_PI;

	for (int i = 0; i < n; ++i)
	{
		grid[i] = first + second * cos(pi * (2 * i + 1) / (2 * n + 2));
	}

	return grid;
}

int main()
{
	cout.precision(16);

	//int n = 3;

	double a = 1.0;
	double b = 5.0;
	double sup = 0.0;

	double* x = new double[21];
	double* err = new double[21];

	cout << "X: " << endl;
	for (int i = 0; i < 21; ++i)
	{
		x[i] = 1.0 + 0.2 * i;
		cout << x[i] << ", ";
	}

	cout << endl;
	for (int n = 1; n <= 50; ++n)
	{
		double* chebGridX = ChebGrid(n, a, b);
		double* chebGridY = new double[n];
		double* chebGridYDer = new double[n];

		for (int i = 0; i < n; ++i)
		{
			chebGridY[i] = f(chebGridX[i]);
			chebGridYDer[i] = fder(chebGridX[i]);
		}

		for (int i = 0; i < 21; ++i)
		{
			err[i] = abs(f(x[i]) - Hermit(x[i], chebGridX, chebGridY, chebGridYDer, n));
		}

		sup = err[0];
		for (int i = 0; i < 21; ++i)
		{
			if (sup < err[i])
			{
				sup = err[i];
			}
		}

		cout << sup << ", ";

		delete[n] chebGridX;
		delete[n] chebGridY;
		delete[n] chebGridYDer;
	}

	

	/*double* x = new double[21];

	cout << "X: " << endl;
	for (int i = 0; i < 21; ++i)
	{
		x[i] = 1.0 + 0.2 * i;
		cout << x[i] << ", ";
	}*/

	/*cout << endl << "Y: " << endl;
	for (int i = 0; i < 21; ++i)
	{
		cout << Hermit(x[i], chebGridX, chebGridY, chebGridYDer, n) << ", ";
	}*/

	/*cout << endl << "XNodes: " << endl;
	for (int i = n - 1; i >= 0; --i)
	{
		cout << chebGridX[i] << ", ";
	}

	cout << endl << "YNodes: " << endl;
	for (int i = n - 1; i >= 0; --i)
	{
		cout << chebGridY[i] << ", ";
	}*/

	//for (int i = 0; i < n; ++i)
	//{
	//	//cout << chebGridX[i] << " " << chebGridY[i] << " " << chebGridYDer[i] << endl;
	//	cout << chebGridX[i] << " " << chebGridY[i] << " " << Hermit(chebGridX[i], chebGridX, chebGridY, chebGridYDer, n) << endl;
	//}

	return 0;
}