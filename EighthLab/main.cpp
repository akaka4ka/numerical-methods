#define _CRT_SECURE_NOWARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>

#define LEFT 1.0
#define RIGHT 2.0
#define MAX_POWER 12
#define TEST_N 16

using namespace std;

double yFunc(double x)
{
	return pow(x, 0.5);
}

double dyFunc(double z)
{
	return z;
}

double d2yFunc_U(double x, double y, double z)
{
	return (pow(x, 0.5) - 2 * y - (2 - x) * z) / (2 * x * (x + 2));
}

double d2yFunc_VW(double x, double y, double z)
{
	return (-2 * y - (2 - x) * z) / (2 * x * (x + 2));
}

void ModEuler(double a, double b, int n, double y0, double z0, double(*Func)(double, double, double), vector<double>& yValue)
{
	double h = (b - a) / n;
	double tmpY = 0, tmpZ = 0;
	double x = a, y = y0, z = z0;
	yValue[0] = y;

	for (int i = 1; i <= n; i++)
	{
		tmpY = y + h * dyFunc(z) / 2.0;
		tmpZ = z + h * Func(x, y, z) / 2.0;

		y = y + h * dyFunc(tmpZ);
		z = z + h * Func(x + h / 2.0, tmpY, tmpZ);

		yValue[i] = y;

		x = x + h;
	}

	return;
}

int main(void)
{
	std::cout.precision(16);

	int n = 8;
	double c1 = 0;
	double c2 = 0;
	double a = LEFT;
	double b = RIGHT;
	double h = (b - a) / n;
	double x = 0;
	double max = 0;
	
	#pragma region TwoStepValues

	vector<double> result;
	vector<double> error;
	
	vector<double> u(n + 1);
	vector<double> v(n + 1);
	vector<double> w(n + 1);

	ModEuler(a, b, n, 0, 0, d2yFunc_U, u);
	ModEuler(a, b, n, 1, 0, d2yFunc_VW, v);
	ModEuler(a, b, n, 0, 1, d2yFunc_VW, w);

	c1 = a;
	c2 = (yFunc(b) - u[n] - c1 * v[n]) / w[n];

	cout << "X: " << endl;
	for (int i = 0; i <= n; i++)
	{
		x = LEFT + i * h;
		cout << x << " ";
		result.push_back(u[i] + c1 * v[i] + c2 * w[i]);
		error.push_back(abs(result[i] - yFunc(x)));
	}

	cout << endl << endl;

	cout << "Y: " << endl;
	for (auto e : result)
	{
		cout << e << " ";
	}

	cout << endl << endl;

	cout << "Error: " << endl;
	for (auto e : error)
	{
		cout << e << " ";
	}

	#pragma endregion

	#pragma region ByStep

	vector<double> error;
	vector<double> delta;
	vector<double> step;

	for (int j = 1; j <= MAX_POWER; j++)
	{
		//n = pow(2, j);
		double result = 0;
		n = TEST_N;
		double _delta = 0;
		_delta = pow(0.1, j);
		h = (b - a) / n;
		max = 0;

		vector<double> u(n + 1);
		vector<double> v(n + 1);
		vector<double> w(n + 1);

		ModEuler(a, b, n, 0, 0, d2yFunc_U, u);
		ModEuler(a, b, n, 1, 0, d2yFunc_VW, v);
		ModEuler(a, b, n, 0, 1, d2yFunc_VW, w);

		c1 = a + _delta;
		c2 = (yFunc(b) - u[n] - c1 * v[n]) / w[n];

		for (int i = 0; i <= n; i++)
		{
			result = u[i] + c1 * v[i] + c2 * w[i];

			if (fabs(result - yFunc(a + i * h)) > max)
			{
				max = fabs(result - yFunc(a + i * h));
			}
		}

		error.push_back(max);
		delta.push_back(_delta);
		step.push_back(h);
	}

	std::cout << "Delta:" << endl;
	for (auto e : delta)
	{
		std::cout << e << " ";
	}

	std::cout << endl << endl;

	std::cout << "Error:" << endl;
	for (auto e : error)
	{
		std::cout << e << " ";
	}

	std::cout << endl << endl;

	std::cout << "Step:" << endl;
	for (auto e : step)
	{
		std::cout << e << " ";
	}

	#pragma endregion
}