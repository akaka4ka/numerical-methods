#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

vector<double> operator-(vector<double> v, vector<double> u)
{
	if (v.size() != u.size())
	{
		return vector<double>();
	}

	vector<double> res;
	for (int i = 0; i < v.size(); ++i)
	{
		res.push_back(v[i] - u[i]);
	}

	return res;
}

double f(double x, double y)
{
	return (-y / (x + 1) - y * y);
}

double fExact(double x)
{
	return (1 / ((x + 1) * log(x + 1)));
}

int iters;
double error;
double anger = 0;
int grid = 16;

void Eiler(double a, double b, double eps)
{
	int n = 1;
	double h = (b - a) / grid;
	iters = 0;
	int curIters = 0;

	vector<double> gridX (grid + 1);
	vector<double> gridY (grid + 1);

	for (int i = 0; i < grid + 1; ++i)
	{
		gridX[i] = a + h * i;
	}

	gridY[0] = 1 / (2 * log(2)) + anger;

	vector<double> x;
	vector<double> y;

	double yLast = 0;

	double errMax = 0;

	for (int i = 1; i < grid + 1; ++i)
	{
		int n1 = 1;
		double h1 = (gridX[i] - gridX[i - 1]) / n1;

		curIters = 0;

		y.push_back(0);

		do
		{
			yLast = y[y.size() - 1];

			h1 = (gridX[i] - gridX[i - 1]) / n1;

			x.clear();
			y.clear();

			x.push_back(gridX[i - 1]);
			y.push_back(gridY[i - 1]);

			for (int i = 1; i < 2 * n1 + 1; ++i)
			{
				x.push_back(x[i - 1] + h1 / 2);
			}

			for (int i = 0; i < 2 * n1 - 1; i += 2)
			{
				y.push_back(y[i] + h1 * f(x[i], y[i]) / 2);
				y.push_back(y[i] + h1 * f(x[i + 1], y[i + 1]));
			}

			n1 *= 2;
			++curIters;

		} while (abs(y[y.size() - 1] - yLast) / (2.0 * 2 - 1) > eps);

		if (iters < curIters)
		{
			iters = curIters;
		}

		gridY[i] = y[y.size() - 1];
	}

	double curError = 0;

	for (int i = 0; i < gridX.size(); ++i)
	{
		if (curError < abs(fExact(gridX[i]) - gridY[i]))
		{
			curError = abs(fExact(gridX[i]) - gridY[i]);
		}
	}

	error = curError;
}

int main()
{
	double a = 1.0;
	double b = 5.0;

	cout.precision(16);

	double eps = 0.1;

	#pragma region DefaultEilerWithEps

	/*ofstream out("out.txt");

	for (int i = 0; i < 12; ++i, eps *= 0.1)
	{
		Eiler(a, b, eps);
		out << eps << " " << iters << " " << error << endl;
	}*/

	#pragma endregion

	anger = 1;
	vector<double> p;
	vector<double> g;

	for (size_t i = 0; i < 12; ++i, anger *= 0.1)
	{
		Eiler(a, b, 1e-3);
		p.push_back(anger);
		g.push_back(error);
	}

	cout << endl;

	for (auto _p : p)
	{
		cout << _p << " ";
	}
	cout << endl;

	cout << endl;
	for (auto _g : g)
	{
		cout << _g << " ";
	}
	cout << endl;

	#pragma region HGraphics

	/*int n = 1;
	double h = (b - a) / grid;
	iters = 0;
	int curIters = 0;

	vector<double> gridX(grid + 1);
	vector<double> gridY(grid + 1);

	for (int i = 0; i < grid + 1; ++i)
	{
		gridX[i] = a + h * i;
	}

	gridY[0] = 1 / (2 * log(2));

	vector<double> x;
	vector<double> y;

	vector<double> err;

	double h1;

	for (int i = 1; i < grid + 1; ++i)
	{
		int n1 = 16;
		h1 = (gridX[i] - gridX[i - 1]) / n1;

		x.clear();
		y.clear();

		x.push_back(gridX[i - 1]);
		y.push_back(gridY[i - 1]);

		for (int i = 1; i < 2 * n1 + 1; ++i)
		{
			x.push_back(x[i - 1] + h1 / 2);
		}

		for (int i = 0; i < 2 * n1 - 1; i += 2)
		{
			y.push_back(y[i] + h1 * f(x[i], y[i]) / 2);
			y.push_back(y[i] + h1 * f(x[i + 1], y[i + 1]));
		}

		gridY[i] = y[y.size() - 1];
	}

	for (int i = 0; i < gridX.size(); ++i)
	{
		err.push_back(abs(fExact(gridX[i]) - gridY[i]));
	}

	std::cout << "x: " << endl;

	for (auto e : gridX)
	{
		std::cout << e << " ";
	}

	std::cout << endl << endl;

	std::cout << "y: " << endl;

	for (auto e : gridY)
	{
		std::cout << e << " ";
	}

	std::cout << endl << endl;

	std::cout << "Error: " << endl;

	for (auto e : err)
	{
		std::cout << e << " ";
	}

	std::cout << endl << endl;

	std::cout << "h: " << h1;*/

	#pragma endregion

	return 0;
}
