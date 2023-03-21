#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

const double integr = -5537379.0 / 250000.0;

double func(double x)
{
	return (pow(x, 5) - 3.2*pow(x, 3) + 1.5*pow(x, 2) - 7*x - 9.5);
}

void ThreeEighths(ofstream& out, double a, double b, double eps)
{
	int n = 1;
	int powerN = 0;
	double h = (b - a) / (n * 3.0);

	double curIntegr = 3 * h / 8 * (func(a) + 3 * (func(a + h) + func(a + 2 * h)) + func(b));
	int p = 4;
	double rungeDown = (pow(2, p) - 1);
	double curPoint = a;
	double prevIntegr = curIntegr;

	do
	{
		prevIntegr = curIntegr;
		curIntegr = 0.0;
		n *= 2;
		++powerN;
		h = (b - a) / (n * 3.0);
		curPoint = a;
		for (int i = 0; i < n; i++) {
			curIntegr += (func(curPoint) + 3 * (func(curPoint + h) + func(curPoint + 2 * h)) + func(curPoint + 3 * h));
			curPoint += 3 * h;
		}
		curIntegr *= 3 * h / 8;
	} while (abs(curIntegr - prevIntegr) / rungeDown >= eps);

	if (out.is_open())
	{
		out << eps << " " << abs(curIntegr - integr) << " " << powerN << " " << h << endl;
	}
}

int main()
{
	ofstream out("out.txt", std::ios::out | std::ios::trunc);
	ofstream& outRef = out;

	double eps[12] = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 };

	for (auto cur : eps)
	{
		ThreeEighths(outRef, -2.5, 1.3, cur);
	}

	out.close();

	return 0;
}
