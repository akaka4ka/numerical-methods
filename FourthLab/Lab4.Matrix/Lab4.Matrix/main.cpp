#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#define DIM 12
#define N 100

using namespace std;

vector<vector<double>> E(DIM);

vector<vector<double>> __LoadMatrixFromFile(string fileName)
{
	vector<vector<double>> A(DIM);

	for (int i = 0; i < DIM; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			A[i].push_back(0);
		}
	}

	string str;

	ifstream in(fileName);
	if (in.is_open())
	{
		for (int i = 0; i < DIM; ++i)
		{
			for (int j = 0; j < DIM; ++j)
			{
				in >> A[i][j];
			}
		}

		in.close();
	}

	return A;
}

vector<double> __LoadVectorFromFile(string fileName)
{
	vector<double> A;

	for (int i = 0; i < DIM; ++i)
	{
		A.push_back(0);
	}

	string str;

	ifstream in(fileName);
	if (in.is_open())
	{
		for (int i = 0; i < DIM; ++i)
		{
			in >> A[i];
		}

		in.close();
	}

	return A;
}

double __VVScalar(vector<double> x, vector<double> y)
{
	double res = 0.0;

	for (int i = 0; i < DIM; i++)
	{
		res += x[i] * y[i];
	}

	return res;
}

vector<double> __MVMultiplication(vector<vector<double>> x, vector<double> y)
{
	vector<double> xy;
	double el = 0.0;
	for (int i = 0; i < DIM; ++i)
	{
		for (int k = 0; k < DIM; ++k)
		{
			el += x[i][k] * y[k];
		}

		xy.push_back(el);
		el = 0.0;
	}

	return xy;
}

vector<vector<double>> __MatrixMinusMatrix(vector<vector<double>> A, vector<vector<double>> B)
{
	vector<vector<double>> C(DIM);

	for (int i = 0; i < DIM; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			C[i].push_back(A[i][j] - B[i][j]);
		}
	}

	return C;
}

vector<double> __VectorMinusVector(vector<double> A, vector<double> B)
{
	vector<double> C;

	for (int i = 0; i < DIM; ++i)
	{
		C.push_back(A[i] - B[i]);
	}

	return C;
}

void __PrintMatrix(vector<vector<double>> m)
{
	for (int i = 0; i < DIM; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			cout << setw(10) << m[i][j] << " |";
		}

		cout << endl;
	}

	cout << endl;
}

void __PrintVector(vector<double> v)
{
	for (int i = 0; i < DIM; ++i)
	{
		cout << setw(10) << v[i] << " |";

		cout << endl;
	}

	cout << endl;
}

vector<double> NormVector(vector<double> v)
{
	double max = 0.0;

	for (int i = 0; i < DIM; ++i)
	{
		if (abs(v[i]) > max)
		{
			max = abs(v[i]);
		}
	}

	for (int i = 0; i < DIM; ++i)
	{
		v[i] = v[i] / max;
	}

	return v;
}

double VectorNorm(vector<double> v)
{
	double norm = 0.0;

	for (int i = 0; i < DIM; ++i)
	{
		norm += v[i] * v[i];
	}

	return sqrt(norm);
}

vector<double> __VectorOnNum(vector<double> v, double u)
{
	vector<double> v1(DIM);
		
	for (int i = 0; i < DIM; ++i)
	{
		v1[i] = v[i] * u;
	}

	return v1;
}

void PowerMethod(vector<vector<double>> A, vector<double> exactX)
{
	vector<double> Xn(DIM);

	Xn[0] = 1;

	vector<double> Xn1(DIM);
	vector<double> diff(DIM);
	vector<double> AX(DIM);
	vector<double> LX(DIM);
	double num = 0;
	double den = 0;
	double Lprev = 0;
	double L = 0;

	double eps = 0.1;
	double iters = 0;

	ofstream outOn("Abs.txt", std::ios::out | std::ios::trunc);
	outOn.precision(16);
	ofstream outTw("VectorNorm.txt", std::ios::out | std::ios::trunc);
	outTw.precision(16);
	ofstream outTh("AXNorm.txt", std::ios::out | std::ios::trunc);
	outTh.precision(16);
	ofstream outFo("Vector.txt", std::ios::out | std::ios::trunc);
	outFo.precision(16);
	ofstream outFi("Iters.txt", std::ios::out | std::ios::trunc);

	while (eps > 1e-13)
	{
		do
		{
			Lprev = L;

			iters++;

			Xn1 = __MVMultiplication(A, Xn);

			double max = 0.0;
			int it = 0;
			for (int i = 0; i < DIM; i++)
			{
				if (abs(Xn[i]) > max)
				{
					max = abs(Xn[i]);
					it = i;
				}
			}

			L = Xn1[it] / Xn[it];

			for (int i = 0; i < DIM; i++)
			{
				Xn[i] = Xn1[i];
			}

			Xn = NormVector(Xn);

		} while (abs(L - Lprev) > eps);

		Xn1 = NormVector(Xn1);

		double diff = exactX[0] / Xn1[0];

		for (int i = 0; i < DIM; ++i)
		{
			Xn1[i] *= diff;
		}

		AX = __MVMultiplication(A, Xn1);

		for (int i = 0; i < DIM; i++)
		{
			LX[i] = Xn1[i] * L;
		}

		cout << L << ", ";
		outOn << abs(12 - L) << ", ";
		outTw << VectorNorm(__VectorMinusVector(Xn1, exactX)) << ", ";
		outTh << VectorNorm(__VectorMinusVector(AX, LX)) << ", ";

		outFo << eps << endl;
		for (int i = 0; i < DIM; i++)
		{
			outFo << Xn1[i] << endl;
		}
		outFo << " " << endl;

		outFi << iters << ", ";

		iters = 0;
		L = 0;

		for (int i = 0; i < DIM; ++i)
		{
			Xn[i] = 0;
		}

		Xn[0] = 1;

		eps *= 0.1;
	}

	/*cout << "L = " << L + mu << endl;
	cout << "Iters = " << iters << endl;
	cout << "Actual eps = " << abs(0.1 - (L + mu)) << endl;
	__PrintVector(NormVector(Xn1));*/
}

int main()
{
	srand(time(NULL));

	ofstream outER("ErrorRate.txt", std::ios::out | std::ios::trunc);
	if (outER.is_open())
	{
		outER.close();
	}

	ofstream outDR("DiscrepancyRate.txt", std::ios::out | std::ios::trunc);
	if (outDR.is_open())
	{
		outDR.close();
	}

	cout.precision(16);

	vector<vector<double>> A(DIM);
	vector<double> X;

	A = __LoadMatrixFromFile("A.txt");
	X = __LoadVectorFromFile("X.txt");

	//__PrintMatrix(A);

	PowerMethod(A, X);

	return 0;
}
