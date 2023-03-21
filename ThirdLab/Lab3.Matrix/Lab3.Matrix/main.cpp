#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#define DIM 12
#define N 20

using namespace std;

vector<vector<double>> E(DIM);

vector<vector<double>> __LoadMatrixFromFile(string fileName, int skip)
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
		if (skip != 0)
		{
			for (int i = 0; i < (skip + (skip / 12)); ++i)
			{
				in.ignore(250, '\n');
			}
		}

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

vector<double> __LoadVectorFromFile(string fileName, int skip)
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
		if (skip != 0)
		{
			for (int i = 0; i < (skip + (skip / 12)); ++i)
			{
				in.ignore(50, '\n');
			}
		}

		for (int i = 0; i < DIM; ++i)
		{
			in >> A[i];
		}

		in.close();
	}

	return A;
}

vector<vector<double>> __MMMultiplication(vector<vector<double>> x, vector<vector<double>> y)
{
	vector<vector<double>> xy(DIM);
	double el = 0.0;
	for (int i = 0; i < DIM; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			for (int k = 0; k < DIM; ++k)
			{
				el += x[i][k] * y[k][j];
			}

			xy[i].push_back(el);
			el = 0.0;
		}
	}

	return xy;
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

vector<double> __VectorPlusVector(vector<double> A, vector<double> B)
{
	vector<double> C;

	for (int i = 0; i < DIM; ++i)
	{
		C.push_back(A[i] + B[i]);
	}

	return C;
}

double __MatrixRate(vector<vector<double>> A)
{
	double curSum = 0.0;
	double rate = 0.0;

	for (int i = 0; i < DIM; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			curSum += abs(A[i][j]);
		}

		if (curSum > rate)
		{
			rate = curSum;
		}

		curSum = 0.0;
	}

	return rate;
}

double __VectorRate(vector<double> A)
{
	double rate = 0.0;

	for (int i = 0; i < DIM; ++i)
	{
		rate += A[i] * A[i];
	}

	rate = sqrt(rate);

	return rate;
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

void ActualErrorRate(vector<double> X, vector<double> exactX)
{
	double rate = 0.0;
	for (int i = 0; i < DIM; ++i)
	{
		rate += (X[i] - exactX[i]) * (X[i] - exactX[i]);
	}

	rate = sqrt(rate);

	ofstream out("ErrorRate.txt", std::ios_base::app);
	if (out.is_open())
	{
		out << rate << "; ";
	}
}

void DiscrepancyRate(vector<vector<double>> A, vector<double> B, vector<double> X)
{
	double rate = 0.0;
	vector<double> AX = __MVMultiplication(A, X);

	for (int i = 0; i < DIM; ++i)
	{
		rate += (AX[i] - B[i]) * (AX[i] - B[i]);
	}

	rate = sqrt(rate);

	ofstream out("DiscrepancyRate.txt", std::ios_base::app);
	if (out.is_open())
	{
		out << rate << "; ";
	}
}

void AFindXJacobi(vector<vector<double>> A, vector<double> B, vector<double> exactX, vector<double> X)
{
	vector<vector<double>> D(DIM);
	vector<vector<double>> DInv(DIM);

	vector<double> XPrev;

	double eps = 0.1;
	double rateC = 0.0;
	int iters = 0;

	for (int i = 0; i < DIM; ++i)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (i != j)
			{
				D[i].push_back(0);
				DInv[i].push_back(0);
			}
			else
			{
				D[i].push_back(A[i][i]);
				DInv[i].push_back(1 / D[i][i]);
			}
		}

		XPrev.push_back(X[i]);
	}

	vector<vector<double>> C = __MatrixMinusMatrix(E, __MMMultiplication(DInv, A));
	vector<double> G = __MVMultiplication(DInv, B);

	double stop = (1 - __MatrixRate(C)) / __MatrixRate(C);

	ofstream out("Results.txt", std::ios::out | std::ios::trunc);
	ofstream iOut("Iters.txt", std::ios::out | std::ios::trunc);

	out.precision(16);

	if (out.is_open() && iOut.is_open())
	{
		while (eps > 1e-13)
		{
			do 
			{
				for (int i = 0; i < DIM; ++i)
				{
					XPrev[i] = X[i];
				}

				X = __VectorPlusVector(__MVMultiplication(C, XPrev), G);

				iters++;
			} while (__VectorRate(__VectorMinusVector(X, XPrev)) >= stop * eps);

			iOut << iters << "; ";

			ActualErrorRate(X, exactX);
			DiscrepancyRate(A, B, X);

			for (int i = 0; i < DIM; i++)
			{
				out << X[i] << endl;
				X[i] = XPrev[i] = 0;
			}

			out << endl << iters << " " << eps << endl;

			out << endl << "___________________" << endl;

			eps *= 0.1;

			iters = 0;
		}
	}

	out.close();
}

void AFindXJacobiZ(vector<double> exactX, vector<double> X)
{
	vector<vector<double>> A(DIM);
	vector<double> B;

	vector<vector<double>> D(DIM);
	vector<vector<double>> DInv(DIM);

	vector<double> XPrev;

	double eps = 1e-2;
	double rateC = 0.0;
	int iters = 0;

	for (int i = 0; i < DIM; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			D[i].push_back(0);
			DInv[i].push_back(0);
		}

		XPrev.push_back(0);
	}

	ofstream out("Results.txt", std::ios::out | std::ios::trunc);
	ofstream iOut("Iters.txt", std::ios::out | std::ios::trunc);

	out.precision(16);

	if (out.is_open() && iOut.is_open())
	{
		for (int i = 0, skip = 0; i < N; ++i, skip += 12)
		{
			A = __LoadMatrixFromFile("A.txt", skip);
			B = __LoadVectorFromFile("B.txt", skip);

			for (int i = 0; i < DIM; ++i)
			{
				for (int j = 0; j < DIM; j++)
				{
					if (i == j)
					{
						D[i][i] = A[i][i];
						DInv[i][i] = 1 / D[i][i];
					}
				}
			}

			for (int i = 0; i < DIM; ++i)
			{
				X[i] = B[i];
			}

			vector<vector<double>> C = __MatrixMinusMatrix(E, __MMMultiplication(DInv, A));
			vector<double> G = __MVMultiplication(DInv, B);

			rateC = __MatrixRate(C);
			double stop = (1 - rateC) / rateC;

			do
			{
				for (int i = 0; i < DIM; ++i)
				{
					XPrev[i] = X[i];
				}

				X = __VectorPlusVector(__MVMultiplication(C, XPrev), G);

				iters++;
			} while (__VectorRate(__VectorMinusVector(X, XPrev)) >= stop * eps);

			iOut << iters << "; ";

			ActualErrorRate(X, exactX);
			DiscrepancyRate(A, B, X);

			for (int i = 0; i < DIM; i++)
			{
				out << X[i] << endl;
				X[i] = XPrev[i] = 0;
			}

			out << endl << iters << " " << eps << endl;

			out << endl << "___________________" << endl;

			iters = 0;
		}
	}

	out.close();
}

int main()
{
	srand(time(NULL));

	ofstream outER("ErrorRate.txt", std::ios::out | std::ios::trunc);
	outER.close();

	ofstream outDR("DiscrepancyRate.txt", std::ios::out | std::ios::trunc);
	outDR.close();

	cout.precision(16);

	vector<vector<double>> A(DIM);

	vector<double> B;
	vector<double> X;
	vector<double> exactX = __LoadVectorFromFile("X.txt", 0);

	A = __LoadMatrixFromFile("A.txt", 0);
	B = __LoadVectorFromFile("B.txt", 0);

	for (int i = 0; i < DIM; ++i)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (i != j)
			{
				E[i].push_back(0);
			}
			else
			{
				E[i].push_back(1);
			}
		}

		X.push_back(0);
	}

	AFindXJacobi(A, B, exactX, X);
	//AFindXJacobiZ(exactX, X);

	return 0;
}
