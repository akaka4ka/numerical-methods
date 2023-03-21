#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#define DIM 12
#define N 100

using namespace std;

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

void __LYB(vector<vector<double>> L, vector<double> B, vector<double>& Y)
{
	double numerator = 0.0;
	for (int i = 0; i < DIM; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			numerator += L[i][j] * Y[j];
		}
		numerator = B[i] - numerator;

		Y[i] = numerator / L[i][i];

		numerator = 0.0;
	}
}

void __UXY(vector<vector<double>> U, vector<double> Y, vector<double>& X)
{
	double numerator = 0.0;
	for (int i = DIM - 1; i >= 0; --i)
	{
		for (int j = i + 1; j < DIM; ++j)
		{
			numerator += U[i][j] * X[j];
		}
		numerator = Y[i] - numerator;

		X[i] = numerator / U[i][i];

		numerator = 0.0;
	}
}

void LU(vector<vector<double>> A, vector<vector<double>>& L, vector<vector<double>>& U)
{
	double sumU = 0.0;
	double sumL = 0.0;

	for (int i = 0; i < DIM; ++i)
	{
		L[i][i] = 1;

		for (int j = 0; j < DIM; ++j)
		{
			if (i <= j)
			{
				for (int k = 0; k < i; ++k)
				{
					sumU += L[i][k] * U[k][j];
				}

				U[i][j] = A[i][j] - sumU;
				sumU = 0.0;
			}
			
			if (i > j)
			{
				for (int k = 0; k < j; ++k)
				{
					sumL += L[i][k] * U[k][j];
				}

				L[i][j] = (A[i][j] - sumL) / U[j][j];
				sumL = 0.0;
			}
		}
	}
}

void AFindXLU(vector<vector<double>> A, vector<double> B, vector<double>& X)
{
	vector<vector<double>> L(DIM);
	vector<vector<double>> U(DIM);

	vector<double> Y;

	for (int i = 0; i < DIM; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			L[i].push_back(0);
			U[i].push_back(0);
		}

		Y.push_back(0);
	}

	LU(A, L, U);

	__LYB(L, B, Y);

	__UXY(U, Y, X);
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

void Check(vector<vector<double>> A, vector<double> exactX, vector<double> B, int condA)
{
	int n = 15;

	vector<double> dX;
	vector<double> dB;
	vector<double> newB;
	vector<double> newX;

	double mag = 1e-3;
	double absDX = 0.0;
	double absExactX = 0.0;
	double absDB = 0.0;
	double absB = 0.0;

	ofstream outL("checkL.txt", std::ios::out | std::ios::trunc);
	ofstream outR("checkR.txt", std::ios::out | std::ios::trunc);
	outL.precision(16);
	outR.precision(16);

	for (int i = 0; i < DIM; ++i)
	{
		dX.push_back(0);
		//dB.push_back((1 + rand() % 2) * mag);
		dB.push_back(mag);
		newB.push_back(B[i] + dB[i]);
		newX.push_back(0);
	}

	for (int i = 0; i < n; i++)
	{
		AFindXLU(A, newB, newX);

		for (int i = 0; i < DIM; ++i)
		{
			dX[i] = newX[i] - exactX[i];

			absDX += dX[i] * dX[i];
			absExactX += exactX[i] * exactX[i];
			absDB += dB[i] * dB[i];
			absB += B[i] * B[i];
		}

		absDX = sqrt(absDX);
		absExactX = sqrt(absExactX);
		absDB = sqrt(absDB);
		absB = sqrt(absB);

		if (outL.is_open() && outR.is_open())
		{
			outL << absDX / absExactX << "; ";
			outR << absDB / absB << "; ";
		}

		mag *= 10;

		for (int i = 0; i < DIM; i++)
		{
			dX[i] = 0;
			dB[i] = mag;
			newB[i] = B[i] + dB[i];
			newX[i] = 0;
		}
	}

	outL.close();
	outR.close();

	/*cout << condA << ": " << check << endl;

	cout << endl;*/
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
	vector<double> exactX = __LoadVectorFromFile("X.txt", 0);;

	for (int i = 0; i < DIM; ++i)
	{
		X.push_back(0);
	}

	/*for (int i = 0, skip = 0, cond = 100; i < N; ++i, skip += 12, cond += 100)
	{
		A = __LoadMatrixFromFile("A.txt", skip);
		B = __LoadVectorFromFile("B.txt", skip);
		AFindXLU(A, B, X);
		ActualErrorRate(X, exactX);
		DiscrepancyRate(A, B, X);
		__PrintVector(X);
		Check(A, X, B, cond);
	}*/

	A = __LoadMatrixFromFile("A.txt", 0);
	B = __LoadVectorFromFile("B.txt", 0);

	Check(A, exactX, B, 100);

	/*__PrintMatrix(A);*/

	/*__PrintVector(B);*/

	/*__PrintVector(X);*/

	return 0;
}
