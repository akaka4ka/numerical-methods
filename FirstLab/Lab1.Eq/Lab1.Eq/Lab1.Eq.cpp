#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cmath>

using namespace std;

double Pow(double x, int n)
{
    double _result = 1.0;

    for (int i = 0; i < n; ++i)
    {
        _result *= x;
    }

    return _result;
}

double Polynom(double x)
{
    return 3 * Pow(x, 4) + 4 * Pow(x, 3) - 12 * Pow(x, 2) - 5;
}

double PolynomDer(double x)
{
    return 12 * Pow(x, 3) + 12 * Pow(x, 2) - 24 * x;
}

double Trans(double x)
{
    return log(x) + Pow((x + 1), 3);
}

double TransDer(double x)
{
    return 1 / x + 3 * Pow((x + 1), 2);
}

void BisectionMethod(double left, double right, double(*Func)(double))
{
    double eps = 0.1;
    double a = left;
    double b = right;
    double c = 0;
    double x = 0;

    int iters = 0;

    if (Func(a) * Func(b) > 0)
    {
        cout << "Чот не то, давай по новой" << endl;
        return;
    }

    FILE* filePtr = fopen("bisectionResult.txt", "w");
    if (filePtr == nullptr)
    {
        cout << "Не удалось открыть файл для вывода..." << endl;
        return;
    }

    while (eps > 1e-15)
    {
        a = left;
        b = right;
        iters = 0;

        while (abs(b - a) > 2 * eps)
        {
            iters++;

            c = (a + b) / 2;

            if (Func(a) * Func(c) < 0)
            {
                b = c;
            }
            else
            {
                a = c;
            }
        }

        x = (a + b) / 2;

        //fprintf(filePtr, "Root: %.16F, epsilon: %.15F, iterations: %d\n", x, eps, iters);
        //fprintf(filePtr, "Root: %.16F, epsilon: %e, iterations: %d\n", x, eps, iters);
        //fprintf(filePtr, "%.15F %d\n", eps, iters);
        //fprintf(filePtr, "%.16F; ", fabs(Func(x)));
        //fprintf(filePtr, "%d; ", iters);

        eps *= 0.1;
    }

    fclose(filePtr);
}

void NewtonMethod(double x, double exactX, double(*Func)(double), double(*FuncDer)(double))
{
    double eps = 1e-5;
    double start = x;
    int iters = 0;

    FILE* filePtr = fopen("newtonResult.txt", "w");
    FILE* fileXPtr = fopen("newtonResultX.txt", "w");
    if ((filePtr == NULL) || (fileXPtr == NULL))
    {
        cout << "Не удалось открыть файл для вывода..." << endl;
        return;
    }

    while (start > exactX + 2 * eps)
    {
        x = start;
        iters = 0;

        while (fabs(Func(x)) > eps)
        {
            x = x - Func(x) / FuncDer(x);
            iters++;
        }

        //fprintf(filePtr, "Root: %.16F, epsilon: %.15F, iterations: %d\n", x, eps, iters);
        //fprintf(filePtr, "Root: %.16F, epsilon: %e, iterations: %d\n", x, eps, iters);
        //fprintf(filePtr, "%.15F %d\n", eps, iters);
        //fprintf(filePtr, "%.16F; ", fabs(Func(x)));
        //fprintf(filePtr, "%d; ", iters);
        //fprintf(fileXPtr, "%llf; ", start);

        start = (exactX + 5 * start) / 6;
    }

    fclose(filePtr);
}

int main()
{
    double exactXP = 1.592088;
    double exactXT = 0.187439;

    //BisectionMethod(0.01, 2, Polynom);
    //BisectionMethod(0.01, 2, Trans);
    
    //NewtonMethod(10000, exactXP, Polynom, PolynomDer);
    //NewtonMethod(10000, exactXT, Trans, TransDer);

    return 0;
}
