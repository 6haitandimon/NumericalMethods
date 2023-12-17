#include "neuton.h"

double** createMatrix(int x)
{
	double **A = new double*[x];
	for (int i = 0; i < x; i++)
		A[i] = new double[x + 1];
	return A;
}
void deleteMatrix(double **X, int x)
{
	for (int i = 0; i < x; i++)
		delete X[i];
	delete[] X;
}
void copyMatrix(double **X, double **copyX, int x)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x + 1; j++)
			copyX[i][j] = X[i][j];
	}
}


double f1(double *uk1, double *uk, double t, double Tau)
{
	return uk1[0] - uk[0] - Tau * (uk1[1] * uk1[2] * (k - a) / a);
}

double f2(double *uk1, double *uk, double t, double Tau)
{
	return uk1[1] - uk[1] - Tau * (uk1[0] * uk1[2] * (k + a) / k);
}
double f3(double *uk1, double *uk, double t, double Tau)
{
	return uk1[2] - uk[2] - Tau * (uk1[0] * uk1[1] * (a - k) / a);
}


double Differential(pf f, double *uk1, double *uk, double t, double Tau, int n)
{
	double dx = 1e-9;
	double* D = new double[3];
	for (int i = 0; i < 3; i++)
	{
		if (i == n - 1)
		{
			D[i] = uk1[i] + dx;
			i++;
		}
		D[i] = uk1[i];
	}

	double F = f(uk1, uk, t, Tau);
	double dF = f(D, uk, t, Tau);

	delete[] D;
	return (dF - F) / dx;
}


double* Neuton(double *yk_plus, double *yk, double tk, double Tau, int n)
{
    double **Jako = createMatrix(n);
	double *An = new double[n];
	double e = 1e-9, b1 = 0, b2 = 0;

	int itr = 1;

	do{
		Jako[0][0] = Differential(f1, yk_plus, yk, tk, Tau, 1);
		Jako[0][1] = Differential(f1, yk_plus, yk, tk, Tau, 2);
		Jako[0][2] = Differential(f1, yk_plus, yk, tk, Tau, 3);
		Jako[0][3] = -f1(yk_plus, yk, tk, Tau);
		Jako[1][0] = Differential(f2, yk_plus, yk, tk, Tau, 1);
		Jako[1][1] = Differential(f2, yk_plus, yk, tk, Tau, 2);
		Jako[1][2] = Differential(f2, yk_plus, yk, tk, Tau, 3);
		Jako[1][3] = -f2(yk_plus, yk, tk, Tau);
		Jako[2][0] = Differential(f3, yk_plus, yk, tk, Tau, 1);
		Jako[2][1] = Differential(f3, yk_plus, yk, tk, Tau, 2);
		Jako[2][2] = Differential(f3, yk_plus, yk, tk, Tau, 3);
		Jako[2][3] = -f3(yk_plus, yk, tk, Tau);

		double **copyF = createMatrix(n);					
		copyMatrix(Jako, copyF, n);
		if (!GaussElimination(An, copyF, n))				
		{
			deleteMatrix(copyF, n);	
			deleteMatrix(Jako, n);
			delete []An;
			return yk_plus;
		}

		yk_plus[0] += An[0];	
		yk_plus[1] += An[1];
		yk_plus[2] += An[2];

		if(fabs(f1(yk_plus, yk, tk, Tau)) > b1)
			b1 = fabs(f1(yk_plus, yk, tk, Tau));
		if(fabs(f2(yk_plus, yk, tk, Tau)) > b1)
			b1 = fabs(f2(yk_plus, yk, tk, Tau));
		if(fabs(f3(yk_plus, yk, tk, Tau)) > b1)
			b1 = fabs(f3(yk_plus, yk, tk, Tau));


		for (int i = 0; i < n; i++)									
		{
			if (fabs(An[i]) < 1)
				b2 = fabs(An[i]);
			else if (fabs(An[i]) >= 1)
				b2 = fabs(An[i] / yk_plus[i]);
		}

		deleteMatrix(copyF, n);
		itr++;
	} while ((b1 > e || b2 > e) && (itr < iter_max));

	deleteMatrix(Jako, n);		
	delete[] An;

	return yk_plus;
}