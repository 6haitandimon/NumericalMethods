#include "neuton.h"

double f1(double *uk1, double *uk, double t, double Tau)
{
	return uk1[0] - uk[0] - Tau*(-uk1[0] * uk1[1] + ((t < 1e-9) ? 0.0 : (sin(t) / t)));
}

double f2(double *uk1, double *uk, double t, double Tau)
{
	return uk1[1] - uk[1] - Tau*(-uk1[1] * uk1[1] + (3.125*t) / (1 + t*t));
}


double Differential(pf f, double *uk1, double *uk, double t, double Tau, int n)
{
	double dx = 1e-9;
	double* D = new double[2];
	for (int i = 0; i < 2; i++)
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
    std::vector<std::vector<double>> Jako(n, std::vector<double>(n));
	std::vector<double> An(n);
	double e = 1e-9, b1 = 0, b2 = 0;

	int itr = 1;

	do{
		Jako[0][0] = Differential(f1, yk_plus, yk, tk, Tau, 1);
		Jako[0][1] = Differential(f1, yk_plus, yk, tk, Tau, 2);
		Jako[0][2] = -f1(yk_plus, yk, tk, Tau);				
		Jako[1][0] = Differential(f2, yk_plus, yk, tk, Tau, 1);
		Jako[1][1] = Differential(f2, yk_plus, yk, tk, Tau, 2);
		Jako[1][2] = -f2(yk_plus, yk, tk, Tau);					

		std::vector<std::vector<double>> copyF(n, std::vector<double>(n));
        std::vector<double> result;
		copyF = Jako;

        

		if (!GaussElimination(copyF, An, result))
            return 0;
		yk_plus[0] += An[0];
		yk_plus[1] += An[1];

		if (fabs(f1(yk_plus, yk, tk, Tau)) > fabs(f2(yk_plus, yk, tk, Tau)))
			b1 = fabs(f1(yk_plus, yk, tk, Tau));
		else
			b1 = fabs(f2(yk_plus, yk, tk, Tau));

		for (int i = 0; i < n; i++)
		{
			if (fabs(An[i]) < 1)
				b2 = fabs(An[i]);
			else if (fabs(An[i]) >= 1)
				b2 = fabs(An[i] / yk_plus[i]);
		}

		itr++;
	} while ((b1 > e || b2 > e) && (itr < iter_max));

	return yk_plus;
}