#include <iostream>
#include <fstream>
#include <iomanip>
#include "gauss/gauss.h"
#include "neutone/neuton.h"
#include "const.h"

double func(double* u, double t, int n){
    switch (n){
        case 0:
            return (k - a) / a * u[1] * u[2];
        case 1:
            return (k + a) / k * u[0] * u[2];
        case 2:
            return (a - k) / k * u[0] * u[1];
    }
    return 0;
}

void ExplicitEulerMethod(double* u, int n){
    double tau = 0, tauMax = 0.01;
    double tk = 0;
    double* yk = new double[n];

    for(int i = 0; i < n; i++)
        yk[i] = u[i];
    
    std::ofstream fout;
    fout.open("data.txt", std::ios::trunc);

    for(int i = 0; i < n; i++){
        if(i == 0){
            fout << "   t" << std::setw(20);
            fout << "   u" << i + 1 << std::setw(10);
        }
        else if(i != n - 1){
            fout << "   u" << i + 1 << std::setw(25);
        }
        else 
            fout <<  "   u" << i + 1 << std::endl;
    }
    int kol = 0;
    do{
        double *tmp = new double[n];
        for(int i = 0; i < n; i++)
            tmp[i] = func(yk, tk, i);
        
        tau = eps1 / (fabs(tmp[0]) + eps1 / tauMax);
        for(int i = 1; i <= n; i++){
            double tauIntermediate = eps1 / (fabs(tmp[i - 1]) + eps1 / tauMax);
            if(tauIntermediate < tau)
                tau = tauIntermediate;
        }

        for(int i = 0; i < n; i++){
            yk[i] += tau * tmp[i];
        }

        tk += tau;

        for(int i = 0; i < n; i++){ 
            if(i != n - 1)
                fout << std::setw(15) << tk << std::setw(15) << yk[i];
            else
                fout << "\n";
        }
        kol++;
        delete[] tmp;
    }while(tk < T);

    fout << std::endl << "Iteration is quantity is " << kol << std::endl;

    std::cout << "Done\n";

    delete[] yk;

    fout.close();
}
void NoExplicitEulerMethod(double* u, int n){
    	double Tau, Tau_minus, Tau_plus;
	double T = 1, TauMax = 0.01, TauMin = 0.01;
	double tk = 0, tk_plus;
	double Eps_k; 
	double* yk = new double[n];               
	double* yk_minus = new double[n];              
	double* yk_plus = new double[n];               
	Tau_minus = Tau = TauMin;

	for (int i = 0; i < n; i++) 
		yk[i] = yk_minus[i] = yk_plus[i] = u[i];

	std::ofstream fout;
	fout.open("EulerNejavny.txt");

	int sposob;
	std::cout << "1 - kvazioptimum, 2 - trjohzon" << std::endl;
	std::cin >> sposob;

	for (int i = 0; i < n; i++)
	{
		if (i == 0)
		{
			fout << "    t" << std::setw(13);
			fout << "    u" << i + 1 << std::setw(13);
		}
		else if (i == n - 1)
			fout << "    u" << i + 1 << std::endl;
	}

	int kol = 0;
	do
	{
		do
		{
			tk_plus = tk + Tau;
			yk_plus = Neuton(yk_plus, yk, tk, Tau, n);

			for (int k = 0; k < n; k++)
				Eps_k = -(Tau / (Tau + Tau_minus)) * (yk_plus[k] - yk[k] - Tau * (yk[k] - yk_minus[k]) / Tau_minus);

			for (int k = 0; k < n; k++) 
			if (fabs(Eps_k) > eps2)
			{
				Tau /= 2;
				tk_plus = tk;
				for (int j = 0; j < n; j++)
					yk_plus[j] = yk[j];
			}
		}
		while (fabs(Eps_k) > eps2);
		
		if (sposob == 1)
			Tau_plus = sqrt(eps2 / abs(Eps_k)) * Tau;

		else if (sposob == 2) 
		{
			for (int i = 0; i < n; i++)  
			{
				if (fabs(Eps_k) > eps2)
					Tau_plus = Tau/2;
				if ((eps2 / 4 < fabs(Eps_k)) && (fabs(Eps_k) <= eps2))
					Tau_plus = Tau;
				if (fabs(Eps_k) <= eps2 / 4)
					Tau_plus = 2 * Tau;
			}
		}
		else break;

		if (Tau_plus > TauMax) 
			Tau_plus = TauMax;

		for (int i = 0; i < n; i++) 
		{
			if (i == 0)
				fout << tk << std::setw(15);
			fout << std::setw(15) << yk[i];
			if (i == n - 1)
				fout << "\n";
		}

		for (int i = 0; i < n; i++)
		{
			yk_minus[i] = yk[i];
			yk[i] = yk_plus[i];
		}
		Tau_minus = Tau;
		Tau = Tau_plus;
		tk = tk_plus;

		kol++;
	} while(tk < T);

    delete[] yk;
    delete[] yk_minus;
    delete[] yk_plus;

	std::cout << std::endl << "Iterations quantity is " << kol << std::endl;
	fout << std::endl << "Iterations quantity is " << kol << std::endl;
	fout.close();
}

int main(){
    int n = 3;
    double* u = new double[n];

    std::cout << "start value of U[0] = 1, U[1] = 1, U[2] = 1\n";
    u[0] = 1;  
    u[1] = 1;  
    u[2] = 1;

    int methods;
    std::cout << "1 - Explicit\n 2 - No Explicit\n";
    std::cin >> methods;
    switch (methods){
        case 1:
            std::cout << "Explicit Euler Method\n";
            ExplicitEulerMethod(u, n);
            break;
        case 2:
            std::cout << "No Explicit Euler Method\n";
            NoExplicitEulerMethod(u, n);
            break;
    default:
        std::cout << "No method\n";
        break; 
    }

    delete[] u;
    return 0;
}