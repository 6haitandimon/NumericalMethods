#include <iostream>
#include <fstream>
#include <iomanip>
#include "neutone/neuton.h"
#include "const.h"

double shag(double e, double hmax, double* uu)
{
	double h1, h2, h3, hmin;
	h1 = e / (abs(uu[0]) + (e / hmax));
	h2 = e / (abs(uu[1]) + (e / hmax));
	h3 = e / (abs(uu[2]) + (e / hmax));

	if (h1 < h2)
		hmin = h1;
	else
		hmin = h2;
	if (h3 < hmin)
		hmin = h3;

	return hmin;
}

double func(double* u, int n){
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
    double tau = 0, tauMax = 0.1;
    double tk = 0;
    double* yk = new double[n];

    for(int i = 0; i < n; i++)
        yk[i] = u[i];
    
    std::ofstream fout;
    fout.open("ExplicitEulerMethod.txt", std::ios::trunc);

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
            tmp[i] = func(yk, i);

        tau = shag(Eps, tauMax, tmp);

        for(int i = 0; i < n; i++){
            yk[i] += tau * tmp[i];
        }

        tk += tau;

        for(int i = 0; i < n; i++){ 
            if(i == 0){
            fout << tk << std::setw(15);
            fout << yk[i] << std::setw(15);
        }
        else if(i != n - 1){
            fout << yk[i] << std::setw(15);
        }
        else 
            fout << yk[i] << std::endl;
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
    double Tau = 1, Tau_minus = 1, Tau_plus = 1;
	double T = 1, TauMax = 0.01, TauMin = 0.01;
	double tk = 0, tk_plus;
	double Eps_k; 
	double* yk = new double[n];               
	double* yk_minus = new double[n];              
	double* yk_plus = new double[n];               
	Tau_minus = Tau = TauMin;
	
	for (int i = 0; i < n; i++) 
		yk[i] = yk_minus[i] = yk_plus[i] = u[i];

	std::ofstream out;
	out.open("EulerNejavny.txt");

	int sposob;
	std::cout << "1 - kvazioptimum, 2 - trjohzon" << std::endl;
	std::cin >> sposob;

	for (int i = 0; i < n; i++)
	{
		if (i == 0)
		{
			out << "    t" << std::setw(13);
			out << "    u" << i + 1 << std::setw(13);
		}
		else if (i == n - 1)
			out << "    u" << i + 1 << std::endl;
	}

	int kol = 0;
	do
	{
		do
		{
			tk_plus = tk + Tau;

			yk_plus = Neuton(yk_plus, yk, tk, Tau, n);
			for (int i = 0; i < n; i++)  
				Eps_k = -(Tau / (Tau + Tau_minus)) * (yk_plus[i] - yk[i] - Tau * (yk[i] - yk_minus[i]) / Tau_minus);
         	for (int i = 0; i < n; i++)  
			if (fabs(Eps_k) > Eps)
			{
				Tau /= 2;
				tk_plus = tk;
				for (int j = 0; j < n; j++)
					yk_plus[j] = yk[j];
			}
		}
		while (fabs(Eps_k) > Eps);
		
		if (sposob == 1)      
			Tau_plus = sqrt(Eps / abs(Eps_k)) * Tau;

		else if (sposob == 2)     
		{
			for (int i = 0; i < n; i++)  
			{
				if (fabs(Eps_k) > Eps)
					Tau_plus = Tau / 2;
				if ((Eps / 4 < fabs(Eps_k)) && (fabs(Eps_k) <= Eps))
					Tau_plus = Tau;
				if (fabs(Eps_k) <= Eps / 4)
					Tau_plus = 2 * Tau;
			}
		}
		else break;

		if (Tau_plus > TauMax) 
			Tau_plus = TauMax;

		for (int i = 0; i < n; i++) 
		{
			if (i == 0)
				out << tk << std::setw(15);
			out << std::setw(15) << yk[i];
			if (i == n - 1)
				out << std::endl;
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

	std::cout << std::endl << "Iterations quantity is " << kol << std::endl;
	out << std::endl << "Iterations quantity is " << kol << std::endl;
	out.close();
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

//jav
//1.00031        1.69082        1.94641      -0.268156



//neJavTrojnoy
//0.99875        1.14911        1.21583      -0.800407

//eps = 10-5
//0.998906        1.15345        1.22281      -0.811857

//kvazyopt
//0.994599        1.15315        1.22155      -0.794629

//eps = 10-5
//0.999912        1.15358        1.22312       -0.81622