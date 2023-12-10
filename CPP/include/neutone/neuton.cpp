#include "neuton.h"

std::vector<std::vector<double>> get_J(double x1, double x2, double k)
{
    return std::vector<std::vector<double>>{{(func1(x1 + x1 * k, x2) - func1(x1, x2)) / k / x1, (func1(x1, x2 + x2 * k) - func1(x1, x2)) / k / x2},
                                            {(func2(x1 + x1 * k, x2) - func2(x1, x2)) / k / x1, (func2(x1, x2 + x2 * k) - func2(x1, x2)) / k / x2}};
}

void Neuton(double x1, double x2, double E0, double E1, int max_iter, double k)
{
    double delta0 = std::max(fabs(func1(x1, x2)),
                             fabs(func2(x1, x2)));
    double delta1 = 1;

    int iter = 0;
   
    while ((delta0 > E0 || delta1 > E1) && iter < max_iter)
    {
        iter++;
        std::cout << iter << ":\t"
                  << "x1: " << x1 << "\t"
                  << "x2: " << x2 << std::endl;

        std::vector<double> F{-func1(x1, x2), -func2(x1, x2)};
        std::vector<std::vector<double>> J = get_J(x1, x2, k);

        std::vector<double> solution;
        GaussElimination(J, F, solution);
        x1 += solution[0];
        x2 += solution[1];

        delta0 = fabs(F[0]);
        for (int i = 1; i < F.size(); i++)
        {
            if (delta0 < fabs(F[i]))
                delta0 = fabs(F[i]);
        }

        double v1 = fabs(x1) < 1 ? solution[0] : solution[0] / x1;
        double v2 = fabs(x2) < 1 ? solution[1] : solution[1] / x2;
        
        delta1 = std::max(v1, v2);
    }
}