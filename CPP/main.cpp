#include <iostream>
#include <vector>
#include <cmath>

double Func(double x){
    return pow(x + pow(x, 3.0), 0.5);
}

double TrapezoidMethod(double a, double b){
    double square = 0;
    double h = (b - a) / 10.0;
    double sum = 0;
    for (int i = 1; i < 10; i++){
        sum += 2.0 * Func(a + i * h);
    }
    sum += Func(a) + Func(b);
    square = h / 2.0 * sum;

    return square;
}

double SimpsonMethod(double a, double b, double eps){
    
    int n = 2;
    
    double h = double((b - a) / n);
    int count = 0;
    double Ih1 = 1;
    double Ih2 = 0;
    
    while(abs(Ih1 - Ih2) > 15.0 * eps){
        Ih1 = Ih2;
        h = double((b - a) / n);

        double square = 0;
        double sum = Func(a) + Func(b);
        
        int m = n / 2;

        for (int i = 1; i <= m; i++){
            sum += 4 * Func(a + (2 * i - 1) * h);
        }

        for(int i = 1; i < m - 1; i++){
            sum += 2 * Func(a + 2 * i * h);
        }

        square = h / 3.0 * sum;
        Ih2 = square;
        n *= 2;
        count++;
    }

    return Ih2;
}

int main(){
    double a = 0.6;
    double b = 1.742;

    double E1 = 1e-4;
    double E2 = 1e-5;

    std::cout << "Метод трапеции: " <<TrapezoidMethod(a, b) << std::endl;
    std::cout << "Метод Симпсона при eps = 1e-4: " << SimpsonMethod(a, b, E1) << std::endl;
    std::cout << "Метод Симпсона при eps = 1e-5: " << SimpsonMethod(a, b, E2) << std::endl;
    return 0;
}