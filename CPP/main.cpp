#include <iostream>
#include <vector>
#include "gauss.hpp"

void Powers(std::vector<int> vector, int m, std::vector<int> &powerX){
    int vectorSize = vector.size();
    for(int i = 1; i <= 2 * m; i++){
        int sum = 0;
        for(int j = 0; j < vectorSize; j++){
            sum += pow(vector[j], i);
        }
        powerX.push_back(sum);
    }
    return;
}

void Sum(std::vector<int> powerX, int m, std::vector<std::vector<double>>& sumX){
    for(int l = 0; l < m + 1; l++){
        for(int j = 1; j < m + 1; j++)
            sumX[l][j] = powerX[l + j - 1];
    }
    return;
}

void Praw(std::vector<int> X, std::vector<int> Y, int m, std::vector<double> &praw){
    int vectorSize = X.size();
    for(int l = 0; l < m + 1; l++){
        int sum = 0;
        int k1 = l;
        for(int i = 0; i < vectorSize; i++)
            sum += Y[i] * pow(X[i], k1);
        
        praw.push_back(sum);
    }
    return;
}

void LeastSquareMethod(std::vector<int> X, std::vector<int> Y, std::vector<double> &coefficient, int N, int m){
    std::vector<int> powerX;
    std::vector<std::vector<double>> sumX(m + 1, std::vector<double>(m + 1));
    std::vector<double> praw;

    Powers(X, m, powerX);

    Sum(powerX, m, sumX);
    sumX[0][0] = N;

    Praw(X, Y, m, praw);

    GaussElimination(sumX, praw, coefficient);
    return;
}

double ResidualVariance(std::vector<int> X, std::vector<int> Y, std::vector<double> coefficient, int m, int N){
    double residualVariance = (1.0 / (N - m - 1.0));
    int vectorSize = X.size();
    double sum = 0;
    for(int i = 0; i < vectorSize; i++){
        sum += pow(Y[i] - coefficient[0] - coefficient[1] * X[i] - coefficient[2] * pow(X[i], 2), 2);
    }

    residualVariance *= sum;

    return residualVariance;
}

double ApproximationError(std::vector<int> X, std::vector<int> Y, std::vector<double> coefficient, int N){
    double approximationError = 1.0 / N;
    int vectorSize = X.size();
    double sum = 0;

    for(int i = 0; i < vectorSize; i++)
        sum = abs((Y[i] - (coefficient[0] + coefficient[1] * X[i] + coefficient[2] * pow(X[i], 2.0))) / Y[i]);
    approximationError *= sum * 100.0;
    return approximationError;
}

int main(){
    std::vector<int> X = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<int> Y = {3, 87, 156, 210, 238, 252, 239, 211, 158, 90, -5};
    std::vector<double> coefficient;

    int N = 11;
    int m = 2;
    double residualVariance = 0;

    LeastSquareMethod(X, Y, coefficient, N, m);
    residualVariance = ResidualVariance(X, Y, coefficient, m, N);

    std::cout << "Coefficients: ";
    for(auto iter: coefficient)
        std::cout << iter << " ";

    std::cout << "\nResidual Variance: " << sqrt(residualVariance);
    //std::cout << "\nApproximation Error: " << ApproximationError(X, Y, coefficient, N);

    std::cout << "\n";

    int vectorSize = X.size();

    // std::cout << vectorSize << "\n";

    for(int i = 0; i < vectorSize; i++)
        std::cout << "x = " << X[i] << ", y = " << coefficient[0] + coefficient[1] * X[i] + coefficient[2] * pow(X[i], 2.0) << "\n";

    return 0;
}