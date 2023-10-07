#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

void Zeros(std::vector<double> &result, int &size){
    for(int index = 0; index < size; index++)
        result.push_back(0);
}

void MatrixMultiplication(std::vector<std::vector<double> > &A, std::vector<double> &x, std::vector<double> &matrixMultiplication){
    int matrixSize = A.size();
    Zeros(matrixMultiplication, matrixSize);
    for(int index = 0; index < matrixSize; index++){
        for(int jIndex = 0; jIndex < matrixSize; jIndex++)
            matrixMultiplication[index] += A[index][jIndex] *  x[jIndex];
    }
    return;
}

void ResidualVectorCalculation(std::vector<std::vector<double> > &A, std::vector<double> &b, std::vector<double> &result, std::vector<double> &residualVector){
    int matrixSize = b.size();
    std::vector<double> matrixMultiplication;
    MatrixMultiplication(A, result, matrixMultiplication);
    for(int index = 0; index < matrixSize; index++)
        residualVector.push_back(b[index] - matrixMultiplication[index]);
    return;
}

void GaussElimination(std::vector<std::vector<double> > &A, std::vector<double> &b, std::vector<double> &result){
    int matrixSize = A.size();
        for(int index = 0; index < matrixSize; index++){
            int indexMax = index;
            for(int jIndex = index + 1; jIndex < matrixSize; jIndex++){
                if(abs(A[jIndex][index]) > abs(A[indexMax][index])){
                    indexMax = jIndex;
                }
            }

            std::swap(A[index], A[indexMax]);
            std::swap(b[index], b[indexMax]);

            for(int jIndex = index + 1; jIndex < matrixSize; jIndex++){
                double factor = A[jIndex][index] / A[index][index];
                for(int kIndex = index; kIndex < matrixSize; kIndex++)
                    A[jIndex][kIndex] -= factor * A[index][kIndex];
                b[jIndex] -= factor * b[index];
            }
        }

    Zeros(result, matrixSize);

    for(int index = matrixSize - 1; index > -1; index--){
        double sum = 0;
        for(int jIndex = index; jIndex < matrixSize; jIndex++){
            sum += A[index][jIndex] * result[jIndex];
        }
        result[index] = (b[index] - sum) / A[index][index];
    }

    return;
}

int main(){
    std::vector<std::vector<double> > A = {
                            {2.30, 5.70, -0.80},
                            {3.50, -2.70, 5.30},
                            {1.70, 2.30, -1.80}
                            };
    std::vector<double> b = {-6.49, 19.20, -5.09};
    std::vector<double> result;

    std::vector<std::vector<double> > A1(A);
    std::vector<double> b1(b);

    std::vector<double> residualVector;

    GaussElimination(A, b, result);
    ResidualVectorCalculation(A1, b1, result, residualVector);

    std::cout << "Решение СЛАУ: ";
    for(auto iter : result)
        std::cout << iter << " ";
    
    std::cout << "\nВектор невязки: ";
    for(auto iter : residualVector)
        std::cout << iter << " ";

    double maxValue = 0;
    for(auto iter : residualVector)
        if(abs(iter) > maxValue)
            maxValue = iter;

    std::cout << "\nНорма вектора невязки: " << maxValue;

    return 0;
}