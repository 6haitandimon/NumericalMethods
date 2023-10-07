#include <iostream>
#include <cmath>
#include <vector>

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
    std::vector<double> result = {0.0, 0.0, 0.0};

    std::vector<std::vector<double> > A1(A);
    std::vector<double> b1(b);

    GaussElimination(A, b, result);

    std::cout << "Решение СЛАУ: ";
    for(auto iter : result)
        std::cout << iter << " ";

    return 0;
}