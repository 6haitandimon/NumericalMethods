#include "gauss.h"

void Zeros(std::vector<double> &result, int &size){
    result.clear();
    for(int index = 0; index < size; index++)
        result.push_back(0);
}

void MatrixMultiplication(std::vector<std::vector<double> > A, std::vector<double> x, std::vector<double> &matrixMultiplication){
    int matrixSize = A.size();
    Zeros(matrixMultiplication, matrixSize);
    for(int index = 0; index < matrixSize; index++){
        for(int jIndex = 0; jIndex < matrixSize; jIndex++)
            matrixMultiplication[index] += A[index][jIndex] * x[jIndex];
    }
    return;
}


void GaussElimination(std::vector<std::vector<double> > A, std::vector<double> b, std::vector<double> &result){
    int matrixSize = A.size();
        for(int index = 0; index < matrixSize; index++){
            int indexMax = index;
            for(int jIndex = index + 1; jIndex < matrixSize; jIndex++){
                if(abs(A[jIndex][index]) > abs(A[indexMax][index])){
                    indexMax = jIndex;
                }
            }
            assert(A[indexMax][index] != 0);
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
