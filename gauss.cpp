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

// int main(){
//     std::vector<std::vector<double> > A = {
//                             {6, 13, -17},
//                             {13, 29, -38},
//                             {-17, -38, 50}
//                             };
//     std::vector<double> b = {2, 4, -5};

//     std::vector<double> result;
//     std::vector<double> resultAx;
//     std::vector<double> resultLDLT;

//     std::vector<double> residualVector;
//     std::vector<double> matrixAXMultiplication;

//     double maxRelativeError = 0.0;
//     int indexError = 0;
    
//     GaussElimination(A, b, result);
//     ResidualVectorCalculation(A, b, result, residualVector);

//     std::cout << "Решение СЛАУ: ";
//     for(auto iter : result)
//         std::cout << std::setprecision(18) << iter << " ";
    
//     std::cout << "\nВектор невязки: ";
//     for(auto iter : residualVector)
//         std::cout << iter << " ";

//     double maxValue = 0;
//     for(auto iter : residualVector)
//         if(abs(iter) > maxValue)
//             maxValue = iter;

//     std::cout << "\nНорма вектора невязки: " << maxValue;

//     MatrixMultiplication(A, result, matrixAXMultiplication);
//     GaussElimination(A, matrixAXMultiplication, resultAx);
    
//     RelativeError(resultAx, result, maxRelativeError, indexError);

//     std::cout << "\n\nРешение СЛАУ Ax = A~x: ";
//     for(auto iter : resultAx)
//         std::cout << iter << " ";
    
//     std::cout << "\nОтносительная погрешность метода Гаусса: " << maxRelativeError;
//     std::cout << "\nМаксимальная разница была дастигнута в элементе: " << indexError;
    
// }

