#include "ldltFactorization.h"

bool solve_TLDL(const std::vector<std::vector<double>> A, const std::vector<double> b, std::vector<double>& x){

  int matrixSize = A.size();
  if (matrixSize == 0 || A[0].size() != matrixSize || b.size() != matrixSize) {
    return false;
  }
  
  std::vector<std::vector<double>> a(matrixSize, std::vector<double>(matrixSize));
  std::vector<std::vector<double>> LT(matrixSize, std::vector<double>(matrixSize));
  std::vector<std::vector<double>> D(matrixSize, std::vector<double>(matrixSize)); 
  std::vector<std::vector<double>> L(matrixSize, std::vector<double>(matrixSize)); 
  
  for(int index = 1; index < matrixSize; index++){
    for(int jIndex = 0; jIndex <= index; jIndex++){
      double sum = 0;
      for(int kIndex = 0; kIndex <= jIndex - 1; kIndex ++)
        sum += a[index][kIndex] * L[jIndex][kIndex];
      a[index][jIndex] = A[index][jIndex] - sum;

      if(jIndex != 0){
        double sum = 0;
        for(int kIndex = 0; kIndex <= jIndex - 1; kIndex ++)
          sum += a[jIndex][kIndex] * L[jIndex][kIndex];
        D[jIndex][jIndex] = A[jIndex][jIndex] - sum;
      }else{
        D[jIndex][jIndex] = A[jIndex][jIndex];
      }

      L[index][jIndex] = a[index][jIndex] / D[jIndex][jIndex];
      L[jIndex][jIndex] = 1;
    }
  }

  LT = L;
  transposition(LT);

  // std::cout << "\n";
  // for(auto iter: a){
  //   for(auto jIter: iter)
  //     std::cout << jIter << " ";
  //   std::cout << "\n";
  // }

  std::cout << "\n";
  for(auto iter: D){
    for(auto jIter: iter)
      std::cout << jIter << " ";
    std::cout << "\n";
  }
  
  std::cout << "\n";
  
  for(auto iter: L){
    for(auto jIter: iter)
      std::cout << jIter << " ";
    std::cout << "\n";
  }
  
  std::cout << "\n";
  
  for(auto iter: LT){
    for(auto jIter: iter)
      std::cout << jIter << " ";
    std::cout << "\n";
  }

  std::vector<double> vectorY;
  Zeros(vectorY, matrixSize);

  //vectorY[0] = b[0];
  for(int index = 0; index < matrixSize; index++){
    double sum = 0.0;
    for(int kIndex = 0; kIndex <= index - 1; kIndex++)
      sum += L[index][kIndex] * vectorY[kIndex];
    vectorY[index] = b[index] - sum;
  }

  std::vector<double> vectorZ;
  Zeros(vectorZ, matrixSize);

  for(int index = 0; index < matrixSize; index++)
    vectorZ[index] = vectorY[index] / D[index][index];
  
  Zeros(x, matrixSize);
 
  //x[matrixSize - 1] = vectorZ[matrixSize - 1];

  for(int index = matrixSize - 1; index >= 0; index--){
    double sum = 0.0;
    for(int kIndex = index + 1; kIndex < matrixSize; kIndex++)
      sum += LT[index][kIndex] * x[kIndex];
    x[index] = vectorZ[index] - sum; 
  }

  

  return true;
}

void transposition(std::vector<std::vector<double>> &matrix){
  int matrixSize = matrix.size();

  for(int index = 0; index< matrixSize; index++)
    for(int jIndex = index; jIndex < matrixSize; jIndex++)
      std::swap(matrix[index][jIndex], matrix[jIndex][index]);
}

void Zeros(std::vector<double> &result, int &size){
    result.clear();
    for(int index = 0; index < size; index++)
        result.push_back(0);
}
/*

  std::cout << "\n";
  for(auto iter: D){
    for(auto jIter: iter)
      std::cout << jIter << " ";
    std::cout << "\n";
  }
  
  std::cout << "\n";
  
  for(auto iter: L){
    for(auto jIter: iter)
      std::cout << jIter << " ";
    std::cout << "\n";
  }


*/