#include <vector>

bool solve_TLDL(const std::vector<std::vector<double>> A, const std::vector<double> b, std::vector<double>& x){

  int n = A.size();
  if (n == 0 || A[0].size() != n || b.size() != n) {
    return false;
  }
  

  std::vector<std::vector<double>> T(n, std::vector<double>(n));
  std::vector<std::vector<double>> L(n, std::vector<double>(n)); 
  std::vector<std::vector<double>> D(n, std::vector<double>(n)); 
  
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      double sum = 0.0;
      for (int k = 0; k < i; k++) {
        sum += T[k][i] * L[k][j] * D[k][k];
      }
      T[i][j] = (A[i][j] - sum) / A[i][i];
    }
    
    for (int j = i + 1; j < n; j++) {
      double sum = 0.0;
      for (int k = 0; k < i; k++) {
        sum += T[k][j] * L[k][i] * D[k][k];
      }
      L[i][j] = (A[j][i] - sum) / A[i][i];
    }
    
    double sum = 0.0;
    for (int k = 0; k < i; k++) {
      sum += T[k][i] * L[k][i] * D[k][k];
    }
    D[i][i] = A[i][i] - sum;
    

    if (D[i][i] == 0.0) {
      return false; 
    }
  }
  
  
  std::vector<double> y(n);
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < i; j++) {
      sum += T[j][i] * y[j];
    }
    y[i] = b[i] - sum;
  }
  
  
  std::vector<double> z(n);
  for (int i = 0; i < n; i++) {
    z[i] = y[i] / D[i][i];
  }
  
 
  std::vector<double> w(n);
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < i; j++) {
      sum += L[j][i] * w[j];
    }
    w[i] = z[i] - sum;
  }
  
  
  x.resize(n);
  for (int i = n - 1; i >= 0; i--) {
    double sum = 0.0;
    for (int j = i + 1; j < n; j++) {
      sum += T[i][j] * x[j];
    }
    x[i] = w[i] - sum;
  }
  

  return true;
}