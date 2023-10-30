

bool solve_TLDL(const std::vector<std::vector<double>> A, const std::vector<double> b, std::vector<double>& x){

  int n = A.size();
  if (n == 0 || A[0].size() != n || b.size() != n) {
    return false;
  }
  
  // Выделяем память для матриц T, L и D
  std::vector<std::vector<double>> T(n, std::vector<double>(n)); // Матрица T - верхнетреугольная с единицами на диагонали
  std::vector<std::vector<double>> L(n, std::vector<double>(n)); // Матрица L - нижнетреугольная с единицами на диагонали
  std::vector<std::vector<double>> D(n, std::vector<double>(n)); // Матрица D - диагональная
  
  // Выполняем T L D L-факторизацию матрицы A
  for (int i = 0; i < n; i++) {
    // Вычисляем элементы матрицы T в i-й строке
    for (int j = i + 1; j < n; j++) {
      double sum = 0.0;
      for (int k = 0; k < i; k++) {
        sum += T[k][i] * L[k][j] * D[k][k];
      }
      T[i][j] = (A[i][j] - sum) / A[i][i];
    }
    
    // Вычисляем элементы матрицы L в i-м столбце
    for (int j = i + 1; j < n; j++) {
      double sum = 0.0;
      for (int k = 0; k < i; k++) {
        sum += T[k][j] * L[k][i] * D[k][k];
      }
      L[i][j] = (A[j][i] - sum) / A[i][i];
    }
    
    // Вычисляем элемент матрицы D в i-й позиции
    double sum = 0.0;
    for (int k = 0; k < i; k++) {
      sum += T[k][i] * L[k][i] * D[k][k];
    }
    D[i][i] = A[i][i] - sum;
    
    // Проверяем, что элемент не равен нулю
    if (D[i][i] == 0.0) {
      return false; // Матрица A не имеет T L D L-факторизации
    }
  }
  
  // Решаем СЛАУ T y = b с помощью прямого хода метода Гаусса
  std::vector<double> y(n); // Вектор y размера n
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < i; j++) {
      sum += T[j][i] * y[j];
    }
    y[i] = b[i] - sum;
  }
  
  // Решаем СЛАУ D z = y с помощью деления на диагональные элементы
  std::vector<double> z(n); // Вектор z размера n
  for (int i = 0; i < n; i++) {
    z[i] = y[i] / D[i][i];
  }
  
  // Решаем СЛАУ L w = z с помощью прямого хода метода Гаусса
  std::vector<double> w(n); // Вектор w размера n
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < i; j++) {
      sum += L[j][i] * w[j];
    }
    w[i] = z[i] - sum;
  }
  
  // Решаем СЛАУ T^T x = w с помощью обратного хода метода Гаусса
  x.resize(n); // Изменяем размер вектора x на n
  for (int i = n - 1; i >= 0; i--) {
    double sum = 0.0;
    for (int j = i + 1; j < n; j++) {
      sum += T[i][j] * x[j];
    }
    x[i] = w[i] - sum;
  }
  

  return true;
}