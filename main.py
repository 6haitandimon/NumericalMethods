import copy

def zeros(n):
    x = []
    for i in range(n):
        x.append(0)
    return x

def matrix_multiplication(a, b):
    lenA = len(a)
    result_vector = zeros(lenA)
    for i in range(lenA):
        for j in range(lenA):
            result_vector[i] += a[i][j] * b[j]
    return result_vector


def residual_vector_calculation(A, x, b): 
    intermediate_vector = matrix_multiplication(A, x)
    residual_vector = []
    for i in range(len(intermediate_vector)):
        residual_vector.append(intermediate_vector[i] - b[i])
    return residual_vector



def gauss_elimination(A, b):
    n = len(A)
    for i in range(n):
        max_row = i
        for j in range(i + 1, n):
            if abs(A[j][i]) > abs(A[max_row][i]):
                max_row = j
        A[i], A[max_row] = A[max_row], A[i]
        b[i], b[max_row] = b[max_row], b[i]
        for j in range(i + 1, n):
            factor = A[j][i] / A[i][i]
            b[j] -= factor * b[i]
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]

    x = zeros(n)
    for i in range(n - 1, -1, -1):
        s = sum(A[i][j] * x[j] for j in range(i, n))
        x[i] = (b[i] - s) / A[i][i]

    return x


def main():
    A = [[2.30, 5.70, -0.80],
        [3.50, -2.70, 5.30], 
        [1.70, 2.30, -1.80]]
    
    b = [-6.49, 19.20, -5.09]

    ACopy = copy.deepcopy(A)
    bCopy = copy.deepcopy(b)

    x = gauss_elimination(A, b)
    residual_vector = residual_vector_calculation(ACopy, x, bCopy)

    max = 0
    for i in residual_vector:
        if abs(i) > max:
            max = abs(i)

    print("Решение СЛАУ: ", x)
    print("Вектор невязки: ", residual_vector)
    print("Норма вектора невязки:", max)
main()







# # Импортируем библиотеку NumPy для работы с матрицами
# import numpy as np

# # Определяем функцию, которая принимает квадратную матрицу A и возвращает ее LDLT-факторизацию
# def ldlt_factorization(A):
#     # Проверяем, что матрица A квадратная и симметричная
#     n = A.shape[0]
#     assert A.shape == (n, n), "Матрица A должна быть квадратной"
#     assert np.allclose(A, A.T), "Матрица A должна быть симметричной"
    
#     # Инициализируем матрицы L и D
#     L = np.eye(n) # Единичная матрица размера n
#     D = np.zeros((n, n)) # Нулевая матрица размера n
    
#     # Проходим по диагональным элементам матрицы A
#     for i in range(n):
#         # Вычисляем диагональный элемент матрицы D
#         D[i, i] = A[i, i] - np.sum(L[i, :i] ** 2 * D[:i, :i])
        
#         # Вычисляем недиагональные элементы матрицы L
#         for j in range(i + 1, n):
#             L[j, i] = (A[j, i] - np.sum(L[j, :i] * L[i, :i] * D[:i, :i])) / D[i, i]
    
#     # Возвращаем матрицы L и D
#     return L, D

# # Определяем функцию, которая принимает матрицу A и вектор b и решает систему Ax = b методом LDLT-факторизации
# def ldlt_solve(A, b):
#     # Получаем LDLT-факторизацию матрицы A
#     L, D = ldlt_factorization(A)
    
#     # Решаем систему Ly = b методом прямой подстановки
#     y = np.linalg.solve(L, b)
    
#     # Решаем систему Dz = y методом прямой подстановки
#     z = np.linalg.solve(D, y)
    
#     # Решаем систему L^Tx = z методом обратной подстановки
#     x = np.linalg.solve(L.T, z)
    
#     # Возвращаем вектор x
#     return x

# # Пример использования функции ldlt_solve для решения системы Ax = b

# # Задаем матрицу A и вектор b
# A = np.array([[4, -2, 2], [-2, 2, -4], [2, -4, 11]])
# b = np.array([1, -2, 3])

# # Решаем систему Ax = b методом LDLT-факторизации
# x = ldlt_solve(A, b)

# # Выводим результат на экран
# print("Решение системы Ax = b:")
# print(x)
