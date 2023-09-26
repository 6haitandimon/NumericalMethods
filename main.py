import numpy as np

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

    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        s = sum(A[i][j] * x[j] for j in range(i, n))
        x[i] = (b[i] - s) / A[i][i]

    return x

A = [ [2.30, 5.70, 3.50], [-2.70, 1.70, 2.30], [-0.80, 5.30, -1.80] ]
b = [-6.49, 19.20, -5.09]

print(gauss_elimination(A, b))