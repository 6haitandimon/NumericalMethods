from math import sqrt
import numpy as np


a = 1
k = 2
t0, T = [0, 1]
u10, u20, u30 = [1, 1, 1]
hmax = 1
u = [1, 1, 1]
h = 1e-1

def f(u2, u3):
    return (k - a / a) * u2 * u3

def g(u1, u3):
    return (a + k / k) * u1 * u3

def h_func(u1, u2):
    return (a - k / a) * u1 * u2

def three_zone_rule(e, e_additional, h):
    if abs(e) < e_additional:
        return h / 2
    if e_additional / 4 < abs(e) and  abs(e) <= abs(e_additional):
        return h
    if abs(e) > abs(e_additional) / 4:
        return h * 2

def explicit_euler(e=1e-3):
    '''Явный метод Эйлера'''
    print("Явный метод Эйлер")
    t, h = t0, hmax
    u = np.array([u10, u20, u30], dtype=float)
    i = 0
    while t < T:
        i += 1
        uu = np.array([f(u[1], u[2]), g(u[0], u[2]), h_func(u[0], u[1])], dtype=float)
        h = e / (np.linalg.norm(uu) + e / hmax) #стр. 20 (методичка)
        u += h * uu
        t += h
        print(f"{i} {u[0]} {u[1]} {u[2]} {t}")


# def implicit_euler(u = np.array([u10, u20, u30]), 
#                    func = lambda uu: uu - u - h * np.array([f(uu[1], uu[2]), g(uu[0], uu[2]), h_func(uu[0], uu[1])]),
#                    e=1e-1):
    
#     '''Неявный метод Эйлера'''
#     print("\nНеявный метод Эйлера\n")
#     t, h = t0, hmax
#     i = 0

#     while t < T:
#         i += 1
#         e_additional = max(u)
#         h = three_zone_rule(e, e_additional, h)
#         u = newton(func, u)
#         t += h
#         print(f"{i} {u[0]} {u[1]} {u[2]} {t}")


def main(): 

    explicit_euler()
    # implicit_euler()

    # 1 = 1, 2 = 2, 3 = 2

    # A = 1/6 * np.array([[10, 2, 2], 
    #                     [2, 10, 2], 
    #                     [2, -2, 10]])

    # b = -1/6 * np.array([-12, -16, 10])

    # u = np.array([10, 22, 9])

    # def f1(uu):
    #     return uu - u - h * (np.dot(A, uu) - b)

    # implicit_euler(func = f1, u = u)

main()