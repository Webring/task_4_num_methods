#!/bin/python3

import datetime
from typing import Tuple
import numpy as np

# Установка порогового значения для сравнения чисел с плавающей запятой
EPS = 1e-15


# Функция для вычисления относительной ошибки
def find_error(x: np.ndarray, x_star: np.ndarray) -> np.float64:
    return np.linalg.norm(x - x_star) / np.linalg.norm(x_star)


# Функция для инициализации значений матрицы A
def start_matrix_value(i: int, j: int) -> np.float64:
    if i == j:
        return np.float64(100.0)
    return np.float64(10.0 + 0.1 * i - 0.3 * j)


# Функция обратной подстановки
def back_row_substitution(A: np.ndarray, f: np.ndarray) -> np.ndarray:
    result = np.copy(f)

    # Проходим по строкам в обратном порядке
    for i in range(f.shape[0] - 1, -1, -1):
        # Проверка на деление на ноль
        if abs(A[i, i]) < EPS:
            raise RuntimeError("Back Row Substitution: A.division by 0... ")

        # Обратная подстановка
        for j in range(i + 1, f.shape[0]):
            result[i] -= A[i, j] * result[j]

        result[i] /= A[i, i]

    return result

"""

# Функция ортогонализации методом вращений Гивенса
def givens_orthogonalization(A: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    # Инициализация матриц Q и R
    Q = np.eye(A.shape[0], dtype=np.float64)
    R = np.copy(A)

    # Проходим по столбцам
    for j in range(R.shape[1] - 1):
        # Проходим по строкам в столбце
        for i in range(j + 1, R.shape[0]):
            # Если элемент под диагональю не равен нулю
            if abs(R[i, j]) > EPS:
                # Вычисление косинуса и синуса угла вращения
                cos_theta = R[j, j] / np.hypot(R[j, j], R[i, j])
                sin_theta = R[i, j] / np.hypot(R[j, j], R[i, j])

                # Применение вращения к матрицам R и Q
                for k in range(j, R.shape[1]):
                    R_jk: np.float64 = R[j, k]
                    R_ik: np.float64 = R[i, k]
                    R[j, k] = cos_theta * R_jk + sin_theta * R_ik
                    R[i, k] = cos_theta * R_ik - sin_theta * R_jk

                for k in range(Q.shape[0]):
                    Q_kj = Q[k, j]
                    Q_ki = Q[k, i]
                    Q[k, j] = cos_theta * Q_kj + sin_theta * Q_ki
                    Q[k, i] = cos_theta * Q_ki - sin_theta * Q_kj

    return Q, R
"""


def givens_rotation(a, b):
    """
    Вычисляет параметры вращения Гивенса для обнуления b.
    """
    if b == 0:
        return 1, 0
    if abs(b) > abs(a):
        tau = -a / b
        s = 1 / np.sqrt(1 + tau ** 2)
        c = s * tau
    else:
        tau = -b / a
        c = 1 / np.sqrt(1 + tau ** 2)
        s = c * tau
    return c, s


def qr_givens(A):
    """
    Выполняет QR-разложение матрицы A с использованием вращений Гивенса.

    Возвращает матрицы Q и R такие, что A = QR.
    """
    m, n = A.shape
    Q = np.eye(m)
    R = np.copy(A)

    for j in range(n):
        for i in range(m - 1, j, -1):
            if R[i, j] != 0:
                c, s = givens_rotation(R[j, j], R[i, j])
                G = np.array([[c, -s], [s, c]])
                R[[j, i], j:] = np.dot(G, R[[j, i], j:])
                Q[:, [j, i]] = np.dot(Q[:, [j, i]], G.T)
    return Q, R

# Основная функция программы
def main() -> None:
    # Ввод размера матрицы от пользователя
    N = int(input())

    # Инициализация матрицы A и вектора x_star
    A = np.zeros((N, N), dtype=np.float64)
    x_star = np.ones(N, dtype=np.float64)

    # Заполнение матрицы A начальными значениями
    for i in range(N):
        for j in range(N):
            A[i, j] = start_matrix_value(i, j)

    # Вычисление вектора f
    f = A.dot(x_star)

    # Засекаем время начала вычислений
    start_time = datetime.datetime.now()

    # Ортогонализация и решение системы уравнений
    Q, R = qr_givens(A)
    y = np.dot(Q.T, f)
    x = back_row_substitution(R, y)

    # Засекаем время окончания вычислений
    end_time = datetime.datetime.now()

    # Вывод результатов
    print((end_time - start_time).microseconds)
    print(find_error(x, x_star))


if __name__ == "__main__":
    main()
