#!/bin/python3

import numpy as np
from datetime import datetime

EPS = 1e-15


def find_error(x: np.ndarray, x_star: np.ndarray):
    return np.linalg.norm(x - x_star) / np.linalg.norm(x_star)


def start_matrix_value(i: int, j: int):
    if i == j:
        return 100
    return 10 + 0.1 * i - 0.3 * j


def find_main_element(A: np.ndarray, column_index: int):
    row_index = column_index
    for i in range(column_index + 1, A.shape[0]):
        if abs(A[i, column_index]) > abs(A[row_index, column_index]):
            row_index = i

    if abs(A[row_index, column_index]) < EPS:
        raise RuntimeError("Gauss_Method: degenerate matrix...")
    return row_index


def direct_way(A: np.ndarray, f: np.ndarray):
    for i in range(A.shape[0] - 1):
        main_element_index = find_main_element(A, i)

        if main_element_index != i:
            A[[i, main_element_index]] = A[[main_element_index, i]]
            f[[i, main_element_index]] = f[[main_element_index, i]]

        for j in range(i + 1, A.shape[0]):
            koef = A[j, i] / A[i, i]
            A[j, i] = 0
            A[j, i + 1:] -= koef * A[i, i + 1:]
            f[j] -= koef * f[i]


def back_row_substitution(A: np.ndarray, f: np.ndarray):
    result = np.copy(f)

    for i in range(f.shape[0] - 1, -1, -1):
        if abs(A[i, i]) < EPS:
            raise RuntimeError("Back Row Substitution: A.division by 0... ")

        for j in range(i + 1, f.shape[0]):
            result[i] -= A[i, j] * result[j]

        result[i] /= A[i, i]

    return result


def main() -> None:
    N: int = int(input())

    A = np.zeros((N, N))
    x = np.zeros(N)
    x_star = np.ones(N)

    for i in range(N):
        for j in range(N):
            A[i, j] = start_matrix_value(i, j)

    f = np.dot(A, x_star)

    start = datetime.now()
    direct_way(A, f)
    x = back_row_substitution(A, f)
    end = datetime.now()

    print((end - start).microseconds)
    print(find_error(x, x_star))


if __name__ == '__main__':
    main()
