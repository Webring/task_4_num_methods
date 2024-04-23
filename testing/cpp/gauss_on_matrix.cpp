#include <iostream>
#include "libs/matrix/matrix.h"
#include "chrono"

#define EPS 1e-15

using namespace std;
using namespace chrono;

double find_error(Matrixd &x, Matrixd &x_star) {
    return (x - x_star).norm() / x_star.norm();
}

double start_matrix_value(int i, int j) {
    if (i == j) return 100;
    return 10 + 0.1 * i - 0.3 * j;
}


int find_main_element(Matrixd A, int column_index) {
    int row_index = column_index;
    for (int i = column_index + 1; i < A.n_rows(); i++) {
        if (abs(A(i, column_index)) > abs(A(row_index, column_index))) {
            row_index = i;
        };
    }

    if (abs(A(row_index, column_index)) < EPS) {
        throw runtime_error("Gauss_Method: degenerate matrix...");
    }
    return row_index;
}

void direct_way(Matrixd &A, Matrixd &f) {
    for (int i = 0; i < A.n_rows() - 1; i++) {
        int main_element_index = find_main_element(A, i);

        if (main_element_index != i) {
            A.swap_row(i, main_element_index);
            f.swap_row(i, main_element_index);
        }

        //для оставшихся строк выполним умножение слева на матрицу преобразований
        for (int j = i + 1; j < A.n_rows(); j++) {
            double koef = A(j, i) / A(i, i);

            //для уменьшения ошибок вычислений обнуляемые компоненты занулим явно
            A(j, i) = 0;

            //вычитаем элементы строки i из строк от i + 1 до rows
            for (int k = i + 1; k < A.n_rows(); k++) {
                A(j, k) -= koef * A(i, k);
            }
            f(j, 0) -= koef * f(i, 0);
        }
    }
}

Matrixd back_row_substitution(Matrixd &A, Matrixd &f) {
    Matrixd result(f);


    for (int i = f.n_rows() - 1; i >= 0; i--) {
        if (abs(A(i, i)) < EPS) {
            throw runtime_error("Back Row Substitution: A.division by 0... ");
        }

        for (int j = i + 1; j < f.n_rows(); j++) {
            result(i, 0) -= A(i, j) * result(j, 0);
        }

        result(i, 0) /= A(i, i);
    }

    return result;
}

int main() {
    int N;

    cin >> N;

    Matrixd A(N, N);
    Matrixd x(N, 1);
    Matrixd x_star(N, 1);
    x_star.to_ones();


    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            A(x, y) = start_matrix_value(x, y);
        }
    }


    Matrixd f = A * x_star;

    auto start = steady_clock::now();

    direct_way(A, f);

    x = back_row_substitution(A, f);

    auto end = steady_clock::now();

    cout << duration_cast<chrono::microseconds>(end - start).count() << endl;
    cout << find_error(x, x_star);
    return 0;
}
