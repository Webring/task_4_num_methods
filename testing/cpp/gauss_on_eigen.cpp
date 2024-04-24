#include <iostream>
#include "libs/Eigen/Dense"
#include "chrono"

#define EPS 1e-20 // для 4 задачи поставить 1e-15


using namespace std;
using namespace chrono;
using namespace Eigen;

double find_error(VectorXd &x, VectorXd &x_star) {
    return (x - x_star).norm() / x_star.norm();
}

/*
Для 4 работы
double start_matrix_value(int i, int j) {
    if (i == j) return 100;
    return 10 + 0.1 * i - 0.3 * j;
}*/

double start_matrix_value(int i, int j) {
    return 1.0 / (1 + 0.1 * i + 2 * j);
}

int find_main_element(MatrixXd A, int column_index) {
    int row_index = column_index;
    for (int i = column_index + 1; i < A.rows(); i++) {
        if (abs(A(i, column_index)) > abs(A(row_index, column_index))) {
            row_index = i;
        };
    }

    if (abs(A(row_index, column_index)) < EPS) {
        throw runtime_error("Gauss_Method: degenerate matrix...");
    }
    return row_index;
}

void direct_way(MatrixXd &A, VectorXd &f) {
    for (int i = 0; i < A.rows() - 1; i++) {
        int main_element_index = find_main_element(A, i);

        if (main_element_index != i) {
            A.row(i).swap(A.row(main_element_index));

            f.row(i).swap(f.row(main_element_index));
        }

        //для оставшихся строк выполним умножение слева на матрицу преобразований
        for (int j = i + 1; j < A.rows(); j++) {
            double koef = A(j, i) / A(i, i);

            //для уменьшения ошибок вычислений обнуляемые компоненты занулим явно
            A(j, i) = 0;

            //вычитаем элементы строки i из строк от i + 1 до rows
            for (int k = i + 1; k < A.rows(); k++) {
                A(j, k) -= koef * A(i, k);
            }
            f(j) -= koef * f(i);
        }
    }
}

VectorXd back_row_substitution(MatrixXd &A, VectorXd &f) {
    VectorXd result(f);

    for (int i = f.rows() - 1; i >= 0; i--) {
        if (abs(A(i, i)) < EPS) {
            throw runtime_error("Back Row Substitution: A.division by 0... ");
        }

        for (int j = i + 1; j < f.rows(); j++) {
            result(i) -= A(i, j) * result(j);
        }

        result(i) /= A(i, i);
    }

    return result;
}

int main() {
    int N;
    cin >> N;

    MatrixXd A(N, N);
    VectorXd x(N);
    VectorXd x_star = VectorXd::Ones(N);


    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            A(x, y) = start_matrix_value(x, y);
        }
    }

    VectorXd f = A * x_star;

    auto start = steady_clock::now();

    direct_way(A, f);
    x = back_row_substitution(A, f);

    auto end = steady_clock::now();

    cout << duration_cast<chrono::microseconds>(end - start).count() << endl;
    cout << find_error(x, x_star) << endl;
    return 0;
}
