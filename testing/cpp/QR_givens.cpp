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


VectorXd back_row_substitution(MatrixXd &A, VectorXd &f) {
    VectorXd result(f);

    for (int i = f.rows() - 1; i >= 0; i--) {
        if (abs(A(i, i)) < EPS) {
            throw runtime_error("Back Row Substitution: A.division by 0... ");
        }

        //двигаемся по столбцам
        for (int j = i + 1; j < f.rows(); j++) {
            result(i) -= A(i, j) * result(j);
        }

        result(i) /= A(i, i);
    }

    return result;
}

void givens_orthogonalization(const MatrixXd &A, MatrixXd &Q, MatrixXd &R) {
    R = A;

    // Алгоритм вращения Гивенса: для каждого столбца
    for (int j = 0; j < R.cols() - 1; j++) {
        // Проходим по строкам в столбце
        for (int i = j + 1; i < R.rows(); i++) {
            // Если очередной элемент под диагональю не нулевой, то требуется поворот вектора
            if (abs(R(i, j)) > EPS) // Пороговое значение EPS заменено на 1e-10
            {
                double cos_theta, sin_theta;
                double denominator = hypot(R(j, j), R(i, j));
                cos_theta = R(j, j) / denominator;
                sin_theta = R(i, j) / denominator;

                // A_new = Gt * A
                for (int k = j; k < R.cols(); k++) {
                    double R_jk = R(j, k);
                    double R_ik = R(i, k);
                    double R_new_jk = cos_theta * R_jk + sin_theta * R_ik;
                    double R_new_ik = cos_theta * R_ik - sin_theta * R_jk;
                    R(j, k) = R_new_jk;
                    R(i, k) = R_new_ik;
                }

                // Перемножаем строки матрицы Q на транспонированную матрицу преобразования Q = Q * G
                for (int k = 0; k < Q.rows(); k++) {
                    double Q_kj = Q(k, j);
                    double Q_ki = Q(k, i);
                    double Q_new_kj = cos_theta * Q_kj + sin_theta * Q_ki;
                    double Q_new_ki = cos_theta * Q_ki - sin_theta * Q_kj;
                    Q(k, j) = Q_new_kj;
                    Q(k, i) = Q_new_ki;
                }
            }
        }
    }
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

    MatrixXd Q = MatrixXd::Identity(N, N);;
    MatrixXd R(N, N);

    auto start = steady_clock::now();

    givens_orthogonalization(A, Q, R);

    x = Q.transpose() * f;
    x = back_row_substitution(R, x);

    auto end = steady_clock::now();

    cout << duration_cast<chrono::microseconds>(end - start).count() << endl;
    cout << find_error(x, x_star) << endl;
    return 0;
}
