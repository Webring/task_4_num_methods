#include <iostream>
#include "libs/Eigen/Dense"
#include "chrono"

#define EPS 1e-20


using namespace std;
using namespace chrono;
using namespace Eigen;

double find_error(VectorXd &x, VectorXd &x_star) {
    return (x - x_star).norm() / x_star.norm();
}

double start_matrix_value(int i, int j) {
    return 1.0 / (1 + 0.1 * i + 2 * j);
}

// Функция для проверки сингулярных значений на положительность и коррекции соответствующих матриц
void abs_singular_values(MatrixXd &Sigma, MatrixXd &U) {
    // Определяем минимальный размер матрицы Sigma
    int min_size = min(Sigma.rows(), Sigma.cols());

    // Проверяем сингулярные значения на положительность
    for (int i = 0; i < min_size; i++) {
        if (Sigma(i, i) < 0) {
            // Если значение отрицательное, делаем его положительным и корректируем соответствующий столбец матрицы U
            Sigma(i, i) = -Sigma(i, i);
            for (int j = 0; j < U.rows(); j++) {
                U(j, i) = -U(j, i);
            }
        }
    }
}

// Функция для сортировки сингулярных значений в порядке убывания и соответствующих коррекций матриц U и V
void sort_singular_values(MatrixXd &Sigma, MatrixXd &U, MatrixXd &V) {
    // Определяем минимальный размер матрицы Sigma
    int min_size = min(Sigma.rows(), Sigma.cols());

    // Сортировка сингулярных значений
    for (int I = 0; I < min_size; I++) {
        double max_elem = Sigma(I, I); // Наибольший элемент
        int index = I; // Индекс наибольшего элемента

        // Поиск наибольшего элемента
        for (int i = I + 1; i < min_size; i++) {
            if (Sigma(i, i) > max_elem) {
                max_elem = Sigma(i, i);
                index = i;
            }
        }

        // Если найден наибольший элемент, меняем местами текущий элемент с наибольшим
        if (I != index) {
            Sigma(index, index) = Sigma(I, I);
            Sigma(I, I) = max_elem;

            // Коррекция соответствующих столбцов матриц U и V
            U.col(index).swap(U.col(I));
            V.col(index).swap(V.col(I));
        }
    }
}

void column_transformation(MatrixXd &A, MatrixXd &U, int i, int j) {
    // Вектор отражения
    VectorXd p(A.rows());

    // Вспомогательные переменные
    double s, beta, mu;

    // Находим квадрат нормы столбца для обнуления
    s = 0;
    for (int I = j; I < A.rows(); I++)
        s += pow(A(I, i), 2);

    // Если ненулевые элементы под диагональю есть:
    // квадрат нормы вектора для обнуления не совпадает с квадратом зануляемого элемента
    if (sqrt(abs(s - A(j, i) * A(j, i))) > EPS) {
        // Выбор знака слагаемого beta = sign(-x1)
        if (A(j, i) < 0)
            beta = sqrt(s);
        else
            beta = -sqrt(s);

        // Вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
        mu = 1.0 / beta / (beta - A(j, i));

        // Формируем вектор p
        p.setZero();
        for (int I = 0; I < A.rows(); I++) {
            if (I >= j)
                p(I) = A(I, i);
        }

        // Изменяем элемент, с которого начнётся обнуление
        p(j) -= beta;

        // Вычисляем новые компоненты матрицы A = ... * U2 * U1 * A
        for (int m = 0; m < A.cols(); m++) {
            // Произведение S = St * p
            s = 0;
            for (int n = j; n < A.rows(); n++)
                s += A(n, m) * p(n);
            s *= mu;

            // S = S - 2 * p * (St * p)^t / ||p||^2
            for (int n = j; n < A.rows(); n++)
                A(n, m) -= s * p(n);
        }

        // Вычисляем новые компоненты матрицы U = ... * H2 * H1 * U
        for (int m = 0; m < A.rows(); m++) {
            // Произведение S = Ut * p
            s = 0;
            for (int n = j; n < A.rows(); n++)
                s += U(m, n) * p(n);
            s *= mu;

            // U = U - 2 * p * (Ut * p)^t / ||p||^2
            for (int n = j; n < A.rows(); n++)
                U(m, n) -= s * p(n);
        }
    }
}

void row_transformation(MatrixXd &A, MatrixXd &V, int i, int j) {
    // Вектор отражения
    VectorXd p(A.cols());

    // Вспомогательные переменные
    double s, beta, mu;

    // Находим квадрат нормы строки для обнуления
    s = 0;
    for (int I = j; I < A.cols(); I++)
        s += pow(A(i, I), 2);

    // Если ненулевые элементы под диагональю есть:
    // Квадрат нормы вектора для обнуления не совпадает с квадратом зануляемого элемента
    if (sqrt(abs(s - A(i, j) * A(i, j))) > EPS) {
        // Выбор знака слагаемого beta = sign(-x1)
        if (A(i, j) < 0)
            beta = sqrt(s);
        else
            beta = -sqrt(s);

        // Вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
        mu = 1.0 / beta / (beta - A(i, j));

        // Формируем вектор p
        p.setZero();
        for (int I = 0; I < A.cols(); I++) {
            if (I >= j)
                p(I) = A(i, I);
        }

        // Изменяем диагональный элемент
        p(j) -= beta;

        // Вычисляем новые компоненты матрицы A = A * H1 * H2 ...
        for (int m = 0; m < A.rows(); m++) {
            // Произведение A * p
            s = 0;
            for (int n = j; n < A.cols(); n++)
                s += A(m, n) * p(n);
            s *= mu;

            // A = A - p * (A * p)^t
            for (int n = j; n < A.cols(); n++)
                A(m, n) -= s * p(n);
        }

        // Вычисляем новые компоненты матрицы V = V * H1 * H2 * ...
        for (int m = 0; m < A.cols(); m++) {
            // Произведение V * p
            s = 0;
            for (int n = j; n < A.cols(); n++)
                s += V(m, n) * p(n);
            s *= mu;

            // V = V - p * (V * p)^t
            for (int n = j; n < A.cols(); n++)
                V(m, n) -= s * p(n);
        }
    }
}

void delete_elem_down_triangle(MatrixXd &A, MatrixXd &U, int I, int J) {
    double help1, help2;
    double c = 0, s = 0; // Косинус и синус для поворота

    // Если элемент не нулевой, требуется выполнить поворот вектора
    if (abs(A(I, J)) > EPS) {
        // Вычисляем косинус и синус угла поворота
        help1 = sqrt(pow(A(I, J), 2) + pow(A(J, J), 2));
        c = A(J, J) / help1;
        s = A(I, J) / help1;

        // Применяем преобразование к матрице A: A_new = Gt * A
        for (int k = 0; k < A.cols(); k++) {
            help1 = c * A(J, k) + s * A(I, k);
            help2 = c * A(I, k) - s * A(J, k);
            A(J, k) = help1;
            A(I, k) = help2;
        }

        // Применяем преобразование к матрице U: U = U * G
        for (int k = 0; k < U.rows(); k++) {
            help1 = c * U(k, J) + s * U(k, I);
            help2 = c * U(k, I) - s * U(k, J);
            U(k, J) = help1;
            U(k, I) = help2;
        }
    }
    // Устанавливаем элемент A[I][J] в нулевое значение
    A(I, J) = 0;
}

void delete_elem_up_triangle(MatrixXd &A, MatrixXd &V, int I, int J) {
    double help1, help2;
    double c = 0, s = 0; // Косинус и синус для поворота

    // Если элемент не нулевой, требуется выполнить поворот вектора
    if (abs(A(I, J)) > EPS) {
        // Вычисляем косинус и синус угла поворота
        help1 = sqrt(pow(A(I, J), 2) + pow(A(I, I), 2));
        c = A(I, I) / help1;
        s = -A(I, J) / help1;

        // Применяем преобразование к матрице A: A_new = A * Gt
        for (int k = 0; k < A.rows(); k++) {
            help1 = c * A(k, I) - s * A(k, J);
            help2 = c * A(k, J) + s * A(k, I);
            A(k, I) = help1;
            A(k, J) = help2;
        }

        // Применяем преобразование к матрице V: V = V * Gt
        for (int k = 0; k < V.rows(); k++) {
            help1 = c * V(k, I) - s * V(k, J);
            help2 = c * V(k, J) + s * V(k, I);
            V(k, I) = help1;
            V(k, J) = help2;
        }
    }
}

void start_svd(MatrixXd &A, MatrixXd &U, MatrixXd &Sigma, MatrixXd &V) {
    int min_size = min(A.rows(), A.cols()); // Определяем минимальный размер матрицы
    int up_size = min_size - 1; // Размер верхней диагонали
    int down_size = min_size - 1; // Размер нижней диагонали

    // Этап I: бидиагонализация
    for (int i = 0; i < min_size - 1; i++) {
        column_transformation(Sigma, U, i, i);
        row_transformation(Sigma, V, i, i + 1);
    }

    // Дополнительные преобразования, если размерность матрицы A больше/меньше по строкам или столбцам
    if (A.rows() > A.cols()) {
        column_transformation(Sigma, U, A.cols() - 1, A.cols() - 1);
        down_size += 1;
    }

    if (A.rows() < A.cols()) {
        row_transformation(Sigma, V, A.rows() - 1, A.rows());
        up_size += 1;
    }

    int count_up_elements;
    VectorXd up(up_size);
    VectorXd down(down_size);

    // Этап II: преследование и приведение к диагональному виду
    do {
        count_up_elements = 0;

        // Процедура преследования для верхней диагонали
        for (int i = 0; i < up_size; i++) {
            if (abs(up[i] - Sigma(i, i + 1)) > EPS) {
                up[i] = Sigma(i, i + 1);
                delete_elem_up_triangle(Sigma, V, i, i + 1);
            } else {
                count_up_elements++;
            }
        }

        // Процедура преследования для нижней диагонали
        for (int i = 0; i < down_size; i++) {
            if (abs(down[i] - Sigma(i + 1, i)) > EPS) {
                down[i] = Sigma(i + 1, i);
                delete_elem_down_triangle(Sigma, U, i + 1, i);
            }
        }
    } while (count_up_elements != up_size);

    // Убираем отрицательные сингулярные значения и сортируем их по возрастанию
    abs_singular_values(Sigma, U);
    sort_singular_values(Sigma, U, V);
}

// Функция для уменьшения размерности матриц Sigma, U и V на основе заданного порога Reduction
void reduction_SVD(MatrixXd &U, MatrixXd &Sigma, MatrixXd &V, double Reduction) {
    // Определяем минимальный размер матрицы Sigma
    int min_size = min(Sigma.rows(), Sigma.cols());

    // Проверка на возможность редукции по сингулярным числам
    for (int i = 0; i < min_size; i++) {
        if (abs(Sigma(i, i)) < Reduction) {
            min_size = i;
            break;
        }
    }

    // Редукция размерности матриц Sigma, U и V
    Sigma.conservativeResize(min_size, min_size);
    U.conservativeResize(U.rows(), min_size);
    V.conservativeResize(V.rows(), min_size);
}

const int mrank(MatrixXd &Sigma) {
    return Sigma.rows();
}

// Функция для вычисления модуля определителя матрицы Sigma
double abs_det(MatrixXd &Sigma) {
    int size = mrank(Sigma); // Получаем ранг матрицы Sigma

    if (size == 0) {
        throw runtime_error("Error in SVD.Rank: SVD is not built ...");
    }

    double det = 1;
    for (int i = 0; i < size; i++) {
        det *= Sigma(i, i); // Вычисляем произведение диагональных элементов
    }
    return det;
}

// Функция для вычисления числа обусловленности матрицы Sigma
double cond(MatrixXd &Sigma) {
    int size = mrank(Sigma); // Получаем ранг матрицы Sigma

    if (size == 0) {
        throw runtime_error("Error in SVD.Rank: SVD is not built ...");
    }

    // Число обусловленности - отношение максимального сингулярного значения к минимальному
    return Sigma(0, 0) / Sigma(size - 1, size - 1);
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


    // cout << A << endl << endl;

    // Создаем матрицы U, Sigma и V, инициализируем Sigma значением матрицы A
    MatrixXd U = MatrixXd::Identity(A.rows(), A.rows());
    MatrixXd Sigma = MatrixXd(A);
    MatrixXd V = MatrixXd::Identity(A.cols(), A.cols());


    auto start = steady_clock::now();

    start_svd(A, U, Sigma, V);
    reduction_SVD(U, Sigma, V, 10e-15);
    x = V * Sigma.inverse() * U.transpose() * f; // Псевдорешение x=V * Sigma^-1 * Ut

    auto end = steady_clock::now();

    /*
    cout << U << endl << endl;
    cout << Sigma << endl << endl;
    cout << V << endl << endl;
    cout << U * Sigma * V.transpose() << endl << endl;
     */

    cout << duration_cast<chrono::microseconds>(end - start).count() << endl;
    cout << "rank A = " << mrank(Sigma) << endl;
    cout << "Cond A = " << cond(Sigma) << endl;
    cout << "det A = " << abs_det(Sigma) << endl;
    cout << find_error(x, x_star) << endl;
    return 0;
}
