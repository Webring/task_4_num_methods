#include "matrix.h"
#include "cmath"

Matrixd::Matrixd() : rows(0), cols(0) {}

Matrixd::Matrixd(int rows, int cols) : rows(rows), cols(cols),
                                       data(vector<vector<double>>(rows, vector<double>(cols))) {}

Matrixd::Matrixd(const Matrixd &other) : rows(other.rows), cols(other.cols),
                                         data(other.data.size(), std::vector<double>(other.cols)) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            data[i][j] = other.data[i][j];
        }
    }
}

int Matrixd::n_rows() {
    return rows;
}

int Matrixd::n_cols() {
    return cols;
}


Matrixd Matrixd::operator+(const Matrixd &other) const {
    if (rows != other.rows || cols != other.cols) {
        throw invalid_argument("Matrix dimensions must be equal for addition.");
    }

    Matrixd result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] + other.data[i][j];
        }
    }
    return result;
}

Matrixd Matrixd::operator-(const Matrixd &other) const {
    if (rows != other.rows || cols != other.cols) {
        throw invalid_argument("Matrix dimensions must be equal for subtraction.");
    }

    Matrixd result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] - other.data[i][j];
        }
    }
    return result;
}

Matrixd Matrixd::operator*(const Matrixd &other) const {
    if (cols != other.rows) {
        throw invalid_argument(
                "Number of columns in the first matrix must be equal to the number of rows in the second matrix for multiplication.");
    }

    Matrixd result(rows, other.cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < other.cols; ++j) {
            for (int k = 0; k < cols; ++k) {
                result.data[i][j] += data[i][k] * other.data[k][j];
            }
        }
    }
    return result;
}

double &Matrixd::operator()(int row, int col) {
    if (row >= rows or col >= cols or col < 0 or row < 0) {
        throw runtime_error("Index error");
    }
    return data[row][col];
}

Matrixd Matrixd::trn() const {
    Matrixd result(cols, rows);
    for (int i = 0; i < cols; ++i) {
        for (int j = 0; j < rows; ++j) {
            result.data[i][j] = data[j][i];
        }
    }
    return result;
}

vector<double> &Matrixd::get_row(int row) {
    return data[row];
}

vector<double> Matrixd::get_column(int col) const {
    vector<double> column(rows);
    for (int i = 0; i < rows; ++i) {
        column[i] = data[i][col];
    }
    return column;
}

void Matrixd::swap_row(int row1, int row2) {
    if (row1 >= rows || row2 >= rows) {
        throw out_of_range("Row index out of range.");
    }
    swap(data[row1], data[row2]);
}

void Matrixd::swap_column(int col1, int col2) {
    if (col1 >= cols || col2 >= cols) {
        throw out_of_range("Column index out of range.");
    }
    for (int i = 0; i < rows; ++i) {
        swap(data[i][col1], data[i][col2]);
    }
}

void Matrixd::print() const {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
}

void Matrixd::to_zeros() {
    fill(0.0);
}

void Matrixd::to_ones() {
    fill(1.0);
}

void Matrixd::fill(double value) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            operator()(i, j) = value;
        }
    }
}

double Matrixd::norm() {
    if (rows != 1 and cols != 1) {
        throw runtime_error("Is not a vector");
    }

    double sum = 0.0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            sum += data[i][j] * data[i][j];
        }
    }
    return sqrt(sum);
}