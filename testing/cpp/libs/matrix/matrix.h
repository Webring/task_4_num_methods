#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>

using namespace std;

class Matrixd {
private:
    int rows;
    int cols;
    vector<vector<double>> data;

public:
    Matrixd();

    Matrixd(int rows, int cols);

    Matrixd(const Matrixd &other);

    int n_rows();

    int n_cols();

    Matrixd operator+(const Matrixd &other) const;

    Matrixd operator-(const Matrixd &other) const;

    Matrixd operator*(const Matrixd &other) const;

    double &operator()(int row, int col);

    Matrixd trn() const;

    vector<double> &get_row(int row);

    vector<double> get_column(int col) const;

    void swap_row(int row1, int row2);

    void swap_column(int col1, int col2);

    void print() const;

    void to_zeros();

    void to_ones();

    void fill(double value);

    double norm();
};


#endif // MATRIX_H
