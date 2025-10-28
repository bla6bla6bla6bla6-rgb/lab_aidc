#ifndef EQUATION_SOLVER_H
#define EQUATION_SOLVER_H

#include "Matrix.h"
#include "Inverse.h"
#include <stdexcept>

template<typename T>
Matrix<T> SolveLinearEquation(const Matrix<T>& A, const Matrix<T>& b) {
    if (!A.IsSquare() || A.GetRowCount() != 3) {
        throw std::invalid_argument("Matrix A must be 3x3");
    }
    if (b.GetRowCount() != 3 || b.GetColCount() != 1) {
        throw std::invalid_argument("Vector b must be 3x1");
    }
    
    // solve: x = A^(-1) * b
    Matrix<T> A_inv = CalculateInverseMatrix3x3(A);
    return A_inv * b;
}

#endif