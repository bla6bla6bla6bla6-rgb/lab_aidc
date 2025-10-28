#ifndef INVERSE_H
#define INVERSE_H

#include "Matrix.h"
#include <cmath>

template<typename T>
T CalculateDeterminant3x3(const Matrix<T>& matrix) {
    if (!matrix.IsSquare() || matrix.GetRowCount() != 3) {
        throw std::invalid_argument("Determinant calculation requires 3x3 matrix");
    }

    return matrix(0,0) * (matrix(1,1) * matrix(2,2) - matrix(1,2) * matrix(2,1)) -
           matrix(0,1) * (matrix(1,0) * matrix(2,2) - matrix(1,2) * matrix(2,0)) +
           matrix(0,2) * (matrix(1,0) * matrix(2,1) - matrix(1,1) * matrix(2,0));
}

template<typename T>
Matrix<T> CalculateCofactorMatrix3x3(const Matrix<T>& matrix) {
    if (!matrix.IsSquare() || matrix.GetRowCount() != 3) {
        throw std::invalid_argument("Cofactor matrix calculation requires 3x3 matrix");
    }

    Matrix<T> cofactor(3, 3, T{0});
    
    cofactor(0, 0) = matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1);
    cofactor(0, 1) = -(matrix(1, 0) * matrix(2, 2) - matrix(1, 2) * matrix(2, 0));
    cofactor(0, 2) = matrix(1, 0) * matrix(2, 1) - matrix(1, 1) * matrix(2, 0);
    
    cofactor(1, 0) = -(matrix(0, 1) * matrix(2, 2) - matrix(0, 2) * matrix(2, 1));
    cofactor(1, 1) = matrix(0, 0) * matrix(2, 2) - matrix(0, 2) * matrix(2, 0);
    cofactor(1, 2) = -(matrix(0, 0) * matrix(2, 1) - matrix(0, 1) * matrix(2, 0));
    
    cofactor(2, 0) = matrix(0, 1) * matrix(1, 2) - matrix(0, 2) * matrix(1, 1);
    cofactor(2, 1) = -(matrix(0, 0) * matrix(1, 2) - matrix(0, 2) * matrix(1, 0));
    cofactor(2, 2) = matrix(0, 0) * matrix(1, 1) - matrix(0, 1) * matrix(1, 0);
    
    return cofactor;
}

template<typename T>
Matrix<T> TransposeMatrix(const Matrix<T>& matrix) {
    Matrix<T> result(matrix.GetColCount(), matrix.GetRowCount(), T{0});
    
    for (size_t i = 0; i < matrix.GetRowCount(); ++i) {
        for (size_t j = 0; j < matrix.GetColCount(); ++j) {
            result(j, i) = matrix(i, j);
        }
    }
    
    return result;
}

template<typename T>
Matrix<T> CalculateInverseMatrix3x3(const Matrix<T>& matrix) {
    if (!matrix.IsSquare() || matrix.GetRowCount() != 3) {
        throw std::invalid_argument("Inverse calculation requires 3x3 matrix");
    }

    T determinant = CalculateDeterminant3x3(matrix);
    
    if (std::abs(determinant) < 1e-10) {
        throw std::runtime_error("Matrix is singular and cannot be inverted");
    }
    
    Matrix<T> cofactor = CalculateCofactorMatrix3x3(matrix);
    Matrix<T> adjugate = TransposeMatrix(cofactor);
    
    return adjugate * (T{1} / determinant);
}

template<typename T>
Matrix<T> CreateIdentityMatrix(size_t size) {
    Matrix<T> identity(size, size, T{0});
    for (size_t i = 0; i < size; ++i) {
        identity(i, i) = T{1};
    }
    return identity;
}

template<typename T>
bool VerifyInverseMatrix(const Matrix<T>& original, const Matrix<T>& inverse, double tolerance = 1e-6) {
    if (!original.IsSquare() || !inverse.IsSquare() || original.GetRowCount() != inverse.GetRowCount()) {
        return false;
    }
    
    Matrix<T> product = original * inverse;
    Matrix<T> identity = CreateIdentityMatrix<T>(original.GetRowCount());
    
    for (size_t i = 0; i < product.GetRowCount(); ++i) {
        for (size_t j = 0; j < product.GetColCount(); ++j) {
            if (std::abs(product(i, j) - identity(i, j)) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

#endif