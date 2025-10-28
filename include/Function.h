#ifndef FUNCTION_H
#define FUNCTION_H

#include "Matrix.h"
#include "Inverse.h"
#include <iostream>

template<typename T>
void PerformLUDecomposition(const Matrix<T>& matrix, Matrix<T>& lower, Matrix<T>& upper) {
    if (!matrix.IsSquare()) {
        throw std::invalid_argument("LU decomposition requires square matrix");
    }

    const size_t size = matrix.GetRowCount();
    lower = Matrix<T>(size, size, T(0));
    upper = Matrix<T>(size, size, T(0));

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = i; j < size; ++j) {
            T sum = T(0);
            for (size_t k = 0; k < i; ++k) {
                sum += lower(i, k) * upper(k, j);
            }
            upper(i, j) = matrix(i, j) - sum;
        }

        for (size_t j = i; j < size; ++j) {
            if (i == j) {
                lower(i, i) = T(1);
            } else {
                T sum = T(0);
                for (size_t k = 0; k < i; ++k) {
                    sum += lower(j, k) * upper(k, i);
                }
                if (upper(i, i) == T(0)) {
                    throw std::runtime_error("Matrix is singular");
                }
                lower(j, i) = (matrix(j, i) - sum) / upper(i, i);
            }
        }
    }
}

void DemonstrateMatrixOperations() {
    try {
        std::cout << "=== Matrix Operations Demonstration ===\n\n";
        
        std::cout << "1. Testing with double type:\n";
        Matrix<double> matrix_a(3, 3, 1.0, 5.0);
        std::cout << "Matrix A:\n" << matrix_a << "\n\n";
        
        std::cout << "Trace of A: " << matrix_a.Trace() << "\n\n";
        
        std::cout << "Inverse of A:\n";
        Matrix<double> inverse_a = CalculateInverseMatrix3x3(matrix_a);
        std::cout << inverse_a << "\n\n";
        
        std::cout << "Verification (A * A^(-1)):\n";
        Matrix<double> verification = matrix_a * inverse_a;
        std::cout << verification << "\n\n";
        
        std::cout << "Inverse verification: " 
                  << (VerifyInverseMatrix(matrix_a, inverse_a) ? "SUCCESS" : "FAILED") << "\n\n";
        
        std::cout << "2. Testing with int type:\n";
        Matrix<int> matrix_b(3, 3, 1, 10);
        std::cout << "Matrix B:\n" << matrix_b << "\n\n";
        
        std::cout << "B + B:\n" << (matrix_b + matrix_b) << "\n\n";
        std::cout << "B * 2:\n" << (matrix_b * 2) << "\n\n";
        
        std::cout << "3. Testing with complex<double> type:\n";
        Matrix<std::complex<double>> complex_matrix(3, 3, 
            std::complex<double>(1.0, 0.0), 
            std::complex<double>(3.0, 2.0)
        );
        std::cout << "Complex matrix:\n" << complex_matrix << "\n\n";
        
        std::cout << "Complex matrix trace: " << complex_matrix.Trace() << "\n\n";
        
        std::cout << "4. LU Decomposition:\n";
        Matrix<double> lower, upper;
        PerformLUDecomposition(matrix_a, lower, upper);
        std::cout << "Lower triangular L:\n" << lower << "\n\n";
        std::cout << "Upper triangular U:\n" << upper << "\n\n";
        std::cout << "L * U (should equal A):\n" << (lower * upper) << "\n\n";
        
    } catch (const std::exception& error) {
        std::cout << "Error: " << error.what() << std::endl;
    }
}

#endif