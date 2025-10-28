#include <iostream>
#include "Matrix.h"
#include "Inverse.h"
#include "EquationSolver.h"

void DemonstrateAllFeatures() {
    try {
        std::cout << "=== Complete Matrix Class Demonstration ===\n\n";
        
        std::cout << "1. Basic operations:\n";
        Matrix<double> A(3, 3, 1.0, 5.0);
        Matrix<double> B(3, 3, 2.0, 6.0);
        std::cout << "Matrix A:\n" << A << "\n\n";
        std::cout << "Matrix B:\n" << B << "\n\n";
        
        std::cout << "A + B:\n" << (A + B) << "\n\n";
        std::cout << "A - B:\n" << (A - B) << "\n\n";
        std::cout << "A * B:\n" << (A * B) << "\n\n";
        std::cout << "A * 2.5:\n" << (A * 2.5) << "\n\n";
        std::cout << "A / 2.0:\n" << (A / 2.0) << "\n\n";
        std::cout << "2.5 * A:\n" << (2.5 * A) << "\n\n";
        

        std::cout << "A == B: " << (A == B ? "true" : "false") << "\n";
        std::cout << "A != B: " << (A != B ? "true" : "false") << "\n\n";
        

        std::cout << "Trace of A: " << A.Trace() << "\n\n";
        Matrix<double> A_inv = CalculateInverseMatrix3x3(A);
        std::cout << "Inverse of A:\n" << A_inv << "\n\n";
   
        
        std::cout << "4. Solving linear equation A x = b:\n";
        Matrix<double> b(3, 1, 1.0, 3.0); 
        std::cout << "Vector b:\n" << b << "\n\n";
        
        Matrix<double> x = SolveLinearEquation(A, b);
        std::cout << "Solution x:\n" << x << "\n\n";
        
        
        Matrix<double> verification = A * x;
        std::cout << "Verification (A * x):\n" << verification << "\n";
        std::cout << "Verification successful: " << (verification == b ? "YES" : "NO") << "\n\n";
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
}

int main() {
    DemonstrateAllFeatures();
    return 0;
}
