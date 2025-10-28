#include <iostream>
#include "Matrix.h"

int main() {
    // Test random double matrix
    Matrix<double> random_matrix(2, 3, 1.0, 10.0);
    std::cout << "Random double matrix (1.0 to 10.0):\n" << random_matrix << "\n";
    
    // Test random int matrix  
    Matrix<int> random_int_matrix(3, 2, 5, 15);
    std::cout << "Random int matrix (5 to 15):\n" << random_int_matrix << "\n";
    
    return 0;
}
