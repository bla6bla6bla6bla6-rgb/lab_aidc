#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <random>
#include <stdexcept>
#include <complex>
#include <type_traits>
#include <cmath>

template <typename T>
class Matrix {
private:
    T* data_ptr_;
    size_t rows_;
    size_t cols_;
    inline static const double kEpsilon = 0.001;

    void InitializeWithValue(T value) {
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            data_ptr_[i] = value;
        }
    }

    void Cleanup() {
        if (data_ptr_) {
            delete[] data_ptr_;
            data_ptr_ = nullptr;
        }
    }

public:
    Matrix() : data_ptr_(nullptr), rows_(0), cols_(0) {}

    // Constructor with fixed value
    Matrix(size_t rows, size_t cols, T value) 
        : data_ptr_(new T[rows * cols]), rows_(rows), cols_(cols) {
        InitializeWithValue(value);
    }

    // Constructor with random values 
    Matrix(size_t rows, size_t cols, T min_val, T max_val) 
        : data_ptr_(new T[rows * cols]), rows_(rows), cols_(cols) {
        
        std::random_device rd;
        std::mt19937 engine(rd());
        
        if constexpr (std::is_integral_v<T>) {
            std::uniform_int_distribution<T> distribution(
                std::min(min_val, max_val), 
                std::max(min_val, max_val)   
            );
            for (size_t i = 0; i < rows_ * cols_; ++i) {
                data_ptr_[i] = distribution(engine);
            }
        } else if constexpr (std::is_floating_point_v<T>) {
            std::uniform_real_distribution<T> distribution(   
                std::min(min_val, max_val), 
                std::max(min_val, max_val)
            );
            for (size_t i = 0; i < rows_ * cols_; ++i) {
                data_ptr_[i] = distribution(engine);
            }
        } else {
            // For other types, use the value initialization
            InitializeWithValue(T{});
        }
    }

    // Constructor for complex numbers
    template<typename U>
    Matrix(size_t rows, size_t cols, std::complex<U> min_val, std::complex<U> max_val)
        : data_ptr_(new T[rows * cols]), rows_(rows), cols_(cols) {
        
        std::random_device rd;
        std::mt19937 engine(rd());

        std::uniform_real_distribution<U> dist_real(
            std::min(min_val.real(), max_val.real()), 
            std::max(min_val.real(), max_val.real())
        );
        std::uniform_real_distribution<U> dist_imag(
            std::min(min_val.imag(), max_val.imag()), 
            std::max(min_val.imag(), max_val.imag())
        );
        
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            data_ptr_[i] = std::complex<U>(dist_real(engine), dist_imag(engine));
        }
    }

    // Copy constructor
    Matrix(const Matrix<T>& other) 
        : data_ptr_(other.data_ptr_ ? new T[other.rows_ * other.cols_] : nullptr),
          rows_(other.rows_), cols_(other.cols_) {
        if (data_ptr_) {
            for (size_t i = 0; i < rows_ * cols_; ++i) {
                data_ptr_[i] = other.data_ptr_[i];
            }
        }
    }

    ~Matrix() {
        Cleanup();
    }

    // Accessors
    size_t GetRowCount() const { return rows_; }
    size_t GetColCount() const { return cols_; }
    bool IsEmpty() const { return rows_ == 0 || cols_ == 0; }
    bool IsSquare() const { return rows_ == cols_; }

    T operator()(size_t row, size_t col) const {
        if (row >= rows_ || col >= cols_) {
            throw std::out_of_range("Matrix indices out of range");
        }
        return data_ptr_[row * cols_ + col];
    }

    T& operator()(size_t row, size_t col) {
        if (row >= rows_ || col >= cols_) {
            throw std::out_of_range("Matrix indices out of range");
        }
        return data_ptr_[row * cols_ + col];
    }

    // Arithmetic operations
    Matrix<T> operator+(const Matrix<T>& rhs) const {
        if (rows_ != rhs.rows_ || cols_ != rhs.cols_) {
            throw std::invalid_argument("Matrices must have same dimensions");
        }
        Matrix<T> result(rows_, cols_, T{0});
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            result.data_ptr_[i] = data_ptr_[i] + rhs.data_ptr_[i];
        }
        return result;
    }

    // Operator subtraction 
    Matrix<T> operator-(const Matrix<T>& rhs) const {
        if (rows_ != rhs.rows_ || cols_ != rhs.cols_) {
            throw std::invalid_argument("Matrices must have same dimensions");
        }
        Matrix<T> result(rows_, cols_, T{0});
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            result.data_ptr_[i] = data_ptr_[i] - rhs.data_ptr_[i];
        }
        return result;
    }

    Matrix<T> operator*(T scalar) const {
        Matrix<T> result(rows_, cols_, T{0});
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            result.data_ptr_[i] = data_ptr_[i] * scalar;
        }
        return result;
    }

    // Operator division by scalar 
    Matrix<T> operator/(T scalar) const {
        if (scalar == T{0}) {
            throw std::invalid_argument("Division by zero");
        }
        Matrix<T> result(rows_, cols_, T{0});
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            result.data_ptr_[i] = data_ptr_[i] / scalar;
        }
        return result;
    }

    Matrix<T> operator*(const Matrix<T>& rhs) const {
        if (cols_ != rhs.rows_) {
            throw std::invalid_argument("Matrix dimensions don't match for multiplication");
        }
        
        Matrix<T> result(rows_, rhs.cols_, T{0});
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < rhs.cols_; ++j) {
                T sum = T{0};
                for (size_t k = 0; k < cols_; ++k) {
                    sum += (*this)(i, k) * rhs(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    // Comparison operators 
    bool operator==(const Matrix<T>& rhs) const {
        if (rows_ != rhs.rows_ || cols_ != rhs.cols_) return false;
        
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            if constexpr (std::is_floating_point_v<T>) {
                if (std::abs(data_ptr_[i] - rhs.data_ptr_[i]) > kEpsilon) return false;
            } else {
                if (data_ptr_[i] != rhs.data_ptr_[i]) return false;
            }
        }
        return true;
    }

    bool operator!=(const Matrix<T>& rhs) const {
        return !(*this == rhs);
    }

    T Trace() const {
        if (!IsSquare()) {
            throw std::invalid_argument("Trace is defined only for square matrices");
        }
        T trace = T{0};
        for (size_t i = 0; i < rows_; ++i) {
            trace += (*this)(i, i);
        }
        return trace;
    }
    };

    template <typename T>
    Matrix<T> operator*(T scalar, const Matrix<T>& matrix) {
        return matrix * scalar;
    }   

    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
        for (size_t i = 0; i < matrix.GetRowCount(); ++i) {
            for (size_t j = 0; j < matrix.GetColCount(); ++j) {
                os << matrix(i, j);
                if (j < matrix.GetColCount() - 1) {
                    os << " ";
                }
            }
            if (i < matrix.GetRowCount() - 1) {
                os << "\n";
            }
        }
        return os;
    }
    #endif