#ifndef GENTRYCRYPTO_MATH_LINALG_HPP
#define GENTRYCRYPTO_MATH_LINALG_HPP

#include <iostream>

struct Col
{
    size_t value = 1;
    explicit Col() {}
    explicit Col(size_t value): value(value) {}
};

struct Row
{
    size_t value = 1;
    explicit Row() {}
    explicit Row(size_t value): value(value) {}
};

template <typename T>
class Matrix
{
    public:
        explicit Matrix();
        Matrix(Row, Col);
        Matrix(Row, Col, T**);
        ~Matrix();
        
        Matrix(const Matrix&);
        Matrix& operator=(const Matrix&);
        
        Matrix(Matrix&&);
        Matrix& operator=(Matrix&&);     

        Matrix operator+=(const Matrix&);
        Matrix operator-=(const Matrix&);
        Matrix operator*=(const Matrix&);
        Matrix operator/=(const Matrix&) = delete;

        Matrix operator*=(const double&);
        Matrix operator/=(const double&);

        Matrix operator+() const;
        Matrix operator-() const;

        Matrix& operator++() = delete;
        Matrix& operator--() = delete;

        Matrix& operator++(int) = delete;
        Matrix& operator--(int) = delete;

        bool operator< (const Matrix&) const = delete;
        bool operator> (const Matrix&) const = delete;
        bool operator<=(const Matrix&) const = delete;
        bool operator>=(const Matrix&) const = delete;

        T& operator()(size_t, size_t);

        Matrix Transpose();
        //Matrix Determinant();

    private:
        void Memory_Alloc();
        void Memory_Dealloc();

        Row rows;
        Col cols;
        T** data = nullptr;

};

template <typename T>
Matrix<T> operator+(const Matrix<T>&, const Matrix<T>&);

template <typename T>
Matrix<T> operator-(const Matrix<T>&, const Matrix<T>&);

template <typename T>
Matrix<T> operator*(const Matrix<T>&, const Matrix<T>&);

template <typename T>
Matrix<T> operator*(const Matrix<T>&, const double&);

template <typename T>
Matrix<T> operator*(const double&, const Matrix<T>&);

template <typename T>
Matrix<T> operator/(const Matrix<T>&, const double&);

#endif // GENTRYCRYPTO_MATH_LINALG_HPP