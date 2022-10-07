#ifndef GENTRYCRYPTO_MATH_MATRIX_HPP
#define GENTRYCRYPTO_MATH_MATRIX_HPP

#include "Vector.hpp"
template <typename T>
class Matrix
{
    private:
        unsigned int m_row = 0;
        unsigned int m_col = 0;
        T* m_data = nullptr;
    
    public:
        Matrix() = delete;    
        explicit Matrix(unsigned int row, unsigned int col);
        explicit Matrix(const Vector<T>);
        ~Matrix();
        
        Matrix(const Matrix&);
        Matrix& operator=(const Matrix&);
        
        Matrix(Matrix&&);
        Matrix& operator=(Matrix&&);     

        Matrix& operator+=(const Matrix&);
        Matrix& operator-=(const Matrix&);
        Matrix& operator*=(const Matrix&);

        Matrix& operator*=(const double&);
        Matrix& operator/=(const double&);

        Matrix operator+() const;
        Matrix operator-() const;

        T& operator()(unsigned int, unsigned int);

        Matrix transpose();
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

#endif // GENTRYCRYPTO_MATH_MATRIX_HPP