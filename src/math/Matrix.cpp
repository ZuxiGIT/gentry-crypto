//#include "../../include/math/linalg.hpp"
#include "Matrix.hpp"
#include <stdexcept>
#include <utility>

template <typename T>
Matrix<T>::Matrix(unsigned int row, unsigned int col)
: m_row(row),
  m_col(col),
  m_data(new T[m_col * m_row])
{}

template <typename T>
Matrix<T>::~Matrix()
{
    delete[] m_data; 
}

template <typename T>
Matrix<T>::Matrix(const Vector<T> v)
: Matrix(v.getSize(), 1)
{
    std::copy(v.getData(), v.getData() + v.getSize(), m_data);
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& rhs)
: Matrix(rhs.m_row, rhs.m_col)
{
    std::copy(rhs.m_data, rhs.m_data + rhs.m_row * rhs.m_col, m_data);
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& rhs)
{
    if (this == &rhs)
        return *this;

    if (m_row != rhs.m_row || m_col != rhs.m_col)
    {   
        delete[] m_data;
        m_row = rhs.m_row;
        m_col = rhs.m_col;
        m_data = new T[m_col * m_row];
    }

    std::copy(rhs.m_data, rhs.m_data + rhs.m_row * rhs.m_col, m_data);
    return *this;
}

template <typename T>
Matrix<T>::Matrix(Matrix<T>&& rhs)
{
    m_data = rhs.m_data;
    m_col = rhs.m_col;
    m_row = rhs.m_row;
    rhs.m_data = nullptr;
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>&& rhs)
{
    delete[] m_data;
    m_data = rhs.m_data;
    m_col = rhs.m_col;
    m_row = rhs.m_row;
    rhs.m_data = nullptr;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& obj)
{
    if (m_row != obj.m_row|| m_col != obj.m_col)
    {
        throw std::domain_error("Error: Matrix spaces don't match. Please fix it and try again.");
    }

    for (int i = 0; i < m_row; i++)
    {
        for (int j = 0; j < m_col; j++)
        {
            m_data[i][j] += obj.m_data[i][j];
        }
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& obj)
{
    if (m_row != obj.m_row || m_col != obj.m_col)
    {
        throw std::domain_error("Error: Matrix spaces don't match. Please fix it and try again.");
    }

    for (int i = 0; i < m_row; i++)
    {
        for (int j = 0; j < m_col; j++)
        {
            m_data[i][j] -= obj.m_data[i][j];
        }
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& obj)
{
    if (m_row != obj.m_col)
    {
        throw std::domain_error("Error: Matrix spaces don't match. Please fix input and try again.");
    }

    Matrix<T> temp(obj.m_row, obj.m_col);
    for (int i = 0; i < temp.m_row; i++)
    {
        for (int j = 0; j < temp.m_col; j++)
        {
            for (int k = 0; k < m_col; k++)
            {
                temp.m_data[i][j] += (m_data[i][k] * obj.m_data[k][j]);
            }
        }
    }
    return (*this = std::move(temp));
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(const double& number)
{
    for (int i = 0; i < m_row; i++)
    {
        for (int j = 0; j < m_col; j++)
        {
            m_data[i][j] *= number;
        }
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator/=(const double& number)
{
    for (int i = 0; i < m_row; i++)
    {
        for (int j = 0; j < m_col; j++)
        {
            m_data[i][j] /= number;
        }
    }
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator+() const
{
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator-() const
{
    Matrix<T> temp(m_row, m_col);
    for (int i = 0; i < m_row; i++)
    {
        for (int j = 0; j < m_col; j++)
        {
            temp.m_data[i][j] = -m_data[i][j];
        }
    }
    return temp;
}

template <typename T>
Matrix<T> operator+(const Matrix<T>& l_obj, const Matrix<T>& r_obj)
{
    Matrix<T> temp(l_obj);
    return (temp += r_obj);   
}

template <typename T>
Matrix<T> operator-(const Matrix<T>& l_obj, const Matrix<T>& r_obj)
{
    Matrix<T> temp(l_obj);
    return (temp -= r_obj);    
}

template <typename T>
Matrix<T> operator*(const Matrix<T>& l_obj, const Matrix<T>& r_obj)
{
    Matrix<T> temp(l_obj);
    return (temp *= r_obj);   
}

template <typename T>
Matrix<T> operator*(const Matrix<T>& obj, const double& number)
{
    Matrix<T> temp(obj);
    return (temp *= number);   
}

template <typename T>
Matrix<T> operator*(const double& number, const Matrix<T>& obj)
{
    return obj * number;   
}

template <typename T>
Matrix<T> operator/(const Matrix<T>& obj, const double& number)
{
    Matrix<T> temp(obj);
    return (temp /= number);
}

template <typename T>
T& Matrix<T>::operator()(unsigned int row, unsigned int col)
{
    if (row == 0 || col == 0) 
    {
        throw std::domain_error("Error: References to indices start with 1. Please fix it and try again.");
    }
    else if (row > m_row || col > m_col) 
    {
        throw std::domain_error("Error: Going out of the matrix boundaries. Please fix it and try again.");
    }
    return m_data[row - 1][col - 1];
}

template <typename T>
Matrix<T> Matrix<T>::transpose()
{
    Matrix<T> temp(m_row, m_col);
    for (int i = 0; i < m_row; i++)
    {
        for (int j = 0; j < m_col; j++)
        {
            temp.m_data[i][j] = m_data[j][i];
        }
    }
    return temp;
}