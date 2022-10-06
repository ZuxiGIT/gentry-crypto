#include "../../include/math/linalg.hpp"
#include <stdexcept>

template <typename T>
Matrix<T>::Matrix() : rows(1), cols(1)
{
    Memory_Alloc();
    data[0][0] = 0;
}

template <typename T>
Matrix<T>::Matrix(Row rows, Col cols) : rows(rows), cols(cols)
{
    Memory_Alloc();
    for (int i = 0; i < rows.value; i++)
    {
        for (int j = 0; j < cols.value; j++)
        {
            data[i][j] = 0  ;
        }
    }
}

template <typename T>
Matrix<T>::Matrix(Row rows, Col cols, T** matrix) : rows(rows), cols(cols)
{
    Memory_Alloc();
    if (matrix != nullptr)
    {
        for (int i = 0; i < rows.value; i++)
        {
            for (int j = 0; j < cols.value; j++)
            {
                data[i][j] = matrix[i][j];
            }
        }
    }
}

template <typename T>
Matrix<T>::~Matrix()
{
    Memory_Dealloc();
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& obj) : rows(obj.rows.value), cols(obj.cols.value)
{
    Memory_Alloc();
    for (int i = 0; i < rows.value; i++)
    {
        std::copy(obj.data[i], obj.data[i] + cols.value, data[i]);
    }
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& obj)
{
    if (this != &obj)
    {
        if (rows.value != obj.rows.value || cols.value != obj.cols.value)
        {   
            // TODO: Replace it with creating an object
            Memory_Dealloc();
            rows.value = obj.rows.value;
            cols.value = obj.cols.value;
            Memory_Alloc();
        }

        for (int i = 0; i < rows.value; i++)
        {
            std::copy(obj.data[i], obj.data[i] + cols.value, data[i]);
        }             
    }
    return *this;
}

template <typename T>
Matrix<T>::Matrix(Matrix<T> && obj)
{
    std::swap<T**>(data, obj.data);
    std::swap<size_t>(cols.value, obj.cols.value);
    std::swap<size_t>(rows.value, obj.rows.value);
    obj.data = nullptr;
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>&& obj)
{
    if (this != &obj)
    {
        Memory_Dealloc();
        std::swap<T**>(data, obj.data);
        std::swap<size_t>(cols.value, obj.cols.value);
        std::swap<size_t>(rows.value, obj.rows.value);    
        obj.data = nullptr;     
    }
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator+=(const Matrix<T>& obj)
{
    if (rows.value != obj.rows.value || cols.value != obj.cols.value)
    {
        throw std::domain_error("Error: Matrix spaces don't match. Please fix it and try again.");
    }

    for (int i = 0; i < rows.value; i++)
    {
        for (int j = 0; j < cols.value; j++)
        {
            data[i][j] += obj.data[i][j];
        }
    }
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator-=(const Matrix<T>& obj)
{
    if (rows.value != obj.rows.value || cols.value != obj.cols.value)
    {
        throw std::domain_error("Error: Matrix spaces don't match. Please fix it and try again.");
    }

    for (int i = 0; i < rows.value; i++)
    {
        for (int j = 0; j < cols.value; j++)
        {
            data[i][j] -= obj.data[i][j];
        }
    }
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator*=(const Matrix<T>& obj)
{
    if (rows.value != obj.cols.value)
    {
        throw std::domain_error("Error: Matrix spaces don't match. Please fix input and try again.");
    }

    Matrix<T> temp(Row(rows.value), Col(obj.cols.value));
    for (int i = 0; i < temp.rows.value; i++)
    {
        for (int j = 0; j < temp.cols.value; j++)
        {
            for (int k = 0; k < cols.value; k++)
            {
                temp.data[i][j] += (data[i][k] * obj.data[k][j]);
            }
        }
    }
    return (*this = temp);
}

template <typename T>
Matrix<T> Matrix<T>::operator*=(const double& number)
{
    for (int i = 0; i < rows.value; i++)
    {
        for (int j = 0; j < cols.value; j++)
        {
            data[i][j] *= number;
        }
    }
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator/=(const double& number)
{
    for (int i = 0; i < rows.value; i++)
    {
        for (int j = 0; j < cols.value; j++)
        {
            data[i][j] /= number;
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
    Matrix<T> temp = Matrix<T>(Row(this->rows.value), Col(this->cols.value));
    for (int i = 0; i < temp.rows.value; i++)
    {
        for (int j = 0; j < temp.cols.value; j++)
        {
            temp.data[i][j] -= this->data[i][j];
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
void Matrix<T>::Memory_Alloc()
{
    data = new T* [rows.value];
    for (int i = 0; i < rows.value; i++)
    {
        data[i] = new T[cols.value];
    }
}

template <typename T>
void Matrix<T>::Memory_Dealloc()
{
    if (data != nullptr)
    {
        for (int i = 0; i < rows.value; i++)
        {
            delete[] data[i];
        }
        delete[] data;
    }
}

template <typename T>
T& Matrix<T>::operator()(size_t row, size_t col)
{
    if (row == 0 || col == 0) 
    {
        throw std::domain_error("Error: References to indices start with 1. Please fix it and try again.");
    }
    else if (row > rows.value || col > cols.value) 
    {
        throw std::domain_error("Error: Going out of the matrix boundaries. Please fix it and try again.");
    }
    return data[row - 1][col - 1];
}

template <typename T>
Matrix<T> Matrix<T>::Transpose()
{
    Matrix<T> temp(Row(cols.value), Col(rows.value));
    for (int i = 0; i < temp.rows.value; i++)
    {
        for (int j = 0; j < temp.cols.value; j++)
        {
            temp.data[i][j] = data[j][i];
        }
    }
    return temp;
}