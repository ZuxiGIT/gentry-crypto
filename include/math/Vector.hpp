#ifndef GENTRYCRYPTO_MATH_VECTOR_HPP
#define GENTRYCRYPTO_MATH_VECTOR_HPP

#include "Matrix.hpp"
#include <stdexcept>
#include <utility>

template <typename T>
class Vector
{
    private:
        unsigned int m_size = 0;
        T* m_data = nullptr;
    
    public:
        Vector() = delete;
        explicit Vector(unsigned int size);
        ~Vector();

        Vector(const Vector&);
        Vector& operator=(const Vector&);
        
        Vector(Vector&&);
        Vector& operator=(Vector&&);     

        Vector& operator+=(const Vector&);
        Vector& operator-=(const Vector&);
        Vector& operator*=(const Vector&);

        Vector& operator*=(const double&);
        Vector& operator/=(const double&);

        Vector operator+() const;
        Vector operator-() const;

        T& operator()(unsigned int) const; 
        unsigned int getSize() const;
};

template <typename T>
Vector<T>::Vector(unsigned int size)
: m_size(size), 
  m_data(new T[m_size])
{}

template <typename T>
Vector<T>::~Vector()
{
    delete[] m_data; 
}

template <typename T>
Vector<T>::Vector(const Vector<T>& rhs)
: Vector(rhs.m_size)
{
    std::copy(rhs.m_data, rhs.m_data + rhs.m_size, m_data);
}

template <typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& rhs)
{
    if (this == &rhs)
        return *this;

    if (m_size != rhs.m_size)
    {   
        delete[] m_data;
        m_size = rhs.m_size;
        m_data = new T[m_size];
    }

    std::copy(rhs.m_data, rhs.m_data + rhs.m_size, m_data);
    return *this;
}

template <typename T>
Vector<T>::Vector(Vector<T>&& rhs)
{
    m_data = rhs.m_data;
    m_size = rhs.m_size;
    rhs.m_data = nullptr;
}

template <typename T>
Vector<T>& Vector<T>::operator=(Vector<T>&& rhs)
{
    delete[] m_data;
    m_data = rhs.m_data;
    m_size = rhs.m_size;
    rhs.m_data = nullptr;
}

template <typename T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& obj)
{
    if (m_size != obj.m_size)
    {
        throw std::domain_error("Error: Vector spaces don't match. Please fix it and try again.");
    }

    for (int i = 0; i < m_size; i++)
    {
        m_data[i] += obj.m_data[i];
    }
    return *this;
}

template <typename T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& obj)
{
    if (m_size != obj.m_size)
    {
        throw std::domain_error("Error: Vector spaces don't match. Please fix it and try again.");
    }

    for (int i = 0; i < m_size; i++)
    {
        m_data[i] -= obj.m_data[i];
    }
    return *this;
}

template <typename T>
Vector<T>& Vector<T>::operator*=(const Vector<T>& obj)
{
    if (m_size != obj.m_size)
    {
        throw std::domain_error("Error: Vector spaces don't match. Please fix input and try again.");
    }

    Vector<T> temp(obj.m_size);
    for (int i = 0; i < m_size; i++)
    {
        temp.m_data[i] = m_data[i] * obj.m_data[i];
    }
    return (*this = std::move(temp));
}

template <typename T>
Vector<T>& Vector<T>::operator*=(const double& number)
{
    for (int i = 0; i < m_size; i++)
    {
        m_data[i] *= number;
    }
    return *this;
}

template <typename T>
Vector<T>& Vector<T>::operator/=(const double& number)
{
    for (int i = 0; i < m_size; i++)
    {
        m_data[i] /= number;
    }
    return *this;
}

template <typename T>
Vector<T> Vector<T>::operator+() const
{
    return *this;
}

template <typename T>
Vector<T> Vector<T>::operator-() const
{
    Vector<T> temp(m_size);
    for (int i = 0; i < m_size; i++)
    {
        temp[i] -= m_data[i];
    }
    return temp;
}

template <typename T>
Vector<T> operator+(const Vector<T>& l_obj, const Vector<T>& r_obj)
{
    Vector<T> temp(l_obj);
    return (temp += r_obj);   
}

template <typename T>
Vector<T> operator-(const Vector<T>& l_obj, const Vector<T>& r_obj)
{
    Vector<T> temp(l_obj);
    return (temp -= r_obj);    
}

template <typename T>
Vector<T> operator*(const Vector<T>& l_obj, const Vector<T>& r_obj)
{
    Vector<T> temp(l_obj);
    return (temp *= r_obj);   
}

template <typename T>
Vector<T> operator*(const Vector<T>& l_obj, const Matrix<T>& r_obj)
{
    if (l_obj.getSize() != r_obj.getColSize())
    {
        throw std::domain_error("Error: Spaces between Vector and Matrix don't match. Please fix input and try again.");
    }

    Vector<T> temp(l_obj.getSize());
    for (unsigned int i = 1; i < r_obj.getRowSize() + 1; i++)
    {
        temp(i) = 0;
        for (unsigned int j = 1; j < l_obj.getSize() + 1; j++)
        {
            temp(i) += l_obj(j) * r_obj(j, i);
        }
    }
    return temp;
}

template <typename T>
Vector<T> operator*(const Matrix<T>& l_obj, const Vector<T>& r_obj)
{
    if (l_obj.getRowSize() != r_obj.getSize())
    {
        throw std::domain_error("Error: Spaces between Matrix and Vector don't match. Please fix input and try again.");
    }

    Vector<T> temp(r_obj.getSize());
    for (unsigned int i = 1; i < l_obj.getColSize() + 1; i++)
    {
        temp(i) = 0;
        for(unsigned int j = 1; j < r_obj.getSize() + 1; j++)
        {
            temp(i) += l_obj(i, j) * r_obj(j);
        }
    }
    return temp; 
}

template <typename T>
Vector<T> operator*(const Vector<T>& obj, const double& number)
{
    Vector<T> temp(obj);
    return (temp *= number);   
}

template <typename T>
Vector<T> operator*(const double& number, const Vector<T>& obj)
{
    return obj * number;   
}

template <typename T>
Vector<T> operator/(const Vector<T>& obj, const double& number)
{
    Vector<T> temp(obj);
    return (temp /= number);
}

template <typename T>
T& Vector<T>::operator()(unsigned int elem) const
{
    if (m_size == 0) 
    {
        throw std::domain_error("Error: References to indices start with 1. Please fix it and try again.");
    }
    else if (elem > m_size) 
    {
        throw std::domain_error("Error: Going out of the Vector boundaries. Please fix it and try again.");
    }
    return m_data[elem - 1];
}

template <typename T>
unsigned int Vector<T>::getSize() const
{
    return m_size;
}

#endif // GENTRYCRYPTO_MATH_VECTOR_HPP
