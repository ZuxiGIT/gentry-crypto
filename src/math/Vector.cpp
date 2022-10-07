#include "Vector.hpp"
#include <stdexcept>
#include <utility>

template <typename T>
Vector<T>::Vector(unsigned int size)
: v_size(size), 
  v_data(new T[v_size])
{}

template <typename T>
Vector<T>::~Vector()
{
    delete[] v_data; 
}

template <typename T>
Vector<T>::Vector(const Vector<T>& rhs)
: Vector(rhs.v_size)
{
    std::copy(rhs.v_data, rhs.v_data + rhs.v_size, v_data);
}

template <typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& rhs)
{
    if (this == &rhs)
        return *this;

    if (v_size != rhs.v_size)
    {   
        delete[] v_data;
        v_size = rhs.v_size;
        v_data = new T[v_size];
    }

    std::copy(rhs.v_data, rhs.v_data + rhs.v_size, v_data);
    return *this;
}

template <typename T>
Vector<T>::Vector(Vector<T>&& rhs)
{
    v_data = rhs.v_data;
    v_size = rhs.v_size;
    rhs.v_data = nullptr;
}

template <typename T>
Vector<T>& Vector<T>::operator=(Vector<T>&& rhs)
{
    delete[] v_data;
    v_data = rhs.v_data;
    v_size = rhs.v_size;
    rhs.v_data = nullptr;
}

template <typename T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& obj)
{
    if (v_size != obj.v_size)
    {
        throw std::domain_error("Error: Vector spaces don't match. Please fix it and try again.");
    }

    for (int i = 0; i < v_size; i++)
    {
        v_data[i] += obj.v_data[i];
    }
    return *this;
}

template <typename T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& obj)
{
    if (v_size != obj.v_size)
    {
        throw std::domain_error("Error: Vector spaces don't match. Please fix it and try again.");
    }

    for (int i = 0; i < v_size; i++)
    {
        v_data[i] -= obj.v_data[i];
    }
    return *this;
}

template <typename T>
Vector<T>& Vector<T>::operator*=(const Vector<T>& obj)
{
    if (v_size != obj.v_size)
    {
        throw std::domain_error("Error: Vector spaces don't match. Please fix input and try again.");
    }

    Vector<T> temp(obj.v_size);
    for (int i = 0; i < v_size; i++)
    {
        temp.v_data[i] = v_data[i] * obj.v_data[i];
    }
    return (*this = std::move(temp));
}

template <typename T>
Vector<T>& Vector<T>::operator*=(const double& number)
{
    for (int i = 0; i < v_size; i++)
    {
        v_data[i] *= number;
    }
    return *this;
}

template <typename T>
Vector<T>& Vector<T>::operator/=(const double& number)
{
    for (int i = 0; i < v_size; i++)
    {
        v_data[i] /= number;
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
    Vector<T> temp(v_size);
    for (int i = 0; i < v_size; i++)
    {
        temp[i] -= v_data[i];
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
T& Vector<T>::operator()(unsigned int elem)
{
    if (v_size == 0) 
    {
        throw std::domain_error("Error: References to indices start with 1. Please fix it and try again.");
    }
    else if (elem > v_size) 
    {
        throw std::domain_error("Error: Going out of the Vector boundaries. Please fix it and try again.");
    }
    return v_data[elem - 1];
}

template <typename T>
unsigned int Vector<T>::getSize()
{
    return v_size;
}

template <typename T>
T* Vector<T>::getData()
{
    return v_data;
}
