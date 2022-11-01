#ifndef GENTRYCRYPTO_MATH_POLYNOMIAL_HPP
#define GENTRYCRYPTO_MATH_POLYNOMIAL_HPP

#include <vector>
#include <algorithm>
#include <math.h> 

#include <stdexcept>
#include <utility>

/* Class Polynomial
 * Created for the mathematical representation of polynomials. Objects - polynomials of a given
 * degree (which can be specified in the constructor or implicitly specified in the passed vector). 
 * 
 * The form of the polynomial notation is written as follows: a_{n-1} * x^{n-1} + ... + a_1 * x + a_0
 * The main purpose of the class is the multiplication/division mod operations of the polynomial.
*/
template <typename T>
class Polynomial
{
    private:
        unsigned int m_degree = 0;
        std::vector<T> m_coeff;
        void syntheticDiv(Polynomial<T>&, const Polynomial<T>&);

    public:
        Polynomial() = delete;
        ~Polynomial();
        explicit Polynomial(unsigned int);
        explicit Polynomial(std::vector<T>);

        Polynomial(const Polynomial&);
        Polynomial& operator=(const Polynomial&);

        Polynomial(Polynomial&&);
        Polynomial& operator=(Polynomial&&); 

        Polynomial& operator+=(const Polynomial&);
        Polynomial& operator-=(const Polynomial&);
        Polynomial& operator*=(const Polynomial&);
        Polynomial& operator/=(const Polynomial&);
        Polynomial& operator%=(const Polynomial&);

        Polynomial operator+();
        Polynomial operator-();

        T& operator[](unsigned int) const; 
        T  operator()(T);
        
	unsigned int size()   const;
	unsigned int degree() const;
};

template <typename T>
void Polynomial<T>::syntheticDiv(Polynomial<T>& dividend, const Polynomial<T>& divider)
{
    unsigned int divider_deg  = divider.m_degree;
    unsigned int dividend_deg = dividend.m_degree;
    unsigned int quotient_deg = dividend_deg - divider_deg + 1;
    
    T normalizer = divider[0];
    T coeff;
    for (int i = 0; i < quotient_deg; i++)
    {
        dividend[i] /= normalizer;
        coeff = dividend[i];
        if (coeff != 0)
        {
            for (int j = 1; j < divider_deg; j++)
            {
                dividend[i + j] -= divider[j] * coeff;
            }
        }
    }
}

template <typename T>
Polynomial<T>::Polynomial(unsigned int degree)
: m_coeff(std::vector<T>(degree + 1)), m_degree(degree)
{}

template <typename T>
Polynomial<T>::Polynomial(std::vector<T> coeff)
: m_coeff(coeff), m_degree(coeff.size() - 1)
{}

template <typename T>
Polynomial<T>::~Polynomial()
{}

template <typename T>
Polynomial<T>::Polynomial(const Polynomial<T>& rhs)
: m_coeff(rhs.m_coeff), m_degree(rhs.m_coeff.size() - 1)
{}

template <typename T>
Polynomial<T>& Polynomial<T>::operator=(const Polynomial<T>& rhs)
{
    if (this != &rhs)
    {
        m_coeff  = rhs.m_coeff;
        m_degree = rhs.m_degree;
    }
    return *this;
}

template <typename T>
Polynomial<T>::Polynomial(Polynomial<T>&& rhs)
: m_coeff(std::move(rhs.m_coeff)), m_degree(std::move(m_degree))
{}

template <typename T>
Polynomial<T>& Polynomial<T>::operator=(Polynomial<T>&& rhs)
{
    m_coeff  = std::move(rhs.m_coeff);
    m_degree = std::move(rhs.m_degree);
    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator+=(const Polynomial<T>& rhs)
{
    unsigned int max_size = (m_degree >= rhs.m_degree ? m_degree : rhs.m_degree) + 1;
    std::vector<T> temp = std::vector<T>(max_size);
    for (unsigned int i = max_size - 1; i >= 0; i--)
    {
        temp[i] = m_coeff[i] + rhs.m_coeff[i];
    }
    return (m_coeff = temp);      
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator-=(const Polynomial<T>& rhs)
{
    unsigned int max_size = (m_degree > rhs.m_degree ? m_degree : rhs.m_degree) + 1;
    std::vector<T> temp = std::vector<T>(max_size);
    for (unsigned int i = max_size - 1; i >= 0; i--)
    {
        temp[i] = m_coeff[i] - rhs.m_coeff[i];
    }
    return (m_coeff = temp);      
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator*=(const Polynomial<T>& rhs)
{ 
    unsigned int lhs_deg = m_degree;
    unsigned int rhs_deg = rhs.m_degree;
    unsigned int prod_deg = lhs_deg + rhs_deg - 1;
    
    Polynomial<T> temp(prod_deg);
    for (int i = lhs_deg; i >= 0; i--)
    {
        for (int j = rhs_deg; j >= 0; j--)
        {
            temp[i + j] += m_coeff[i] * rhs.m_coeff[j];
        }
    }
    return (*this = std::move(temp));
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator/=(const Polynomial<T>& rhs)
{    
    if (m_degree < rhs.m_degree)
        return *this;
    
    if ((rhs.m_degree == 0) && (rhs.m_coeff[0] == 0))    
        throw std::domain_error("Error: Division by zero. Please fix it and try again.");
	    
    syntheticDiv((*this), rhs);
    unsigned int quotient_deg = m_degree - rhs.m_degree + 1;
    for (; m_degree > quotient_deg; m_degree--)
    {
        m_coeff.pop_back();
    }
    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator%=(const Polynomial<T>& rhs)
{
    if (m_degree < rhs.m_degree)
        return *this;
    
    if ((rhs.m_degree == 0) && (rhs.m_coeff[0] == 0))    
        throw std::domain_error("Error: Division by zero. Please fix it and try again.");

    syntheticDiv((*this), rhs);
    unsigned int rhs_deg = rhs.m_degree;
    for (; m_degree >= rhs_deg; m_degree--)
    {
        m_coeff.erase(m_coeff.begin());
    }
    return *this;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator+()
{
    return *this;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator-()
{
    if (m_degree == 0)
        return *this;

    for (int i = 1; i < m_degree + 1; i += 2) 
    {
        m_coeff[i] = -m_coeff[i];
    }
    return *this;
}

template <typename T>
Polynomial<T> operator+(const Polynomial<T>& lhs, const Polynomial<T>& rhs)
{
    Polynomial<T> temp(lhs);
    return (temp += rhs);   
}

template <typename T>
Polynomial<T> operator-(const Polynomial<T>& lhs, const Polynomial<T>& rhs)
{
    Polynomial<T> temp(lhs);
    return (temp -= rhs);   
}


template <typename T>
Polynomial<T> operator*(const Polynomial<T>& lhs, const Polynomial<T>& rhs)
{
    Polynomial<T> temp(lhs);
    temp *= rhs;
    return temp;   
}

template <typename T>
Polynomial<T> operator/(const Polynomial<T>& lhs, const Polynomial<T>& rhs)
{
    Polynomial<T> temp(lhs);
    return (temp /= rhs);   
}

template <typename T>
Polynomial<T> operator%(const Polynomial<T>& lhs, const Polynomial<T>& rhs)
{
    Polynomial<T> temp(lhs);
    return (temp %= rhs);   
}

template <typename T>
T& Polynomial<T>::operator[](unsigned int coeff_num) const
{
    return const_cast<T&>(m_coeff[coeff_num]);
}

template <typename T>
T Polynomial<T>::operator()(T x)
{
    T value = 0;
    for (int i = m_degree; i >= 0; i--) 
    {
        value += pow(x, i) * m_coeff[i];
    }
    return value;
}

template <typename T>
unsigned int Polynomial<T>::degree() const
{
    return m_degree;
}
#endif //GENTRYCRYPTO_MATH_POLYNOMIAL_HPP