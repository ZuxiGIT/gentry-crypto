#ifndef GENTRYCRYPTO_MATH_VECTOR_HPP
#define GENTRYCRYPTO_MATH_VECTOR_HPP

template <typename T>
class Vector
{
    private:
        unsigned int v_size = 0;
        T* v_data = nullptr;
    
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

        T& operator()(unsigned int);

        unsigned int getSize();
        T* getData();
};

template <typename T>
Vector<T> operator+(const Vector<T>&, const Vector<T>&);

template <typename T>
Vector<T> operator-(const Vector<T>&, const Vector<T>&);

template <typename T>
Vector<T> operator*(const Vector<T>&, const Vector<T>&);

template <typename T>
Vector<T> operator*(const Vector<T>&, const double&);

template <typename T>
Vector<T> operator*(const double&, const Vector<T>&);

template <typename T>
Vector<T> operator/(const Vector<T>&, const double&);

#endif // GENTRYCRYPTO_MATH_VECTOR_HPP