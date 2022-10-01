#ifndef LINALG_HPP
#define LINALG_HPP

template <int Length, class Object_Type>
struct Vector
{
	Object_Type * Values;

	Vector(Object_Type * values = NULL);
	~Vector();

	Vector operator - ();

	Vector operator + (Vector vec);
	Vector operator +=(Vector vec);

	Vector operator - (Vector vec);
	Vector operator -=(Vector vec);

	Vector operator * (int number);
	Vector operator *=(int number);

	Vector operator / (int number);
	Vector operator /=(int number);
};

template <int Width, int Height, class Object_Type>
struct Matrix
{
	Object_Type ** Values;

	Matrix(Object_Type ** values = NULL);
	~Matrix();

	Matrix operator - ();

	Matrix operator + (Matrix mat);
	Matrix operator +=(Matrix mat);

	Matrix operator - (Matrix mat);
	Matrix operator -=(Matrix mat);

	Matrix operator * (int number);
	Matrix operator *=(int number);

	Matrix operator / (int number);
	Matrix operator /=(int number);
};

#endif