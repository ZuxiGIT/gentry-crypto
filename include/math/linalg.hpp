#ifndef LINALG_HPP
#define LINALG_HPP

template <class OT, int L>
struct Vector
{
	OT * Values;

	Vector(OT * values = NULL);
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

template <class OT, int W, int H>
struct Matrix
{
	OT ** Values;

	Matrix(OT ** values = NULL);
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

template <int OT, int L>
Vector<OT, L>::Vector(ObjectT * values)
{
	Values = new OT[L];

	for (int i = 0; i < L; i++)
	{
		Values[i] = values[i];
	}
}

template <int OT, int L>
Vector<OT, L>::~Vector()
{
	delete[] Values;
}

#endif