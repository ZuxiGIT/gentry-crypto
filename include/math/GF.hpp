#ifndef GF_HPP
#define GF_HPP

template <int BM, int FM, int C = 0>
struct DegP
{
	static const int Value = DegP<BM, FM / BM, C + 1>::Value;
};

template <int BM, int C>
struct DegP<BM, 1, C>
{
	static const int Value = C;
};

template <int BM, int FM> 
struct GF
{
	int Value;

	GF(int value = 0);
	GF(Vector<Zn, DegP<BM, FM>::Value> vector);

	GF operator - ();
	GF operator ~ ();

	GF operator + (GF number);
	GF operator +=(GF number);

	GF operator - (GF number);
	GF operator -=(GF number);
	
	GF operator * (GF number);
	GF operator *=(GF number);

	GF operator + (Zn number);
	GF operator +=(Zn number);

	GF operator - (Zn number);
	GF operator -=(Zn number);
	
	GF operator * (Zn number);
	GF operator *=(Zn number);		

	GF operator + (int number);
	GF operator +=(int number);

	GF operator - (int number);
	GF operator -=(int number);
	
	GF operator * (int number);
	GF operator *=(int number);

	static int Congruence(int number);
	static int Congruence(Vector<Zn, DegP<BM, FM>::Value> vector);

	static Vector<Zn, DegP<BM, FM>::Value> NumberToVector(int number);
	static int VectorToNumber(Vector<Zn, DegP<BM, FM>::Value> vector);
};

template <int BM, int FM>
int GF<BM, FM>::Congruence(int number)
{
	return GF<BM, FM>::Congruence(GF<BM, FM>::NumberToVector(number));
}

template <int BM, int FM>
int GF<BM, FM>::Congruence(Vector<Zn, DegP<BM, FM>::Value> vector)
{
	Vector<Zn, DegP<BM, FM>::Value> field_module = GF<BM, FM>::NumberToVector(FM);

	for (int i = 0; i < DegP<BM, FM>::Value; i++)
	{
		
	}
}

template <int BM, int FM>
Vector<Zn, DegP<BM, FM>::Value> GF<BM, FM>::NumberToVector(int number)
{
	Vector<Zn, DegP<BM, FM>::Value> result;

	int digit = 0;

	for (int i = 0; number; i++)
	{
		digit = number % BM;

		result[DegP<BM, FM>::Value - i] = digit;

		number /= BM;
	}

	return result;
}

template <int BM, int FM>
int GF<BM, FM>::VectorToNumber(Vector<Zn, DegP<BM, FM>::Value> vector)
{
	int result = 0;

	int coeff = 1;

	for (int i = 0; i < DegP<BM, FM>::Value; i++)
	{
		result += vector[i];
		result *= coeff;
		coeff *= BM;
	}

	return result;
}

template <int BM, int FM> 
GF<BM, FM>::GF(int value)
{
	Value = GF<BM, FM>::Congruence(value);
}

template <int BM, int FM> 
GF<BM, FM>::GF(Vector<Zn, DegP<BM, FM>::Value> vector)
{
	Value = GF<BM, FM>::VectorToNumber(value);
}

template <int BM, int FM>
GF GF<BM, FM>::operator +(GF number)
{
	return GF(GF<BM, FM>::NumberToVector(Value) + GF<BM, FM>::NumberToVector(number.Value));
}

template <int BM, int FM>
GF GF<BM, FM>::operator +=(GF number)
{
	Value = GF<BM, FM>::Congruence(GF<BM, FM>::NumberToVector(Value) + GF<BM, FM>::NumberToVector(number.Value));
	return *this;
}

template <int BM, int FM>
GF GF<BM, FM>::operator - (GF number)
{
	return GF(GF<BM, FM>::NumberToVector(Value) - GF<BM, FM>::NumberToVector(number.Value));
}

template <int BM, int FM>
GF GF<BM, FM>::operator -=(GF number)
{
	Value = GF<BM, FM>::Congruence(GF<BM, FM>::NumberToVector(Value) + GF<BM, FM>::NumberToVector(number.Value));
	return *this;
}

template <int BM, int FM>
GF GF<BM, FM>::operator * (GF number)
{

}

template <int BM, int FM>
GF GF<BM, FM>::operator *=(GF number)
{

	return *this;
}

#endif
