#ifndef GF_HPP
#define GF_HPP

template <int BM, int FM> 
struct GF
{
	int Value;

	GF(int value = 0);
	GF(Vector vector);

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

	static Vector NumberToVector(int number);
	static int VectorToNumber(Vector vector);
};

template <int BM, int FM>
int GF<int BM, int FM>::Congruence(int number)
{
	//Not implemented
}

template <int BM, int FM>
Vector GF<int BM, int FM>::NumberToVector(int number)
{
	//Not implemented
}

template <int BM, int FM>
int GF<int BM, int FM>::VectorToNumber(Vector vector)
{
	//Not implemented
}

template <int BM, int FM> 
GF<BM, FM>::GF(int value)
{
	Value = GF<BM, FM>::Congruence(value);
}
template <int BM, int FM>
GF GF<BM, FM>::operator +(GF number)
{
	return GF(GF<BM, FM>::NumberToVector(Value) + GF<BM, FM>::NumberToVector(number.Value));
}

GF operator +=(GF number);

GF operator - (GF number);
GF operator -=(GF number);

GF operator * (GF number);
GF operator *=(GF number);

#endif
