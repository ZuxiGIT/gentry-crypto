#ifndef GF_HPP
#define GF_HPP

template <int Base_Modulus, int Field_Modulus> 
struct GF
{
	int Value;

	GF(int value = 0)
	{
		Value = value;
	}

	GF operator - ();

	GF operator + (GF number);
	GF operator +=(GF number);

	GF operator - (GF number);
	GF operator -=(GF number);

	GF operator * (int number);
	GF operator *=(int number);

	GF operator ~ ();
	
	GF operator * (GF number);
	GF operator *=(GF number);	
};

#endif
