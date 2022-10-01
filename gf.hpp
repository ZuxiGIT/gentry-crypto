#ifndef GF_H
#define GF_H

template <int base_modulus, int field_modulus> 
class GF
{
	int number;

	GF(int number = 0);
	GF(Vector<int> polynomial);
	~GF();

	GF operator -();

	GF operator +(GF number);
	GF operator -(GF number);

	GF operator *(int number);

	GF operator *(GF number);	
};

#endif
