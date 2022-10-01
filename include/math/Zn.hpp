#ifndef ZN_HPP
#define ZN_HPP

template <int M>
struct Zn
{
	int Value;

	Zn(int value = 0);

	Zn operator - ();
	Zn operator ++();

	Zn operator + (Zn number);
	Zn operator +=(Zn number);

	Zn operator - (Zn number);
	Zn operator -=(Zn number);
	
	Zn operator * (Zn number);
	Zn operator *=(Zn number);	

	Zn operator + (int number);
	Zn operator +=(int number);

	Zn operator - (int number);
	Zn operator -=(int number);
	
	Zn operator * (int number);
	Zn operator *=(int number);	

	static int Congruence(int number);
};

template <int M>
int Zn<M>::Congruence(int number)
{
	int result = number % M;
	if (number < 0) result += M;

	return result;
}

template <int M>
Zn<M>::Zn(int value)
{
	Value = Zn<M>::Congruence(value);
}

template <int M>
Zn<M> Zn<M>::operator -()
{
	return Zn<M>(M - Value);
}

template <int M>
Zn<M> Zn<M>::operator ++()
{
	return Zn<M>(Value + 1);
}

template <int M>
Zn<M> Zn<M>::operator + (Zn number)
{
	return Zn<M>(Value + number.Value);
}

template <int M>
Zn<M> Zn<M>::operator +=(Zn number)
{
	Value = Zn<M>::Congruence(Value + number.Value);
	return *this;
}

template <int M>
Zn<M> Zn<M>::operator - (Zn number)
{
	return Zn<M>(Value - number.Value);
}

template <int M>
Zn<M> Zn<M>::operator -=(Zn number)
{
	Value = Zn<M>::Congruence(Value - number.Value);
	return *this;
}

template <int M>
Zn<M> Zn<M>::operator * (Zn number)
{
	return Zn<M>(Value * number.Value);
}

template <int M>
Zn<M> Zn<M>::operator *=(Zn number)
{
	Value = Zn<M>::Congruence(Value * number.Value);
	return *this;
}	

template <int M>
Zn<M> Zn<M>::operator + (int number)
{
	return Zn<M>(Value + number);
}

template <int M>
Zn<M> Zn<M>::operator +=(int number)
{
	Value = Zn<M>::Congruence(Value + number);
	return *this;
}

template <int M>
Zn<M> Zn<M>::operator - (int number)
{
	return Zn<M>(Value - number);
}

template <int M>
Zn<M> Zn<M>::operator -=(int number)
{
	Value = Zn<M>::Congruence(Value - number);
	return *this;
}

template <int M>
Zn<M> Zn<M>::operator * (int number)
{
	return Zn<M>(Value * number);
}

template <int M>
Zn<M> Zn<M>::operator *=(int number)
{
	Value = Zn<M>::Congruence(Value * number);
	return *this;
}	

#endif