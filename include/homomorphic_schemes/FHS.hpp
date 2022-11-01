#ifndef GENTRYCRYPTO_HOMOMORPHIC_FHS_HPP
#define GENTRYCRYPTO_HOMOMORPHIC_FHS_HPP

class FHS
{
    private:
        unsigned int sec_parameter  = 0; // lambda
	unsigned int subset_size    = 0; // s 
	unsigned int prec_parameter = 0; // p
	unsigned int polynomial_deg = 0; // d
	unsigned int coeff_bitsize  = 0; // t
	unsigned int lattice_dim    = 0; // n
	unsigned int bigset_size    = 0; // S
        unsigned int bigset_ratio   = 0; // R
        float        bdd_parameter  = 0; // mu
    public:
	/*
	 * There will be key generation, encryption, decryption functions
	 */

};

#endif //GENTRYCRYPTO_HOMOMORPHIC_FHS_HPP
