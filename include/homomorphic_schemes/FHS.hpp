#ifndef GENTRYCRYPTO_HOMOMORPHIC_FHS_HPP
#define GENTRYCRYPTO_HOMOMORPHIC_FHS_HPP

/*
 * This project is made for the learning purpose of studying the
 * Smart-Vercauteren and Gentry-Halevi schemes for the subject of
 * Information Security. There is no purpose for commercial use.
 *
 * Code is based Craig Gentry's IBM implementation and the article
 * "Implementing a fully homomorphic Gentry-Halevi encryption scheme"
 * by Gentry-Halevi, IBM Research February 4, 2011
 */

#include "SWHS.hpp"
#include <vector>

class PKblock
{
    public:
        NTL::ZZ x;         // the sequence of elements is xi = x*2^i mod d
        long idx;          // index of the xi that belongs to the sparse subset
        long operator==(const PKblock& other) const
        {
            return (x == other.x && idx == other.idx);
        }
};

typedef std::vector<PKblock> PKblockv;
class FHS: private SWHS
{
    private:
        unsigned int sec_parameter    = 0; // lambda
        unsigned int subset_size      = 0; // s
        unsigned int prec_parameter   = 0; // p
        unsigned int squash_poly_deg  = 0; // d
        unsigned int bigset_size      = 0; // S
        unsigned int bigset_log_ratio = 0; // R
        float        bdd_parameter    = 0; // mu

        PKblockv     pkBlocks;
        NTL::vec_ZZ  ctxts;

        unsigned long mChoose2(unsigned long);
        void encodeIndex(unsigned long&, unsigned long&, unsigned long, unsigned long);

    public:
        using SWHS::SWHS;
        FHS() = delete;
        ~FHS();
        explicit FHS::FHS(unsigned int, unsigned int, unsigned int, unsigned int,
                          unsigned int, unsigned int, unsigned int);

        FHS(const FHS&);
        FHS& operator=(const FHS&);

        FHS(FHS&&);
        FHS& operator=(FHS&&);

        void keyGen();
};

unsigned long FHS::mChoose2(unsigned long S)
{
    return static_cast<unsigned long>(ceil(2 * sqrt(static_cast<double>(S))));
}

void FHS::encodeIndex(unsigned long& j1, unsigned long& j2,
                      unsigned long  i,  unsigned long  m)
{}

FHS::FHS(unsigned int n, unsigned int t, unsigned int s, unsigned int p,
         unsigned int d, unsigned int S, unsigned int R)
: SWHS(n, t),
  subset_size(s),
  prec_parameter(p),
  squash_poly_deg(d),
  bigset_size(S),
  bigset_log_ratio(R)
{}

FHS::~FHS()
{}

FHS::FHS(const FHS& rhs)
: lattice_log_dim(rhs.lattice_log_dim),
  coeff_bitsize(rhs.coeff_bitsize)
  subset_size(rhs.subset_size),
  prec_parameter(rhs.prec_parameter),
  squash_poly_deg(rhs.squash_poly_deg),
  bigset_size(rhs.bigset_size),
  bigset_log_ratio(rhs.bigset_log_ratio)
{}

FHS& FHS::operator=(const FHS& rhs)
{
    if (this != &rhs)
    {
        lattice_log_dim  = rhs.lattice_log_dim;
        coeff_bitsize    = rhs.coeff_bitsize;
        subset_size      = rhs.subset_size;
        prec_parameter   = rhs.prec_parameter;
        squash_poly_deg  = rhs.squash_poly_deg;
        bigset_size      = rhs.bigset_size;
        bigset_log_ratio = rhs.bigset_log_ratio;
    }
    return *this;
}

FHS::FHS(FHS&& rhs)
: lattice_log_dim(std::move(rhs.lattice_log_dim)),
  coeff_bitsize(std::move(rhs.coeff_bitsize)),
  subset_size(std::move(rhs.subset_size)),
  prec_parameter(std::move(rhs.prec_parameter),
  squash_poly_deg(std::move(rhs.squash_poly_deg),
  bigset_size(std::move(rhs.bigset_size),
  bigset_log_ratio(std::move(rhs.bigset_log_ratio),
{}

FHS& FHS::operator=(FHS&& rhs)
{
    lattice_log_dim  = std::move(rhs.lattice_log_dim);
    coeff_bitsize    = std::move(rhs.coeff_bitsize);
    subset_size      = std::move(rhs.subset_size);
    prec_parameter   = std::move(rhs.prec_parameter);
    squash_poly_deg  = std::move(rhs.squash_poly_deg);
    bigset_size      = std::move(rhs.bigset_size);
    bigset_log_ratio = std::move(rhs.bigset_log_ratio);
    return *this;
}

void FHS::keyGen()
{
    SWHS::keyGen();
    {
        NTL::ZZ sum, factor;
        pkBlocks.resize(subset_size);

        sum = 0;
        for (unsigned int i = 0; i < subset_size - 1; i++)
        {
            /*
             * Selecting a random element (within [0, d - 1])
             * and a random index (within [0, S - 1])
             */

            NTL::RandomBnd(pkBlocks[i].x, lattice_det);
            pkBlocks[i].idx = NTL::RandomBnd(bigset_size);

            // Adding R^{idx} * x to the sum
            NTL::power2(factor, pkBlocks[i].idx * bigset_log_ratio); // R as a power of index
            if (factor > lattice_det)
            {
                factor %= lattice_det;
            }
            NTL::MulMod(factor, pkBlocks[i].x, factor, lattice_det);
            NTL::AddMod(sum, sum, factor, lattice_det);
        }

        /*
         * Choose the last progression so that the
         * sum equals (\sum sigma_i * x_i = w)
         */

        int i = subset_size - 1;
        pkBlocks[i].idx = NTL::RandomBnd(bigset_size);
        sum = NTL::SubMod(private_key, sum, lattice_det);
        NTL::power2(factor, pkBlocks[i].idx * bigset_log_ratio);
        if (factor > lattice_det)
        {
            factor %= lattice_det;
        }

        factor = NTL::InvMod(factor, lattice_det);
        pkBlocks[i].x = NTL::MulMod(sum, factor, lattice_det);

        // Check \sum sigma_i * x_i = w
        sum = 0;
        for (unsigned int i = 0; i < subset_size; i++)
        {
            NTL::power2(factor, pkBlocks[i].idx * bigset_log_ratio);
            if (factor > lattice_det)
            {
                factor %= lattice_det;
            }

            NTL::MulMod(factor, pkBlocks[i].x, factor, lattice_det);
            NTL::AddMod(sum, sum, factor, lattice_det);
        }
        if (sum != private_key) NTL::Error("[Error] The amount is not equal to the private key");
    }

    // Encrypting the characteristic vector
    unsigned long nCtxts = mChoose2(bigset_size);

    // Initialize bits to all zero, then set some of them to one
    unsigned int bits[nCtxts * subset_size];
    memset(bits, 0, sizeof(bits));
    for (unsigned int i = 0; i < subset_size; i++)
    {
        unsigned long j1, j2;
        // Let j1,j2 be the idx'th pair
        encodeIndex(j1, j2, pkBlocks[i].idx, nCtxts);
        bits[i * nCtxts + j1] = bits[i*nCtxts +j2] = 1; // set these two bits to one
    }
    //encrypt(ctxts, bits, nCtxts * subset_size);
}

#endif //GENTRYCRYPTO_HOMOMORPHIC_FHS_HPP
