#ifndef GENTRYCRYPTO_HOMOMORPHIC_SWHS_HPP
#define GENTRYCRYPTO_HOMOMORPHIC_SWHS_HPP

/*
 * This project is made for the learning purpose of studying the
 * Smart-Vercauteren and Gentry-Halevi schemes for the subject of
 * Information Security. There is no purpose for commercial use.
 *
 * Code is based Craig Gentry's IBM implementation and the article
 * "Implementing a fully homomorphic Gentry-Halevi encryption scheme"
 * by Gentry-Halevi, IBM Research February 4, 2011
 */

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>

/*
 * Somewhat homomorphic scheme, which was built on the idea of Smart-Vercauteren.
 * Unlike the first, key generation in SWHS is done in O(n^{1.5}).
 */
class SWHS
{
    protected:
        unsigned int lattice_log_dim = 0;   // n
        unsigned int coeff_bitsize   = 0;   // t

        NTL::ZZ root;                       // r
        NTL::ZZ lattice_det;                // d
        NTL::ZZ private_key;                // w

        // Algorithms for key generation, the procedure for key generation is described above
        void polyCompress(NTL::ZZX&, unsigned int);
        void ratioCoeffs (NTL::ZZ&, NTL::ZZ&,           const NTL::ZZX&, unsigned int);
        int  invModFx    (NTL::ZZ&, NTL::ZZ&, NTL::ZZ&, const NTL::ZZX&, unsigned int);

        /*
         * Algorithms for generating random numbers and random bits
         */
        void random_ZZX  (NTL::ZZX&,    long, long);
        void random_vecZZ(NTL::vec_ZZ&, long, long);
        int  randomBit   (double);

    public:
        SWHS() = delete;
        ~SWHS();
        explicit SWHS(unsigned int, unsigned int);

        SWHS(const SWHS&);
        SWHS& operator=(const SWHS&);

        SWHS(SWHS&&);
        SWHS& operator=(SWHS&&);

        void keyGen();
};

void SWHS::random_vecZZ(NTL::vec_ZZ &vec, long dim, long bitsize)
{
    vec.SetLength(dim);
    for (unsigned int i = 0; i < vec.length(); i++)
    {
        NTL::RandomBits(vec[i], bitsize);
        vec[i] -= NTL::power2_ZZ(bitsize - 1);
    }
}
void SWHS::random_ZZX(NTL::ZZX &poly, long dim, long bitsize)
{
    clear(poly);
    random_vecZZ(poly.rep, dim, bitsize);
    poly.normalize();
}

/*
 * A random bit is generated equal to 0 with some probability q,
 * and is generated equal to +-1 with probability (1 - q) / 2.
 */
int SWHS::randomBit(double p = 0.5)
{
    static unsigned long f = (unsigned long) - 1;
    unsigned long w = NTL::RandomWord();
    return (w <= p * f) ? ((w & 1) ? 1 : -1) : 0;
}

void SWHS::polyCompress(NTL::ZZX& P, unsigned int init_coeff)
{
    unsigned int i = init_coeff;
    for (; i <= deg(P) / 2; i++)
    {
        P.rep[i] = P.rep[2 * i];
    }
    for (; i <= deg(P); i++)
    {
        NTL::clear(P.rep[i]);
    }
    P.normalize();
}

/*
 * In fact, this function finds the first two coefficients
 * g_0 = res(v, f_n) = d, g_1 = n * w_0 of the function g(z)
 * (defined over a complex field) by the modulo z^2.
 */
void SWHS::ratioCoeffs(NTL::ZZ& g0, NTL::ZZ& g1, const NTL::ZZX& q, unsigned int n)
{
    unsigned int N = 1 << n;

    //V = q, V2 = V(-x)
    NTL::ZZX V(q);
    NTL::ZZX V2(NTL::INIT_SIZE, N);

    //U = 1
    NTL::ZZX U(NTL::INIT_SIZE, N);
    SetCoeff(U, 0);

    //F = x^N + 1
    NTL::ZZX F(NTL::INIT_SIZE, N+1);
    SetCoeff(F, 0);
    SetCoeff(F, N);


    while (N > 1)
    {
        //Inverting odd degrees of a polynomial
        V2 = V;

        for (unsigned int i = 1; i <= deg(V2); i += 2)
        {
            NTL::negate(V2.rep[i], V2.rep[i]);
        }
        V2.normalize();

        //V = V(x) * V(-x) mod F = V * V2 mod F
        //U = V(x) * U(x)  mod F
        V = MulMod(V, V2, F);
        U = MulMod(U, V2, F);
        for (unsigned int i = 1; i <= deg(V); i += 2)
        {
            if (!NTL::IsZero(V.rep[i]))
            {
                 NTL::Error("[Error] One of the odd coeffs V is zero");
            }
        }

        //Polynomial "compression"
        polyCompress(V, 1);
        polyCompress(U, 0);

        /*
         * Obviously we are using computer tools, so we
         * do a bit shift of the variable N, in fact the
         * action takes place on the variable nj = n / 2^j
         */
        SetCoeff(F, N, 0);
        N >>= 1;
        SetCoeff(F, N);
        F.normalize();
    }

    g0 = ConstTerm(V);
    g1 = ConstTerm(U);
}

/*
 * In fact, this function creates a private key, ratio
 * r = w_1 / w_0 mod d and makes a lattice determinant check
 */
int SWHS::invModFx(NTL::ZZ& d, NTL::ZZ& root, NTL::ZZ& wi, const NTL::ZZX& q, unsigned int n)
{
    unsigned int N = 1 << n;
    unsigned int ind = 0;

    NTL::ZZ w0, w1;

    // Calculation of the coefficient w_0
    ratioCoeffs(d, w0, q, n);

    /*
     * From the Gentry-Halevi article, it is sufficient
     * for the HNF determinant to be odd, not prime.
     */
    if (!IsOdd(d))
    {
        return 0;
    }

    {
        // Calculation of the coefficient w_1
        NTL::ZZX qx(NTL::INIT_SIZE, N);
        for (unsigned int i = N - 1; i > 0; i--)
        {
            SetCoeff(qx, i, coeff(q, i - 1));
        }
        NTL::negate(qx.rep[0], coeff(q, N - 1));
        qx.normalize();
        ratioCoeffs(d, w1, qx, n);
    }

    //Mathematical operations for positive sign of numbers
    if (sign(d) == -1)
    {
        NTL::negate(d, d);
        NTL::negate(w0, w0);
        NTL::negate(w1, w1);
    }
    if (sign(w0)<0)
    {
        w0 += d;
    }
    if (sign(w1)<0)
    {
        w1 += d;
    }

    {
        NTL::ZZ tmp;

        /*
         * Check for the possibility of inverting w_0 mod d,
         * because r = w_1 / w_0 mod d = w_1 * w_0^{-1} mod d
         */

        if (InvModStatus(tmp, w0, d) != 0)
        {
            return 0;
        }

        //Check for a solution r^N = -1 over the set of complex numbers.
        root = MulMod(tmp, w1, d);
        tmp  = PowerMod(root, N, d);
        if (tmp + 1 != d)
        {
            return 0;
        }
    }

    // Check for odd values and their ranges
    if (((w0 <= d / 2) && IsOdd(w0)) || ((w0 > d / 2) && !IsOdd(w0)))
    {
        wi = w0;
    }
    else if (((w1 <= d / 2) && IsOdd(w1)) || ((w1 > d / 2) && !IsOdd(w1)))
    {
        wi = w1;
    }
    else
    {
        //Find the first odd w_ind to assign it private key status
        for (ind = 2; ind < N; ind++)
        {
            MulMod(w1, w1, root, d);
            if (((w1 <= d / 2) && IsOdd(w1)) || ((w1 > d / 2) && !IsOdd(w1)))
            {
                wi = w1;
                break;
            }
        }
    }
    return ((ind == N) ? 0 : 1);
}

SWHS::SWHS(unsigned int n, unsigned int t)
: lattice_log_dim(n), coeff_bitsize(t)
{}

SWHS::~SWHS()
{}

SWHS::SWHS(const SWHS& rhs)
: lattice_log_dim(rhs.lattice_log_dim),
  coeff_bitsize(rhs.coeff_bitsize)
{}

SWHS& SWHS::operator=(const SWHS& rhs)
{
    if (this != &rhs)
    {
        lattice_log_dim = rhs.lattice_log_dim;
        coeff_bitsize   = rhs.coeff_bitsize;
    }
    return *this;
}

SWHS::SWHS(SWHS&& rhs)
: lattice_log_dim(std::move(rhs.lattice_log_dim)),
  coeff_bitsize(std::move(rhs.coeff_bitsize))
{}

SWHS& SWHS::operator=(SWHS&& rhs)
{
    lattice_log_dim = std::move(rhs.lattice_log_dim);
    coeff_bitsize = std::move(rhs.coeff_bitsize);
    return *this;
}

void SWHS::keyGen()
{
    unsigned long bitsize = coeff_bitsize;
    unsigned long n = 1 << lattice_log_dim;
    long LSB = 0;

    NTL::ZZX v(NTL::INIT_SIZE, n);
    while (!invModFx(lattice_det, root, private_key, v, lattice_log_dim))
    {
        //Polynomial generation with odd coefficients
        random_ZZX(v, n, bitsize);
        LSB = NTL::trunc_long(coeff(v, 0), 1);

        for (int i = 1; i <= deg(v); i++)
        {
            LSB ^= NTL::trunc_long(coeff(v, i), 1);
        }
        if (LSB & 1 == 0)
        {
            v.rep[1]++;
        }
    }
}

#endif //GENTRYCRYPTO_HOMOMORPHIC_SWHS_HPP
