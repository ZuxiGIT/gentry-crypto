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
#include <NTL/RR.h>
#include <vector>
#include <stack>

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

// void AddMod(NTL::ZZ& x, const NTL::ZZ& a, long b, const NTL::ZZ& n)
// {
//    static NTL::ZZ B;
//    NTL::conv(B, b);
//    NTL::AddMod(x, a, B, n);
// }
//
// void AddMod(NTL::ZZ& x, const NTL::ZZ& a, const NTL::ZZ& b, const NTL::ZZ& n)
// {
//     NTL_zaddmod(a.rep, b.rep, n.rep, &x.rep);
// }

typedef std::vector<PKblock> PKblockv;
typedef std::stack<NTL::ZZ> ZZstack;

class FHS
{
    private:
        unsigned int lattice_log_dim  = 0;   // n
        unsigned int coeff_bitsize    = 0;   // t
        unsigned int sec_parameter    = 0;   // lambda
        unsigned int subset_size      = 0;   // s
        unsigned int prec_parameter   = 4;   // p
        unsigned int squash_poly_deg  = 0;   // d
        unsigned int bigset_size      = 0;   // S
        unsigned int bigset_log_ratio = 0;   // logR
        float        bdd_parameter    = 0;   // mu

        double       noise            = 256; // noise

        NTL::ZZ root;                        // r
        NTL::ZZ lattice_det;                 // d
        NTL::ZZ private_key;                 // w

        PKblockv     pkBlocks;
        NTL::vec_ZZ  ctxts;

        // Algorithms for key generation
        void polyCompress(NTL::ZZX&, unsigned int);
        void ratioCoeffs (NTL::ZZ&, NTL::ZZ&,           const NTL::ZZX&, unsigned int);
        int  invModFx    (NTL::ZZ&, NTL::ZZ&, NTL::ZZ&, const NTL::ZZX&, unsigned int);

        // Algorithms for encrypting
        bool aboveThreshold(long, long) const;
        void evalRandPoly (NTL::vec_ZZ&, NTL::ZZ&, long, long, double, const NTL::ZZ&, const NTL::ZZ&) const;
        void basicRandPoly(NTL::vec_ZZ&, NTL::ZZ&, long, long, double, const NTL::ZZ&, const NTL::ZZ&) const;

        // Algorithms for generating random numbers and random bits
        void random_ZZX  (NTL::ZZX&,    long, long);
        void random_vecZZ(NTL::vec_ZZ&, long, long);
        int  randomBit   (double) const;

        unsigned long mChoose2(unsigned long);
        void encodeIndex(unsigned long&, unsigned long&, unsigned long, unsigned long);
        //void gradeSchoolAdd(NTL::ZZ&, const NTL::mat_ZZ&, const NTL::ZZ&);

        unsigned long getBinaryRep(NTL::ZZ,      const NTL::ZZ&, long);
        void processBlock         (NTL::vec_ZZ&, const NTL::ZZ&, long) const;
        void evalSymPolys         (NTL::vec_ZZ&, ZZstack&, long, const NTL::ZZ&);
        bool verifyVector   (const NTL::vec_ZZ&, const NTL::ZZ&, const PKblock&);

    public:
        FHS() = delete;
        ~FHS();
        explicit FHS(unsigned int, unsigned int, unsigned int, unsigned int,
                     unsigned int, unsigned int, unsigned int);

        FHS(const FHS&);
        FHS& operator=(const FHS&);

        FHS(FHS&&);
        FHS& operator=(FHS&&);

        void keyGen();
        bool encrypt(NTL::vec_ZZ&, unsigned int      ) const;
        bool encrypt(NTL::vec_ZZ&, unsigned int*, int) const;
        void recrypt(NTL::ZZ&) const;
};

void FHS::random_vecZZ(NTL::vec_ZZ &vec, long dim, long bitsize)
{
    vec.SetLength(dim);
    for (unsigned int i = 0; i < vec.length(); i++)
    {
        NTL::RandomBits(vec[i], bitsize);
        vec[i] -= NTL::power2_ZZ(bitsize - 1);
    }
}
void FHS::random_ZZX(NTL::ZZX &poly, long dim, long bitsize)
{
    clear(poly);
    random_vecZZ(poly.rep, dim, bitsize);
    poly.normalize();
}

/*
 * A random bit is generated equal to 0 with some probability q,
 * and is generated equal to +-1 with probability (1 - q) / 2.
 */
int FHS::randomBit(double p = 0.5) const
{
    static unsigned long f = (unsigned long) - 1;
    unsigned long w = NTL::RandomWord();
    return (w <= p * f) ? ((w & 1) ? 1 : -1) : 0;
}

void FHS::polyCompress(NTL::ZZX& P, unsigned int init_coeff)
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
void FHS::ratioCoeffs(NTL::ZZ& g0, NTL::ZZ& g1, const NTL::ZZX& q, unsigned int n)
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
int FHS::invModFx(NTL::ZZ& d, NTL::ZZ& root, NTL::ZZ& wi, const NTL::ZZX& q, unsigned int n)
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

bool FHS::aboveThreshold(long n, long m) const
{
    return (n + 1 + (m & 1) < m / 2);
}

void FHS::evalRandPoly(NTL::vec_ZZ& vals, NTL::ZZ& r_to_m, long n, long m, double p, const NTL::ZZ& r, const NTL::ZZ& M) const
{
    if (aboveThreshold(n,m))
    {
        evalRandPoly(vals, r_to_m, 2*n, m/2, p, r, M);
        NTL::ZZ tmp;
        int s;
        for (int i=0; i<n; i++)
        {
            if ((m & 1) && (s = randomBit(p)))
            if (s==1) AddMod(vals[i+n],vals[i+n],r_to_m,M);
            else      SubMod(vals[i+n],vals[i+n],r_to_m,M);
            MulMod(tmp, vals[i + n], r_to_m, M);
            vals[i] += tmp;
        }
        for (int i=0; i<n; i++) vals[i+n].kill();
        vals.SetLength(n);

        SqrMod(r_to_m, r_to_m, M);
        if (m & 1) MulMod(r_to_m,r_to_m,r,M);
    }
    else basicRandPoly(vals, r_to_m, n, m, p, r, M);
}

void FHS::basicRandPoly(NTL::vec_ZZ& vals, NTL::ZZ& r_to_m, long n, long m, double p, const NTL::ZZ& r, const NTL::ZZ& M) const
{
    int i, j, k, s;
    vals.SetLength(n);
    if (m <= 0)
    {
        r_to_m = NTL::to_ZZ(1);
        return;
    }
    for (i=0; i<n; i++) vals[i] = randomBit(p); // the free term (0/1)
    if (m==1)
    {
        r_to_m = r;
        return;
    }
    NTL::ZZ rSqr;
    SqrMod(rSqr,r,M);
    for (i = 0; i < n; i++)
    {
        if (s=randomBit(p))
        if (s==1) vals[i] += r;
        else      SubMod(vals[i], vals[i], r, M);
        if (m>2 && (s=randomBit(p)))
        if (s==1) NTL::AddMod(vals[i], vals[i], rSqr, M);
        else      NTL::SubMod(vals[i], vals[i], rSqr, M);
    }
    if (m>4) {
        r_to_m = rSqr;
        for (j=4; j<m; j*=2)
	{
            SqrMod(r_to_m, r_to_m, M); // r^j := (previous-r^j)^2
            for (i=0; i<n; i++)        // vals[i] += b_{i,j} * r^j mod M
            if (s=randomBit(p))
            if (s==1) NTL::AddMod(vals[i],vals[i],r_to_m,M);
            else      NTL::SubMod(vals[i],vals[i],r_to_m,M);
        }
    }
    else if (m<4) { // if m==2 or 3 we're done, just return the correct r_to_m
        if (m==2) r_to_m = rSqr;
        else      NTL::MulMod(r_to_m,rSqr,r,M);
        return;
    }

    NTL::ZZ r_odd_pwr = r;
    for (int j = 3; j < m; j += 2)
    {
        MulMod(r_odd_pwr, r_odd_pwr, rSqr, M); // next odd power of r
        r_to_m = r_odd_pwr;
        k = j;
        while (true) {
        for (i=0; i<n; i++)     // vals[i] += b_{i,j} * r^k mod M
        if (s=randomBit(p))
        if (s==1) NTL::AddMod(vals[i],vals[i],r_to_m,M);
        else      NTL::SubMod(vals[i],vals[i],r_to_m,M);
        k *= 2;
        if (k >= m) break;
        NTL::SqrMod(r_to_m, r_to_m, M); // r^k := (previous-r^k)^2 mod M
        }
    }

    // r_odd_power is r^{m-1} or r^{m-2}, depending  on whether m is even or odd
    if (m&1) NTL::MulMod(r_to_m, r_odd_pwr, rSqr, M);
    else     NTL::MulMod(r_to_m, r_odd_pwr, r, M);
}

unsigned long FHS::mChoose2(unsigned long S)
{
    return static_cast<unsigned long>(ceil(2 * sqrt(static_cast<double>(S))));
}

void FHS::encodeIndex(unsigned long& j1, unsigned long& j2, unsigned long  i,  unsigned long  m)
{
    if (i > m * (m - 1) / 2)
    {
        std::cerr << "encodeIndex: index " << i << " can't be encoded in (" << m << " choose 2)\n";
        abort();
    }
    for (j1 = 0; j1 < m - 1; j1++) {
        unsigned long pairsFor_j1 = m - j1 - 1;
        if (pairsFor_j1 > i)
        {
            j2 = i + j1 + 1;
            return;
        }
        else i -= pairsFor_j1; // we already counted the pairs for this j1
    }
    NTL::Error("encodeIndex: something went wrong, maybe an overflow?");
}

FHS::FHS(unsigned int n, unsigned int t, unsigned int s, unsigned int p,
         unsigned int d, unsigned int S, unsigned int R)

: lattice_log_dim(n),
  coeff_bitsize(t),
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
  coeff_bitsize(rhs.coeff_bitsize),
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
  prec_parameter(std::move(rhs.prec_parameter)),
  squash_poly_deg(std::move(rhs.squash_poly_deg)),
  bigset_size(std::move(rhs.bigset_size)),
  bigset_log_ratio(std::move(rhs.bigset_log_ratio))
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
    encrypt(ctxts, bits, nCtxts * subset_size);
}

bool FHS::encrypt(NTL::vec_ZZ& c, unsigned int b[], int num) const
{
    int i;
    NTL::ZZ tmp;
    unsigned long n = 1UL<<(lattice_log_dim);
    double p = noise / n;
    if (p > 0.5)
    {
        p = 0.5;
    }

    evalRandPoly(c, tmp, num, lattice_log_dim, prec_parameter, root, lattice_det);
    for (i=0; i < num; i++)
    {
        c[i] <<= 1;
        c[i] += b[i];
        if (c[i] >= lattice_det)
        {
            c[i] -= lattice_det;
        }
    }
    return true;
}

bool FHS::encrypt(NTL::vec_ZZ& c, unsigned int m) const
{
    NTL::vec_ZZ vc(NTL::INIT_SIZE, 1);
    vc = c;
    if (!encrypt(vc, &m, 1)) return false;
    c = vc;
    return true;
}

void FHS::recrypt(NTL::ZZ& c) const
{
    clear(fhe_recrypt_sum);
    if (testCtxt(c) < prec_parameter + 1) std::cerr << "recrypt: c * w larger than det / 2^{p+1}!!\n";

    mat_ZZ vars(NTL::INIT_SIZE, subset_size, prec_parameter +1);
    for (long i=0; i<subset_size; i++) processBlock(vars[i], c, i);

    NTL::ZZ tmp;
    NTL::MulMod(tmp, c, w, det);
    if (tmp != fhe_recrypt_sum)
        NTL::Error("Recrypt: sum from blocks differ from w*c");

    gradeSchoolAdd(c, vars, det);
    timer += GetTime();
    std::cerr << "Grade-school addition algorithm in time "<<timer<<std::endl;
}

unsigned long FHS::getBinaryRep(NTL::ZZ n, const NTL::ZZ& d, long nBits)
{
    NTL::RR ratio = to_RR(n)/to_RR(d);
    n <<= nBits;
    n /= d;         // integer division implies truncation

    unsigned long sn = to_long(n);
    sn = (sn >> 1) + (sn & 1);
    long twoToP = 1L << (nBits - 1);
    double approx = ((double)sn) / ((double)twoToP);

    NTL::RR error = abs(ratio-approx);  // check that |ratio-approx| is small enough
    if (error*(twoToP * 2) > 1)
    {
        NTL::Error("getBinaryRep: sn/2^p too far from n/d");
    }
    return sn;
}

void FHS::processBlock(NTL::vec_ZZ& vars, const NTL::ZZ& c, long i) const
{
    double ptimer = 0.0;
    double atimer = 0.0;
    double mtimer = 0.0;

    long nCtxts = mChoose2(subset_size);
    unsigned long baseIdx = i * nCtxts;

    int k;
    for (k = 0; k < vars.length(); k++) clear(vars[k]); // initialize to zero

    unsigned long j, j1, j2;
    NTL::ZZ factor = pkBlocks[i].x;
    NTL::MulMod(factor, factor, c, det);
    NTL::vec_ZZ psums(NTL::INIT_SIZE, vars.length()); // keep partial sums
    for (j = j1 = 0; j1 < nCtxts-1; j1++) {       // sk-bits indexed by (j1,*) pairs
        for (k = 0; k < psums.length(); k++) clear(psums[k]);  // initialize to zero

        for (j2 = j1+1; j2 < nCtxts; j2++)
        {
            unsigned long binary = NTL::getBinaryRep(factor, det, vars.length());
            if (IsOdd(factor))     // "xor" the LSB to column 0
            binary ^= (1UL << prec_parameter);

            for (k=0; k<psums.length(); k++) if (NTL::bit(binary, k) == 1) {
            long k2 = psums.length() -k-1;
            AddMod(psums[k2], psums[k2], ctxts[baseIdx+j2], det);
        }

        j++;              // done with this element
        if (j < subset_size) { // compute next element = current * R mod det

        factor <<= bigset_log_ratio;
        factor %= det;
        }
            else break;
        }

        for (k = 0; k < vars.length(); k++)
        {
            NTL::MulMod(psums[k], psums[k], ctxts[baseIdx+j1], det);
            NTL::AddMod(vars[k],vars[k],psums[k], det);
        }
        if (j >= subset_size) break;
    }
    // Sanity-check: j should be at least S, else we've missed some terms
    if (j < subset_size) NTL::Error("FHEkeys::processBlock: loop finished with j<S");

    if (!verifyVector(vars, c, pkBlocks[i]))
        NTL::Error("FHEkeys::processBlock: decrypted vector does not match");
}

void FHS::gradeSchoolAdd(NTL::ZZ& c, const mat_ZZ& vars, const ZZ& M)
{
    long i,j;
    NTL::vec_long zCols(NTL::INIT_SIZE, vars.NumCols());
    for (j=0; j<vars.NumCols(); j++) for (i=0; i<vars.NumRows(); i++)
        zCols[j] += fhe_recrypt_curKey->decrypt(vars[i][j]);
    std::cerr << "  Column-weight in input to gradeSchoolAdd: " << zCols << std::endl;

    // Below it is more convenient to have each column of the matrix in
    // a separate stack (since we would want to push carry bits on top)
    std::vector<ZZstack> columns(vars.NumCols());
    for (j=0; j<vars.NumCols(); j++) for (i=vars.NumRows()-1; i>=0; i--)
        columns[j].push(vars[i][j]);

    NTL::vec_ZZ sp; // space to store symmetric polynomials

    // add columns from right to left, upto column -1
    for (j=vars.NumCols()-1; j>0; j--)
    {
        long s = columns[j].size();
        long log = NTL::NextPowerOfTwo(s); // (log of) # of carry bits to compute
        if (log > j)       log = j;     // no more carry than what can reach col 0
        if ((1L<<log) > s) log--; // no more carry than what s bits can produce

        evalSymPolys(sp, columns[j], 1L << log, M);
        zCols[j] = fhe_recrypt_curKey->decrypt(sp[1]);
        long k = 2;
        for (long j2 = j - 1; j2 >= 0 && k < sp.length(); j2--)
        {
            columns[j2].push(sp[k]);   // push carry bits on top of their column
            k <<= 1;
        }
    }

    c = columns[0].top();
    columns[0].pop();
    while (!columns[0].empty())
    {
        NTL::AddMod(c, c, columns[0].top(), M);
        columns[0].pop();
    }
    zCols[0] = fhe_recrypt_curKey->decrypt(c);
    std::cerr << " output from gradeSchoolAdd: " << zCols << std::endl;
    NTL::AddMod(c, c, sp[1], M);
}

void FHS::evalSymPolys(NTL::vec_ZZ& out, ZZstack& vars, long deg, const NTL::ZZ& M)
{
    long i, j;
    out.SetLength(deg + 1);
    set(out[0]);
    for (i=1; i<=deg; i++) clear(out[i]);

    NTL::ZZ tmp;
    for (i=1; !vars.empty(); i++)
    {
        for (j=NTL::min(i,deg); j>0; j--)
        {
            NTL::MulMod(tmp, out[j-1], vars.top(), M);
            NTL::AddMod(out[j], out[j], tmp, M);
        }
        vars.pop();
    }
}

bool FHS::verifyVector(const NTL::vec_ZZ& vars, const NTL::ZZ& c, const PKblock& block)
{
    unsigned int v[vars.length()];
    fhe_recrypt_curKey->decrypt(v, vars);

    NTL::ZZ factor;
    long logR = fhe_recrypt_curKey->getPrms().logR;
    long p = fhe_recrypt_curKey->getPrms().p;
    const NTL::ZZ& det = fhe_recrypt_curKey->getDet();

    NTL::power2(factor, block.idx * logR);     // 2^{idx*logR} = R^{idx}
    if (factor>det) factor %= det;
    NTL::MulMod(factor, block.x, factor, det); // factor *= x mod det
    NTL::MulMod(factor, c, factor, det);       // factor *= c mod det

    NTL::AddMod(fhe_recrypt_sum, fhe_recrypt_sum, factor, det);
    unsigned long binary = getBinaryRep(factor, det, p+1);
    if (IsOdd(factor)) v[0] ^= 1; // xor the LSB

    for (long j=0; j<=p; j++) if (NTL::bit(binary,j) != v[p-j])
    {
        std::cerr << "binary = "<<binary<<" but v[0-4]= ["
        <<v[0]<<","<<v[1]<<","<<v[2]<<","<<v[3]<<","<<v[4]<<"]\n";
        return false;
    }
    std::cerr << "verifyVector: binary = "<<binary<<" verified\n";
    return true;
}

#endif //GENTRYCRYPTO_HOMOMORPHIC_FHS_HPP
