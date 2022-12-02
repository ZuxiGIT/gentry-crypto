#ifndef GENTRYCRYPTO_HOMOMORPHIC_SWHS_HPP
#define GENTRYCRYPTO_HOMOMORPHIC_SWHS_HPP

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <string>
#include <bitset>
#include <time.h>

/*
 * Somewhat homomorphic scheme, which was built on the idea of Smart-Vercauteren.
 * Unlike the first, key generation in SWHS is done in O(n^{1.5}).
 */
class SWHS
{
    protected:
        unsigned int lattice_log_dim = 0;   // n
        unsigned int coeff_bitsize   = 0;   // t
        unsigned int prec_parameter  = 4;   // p
        double       noise           = 8;   // noise

        NTL::ZZ root;                       // r
        NTL::ZZ lattice_det;                // d
        NTL::ZZ private_key;                // w

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

    public:
        ~SWHS();
        explicit SWHS(unsigned int, unsigned int);

        SWHS(const SWHS&);
        SWHS& operator=(const SWHS&);

        SWHS(SWHS&&);
        SWHS& operator=(SWHS&&);

        void keyGen();
        bool encrypt(NTL::vec_ZZ&, unsigned int      ) const;
        bool encrypt(NTL::vec_ZZ&, unsigned int*, int) const;
        bool encrypt(NTL::vec_ZZ&, std::string       ) const;

        unsigned int decrypt(const NTL::ZZ&) const;
        int decrypt(unsigned int*, const NTL::vec_ZZ&) const;
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
int SWHS::randomBit(double p = 0.5) const
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

bool SWHS::aboveThreshold(long n, long m) const
{
    return (n + 1 + (m & 1) < m / 2);
}

void SWHS::evalRandPoly(NTL::vec_ZZ& vals, NTL::ZZ& r_to_m, long n, long m, double p, const NTL::ZZ& r, const NTL::ZZ& M) const
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

void SWHS::basicRandPoly(NTL::vec_ZZ& vals, NTL::ZZ& r_to_m, long n, long m, double p, const NTL::ZZ& r, const NTL::ZZ& M) const
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

std::string toBinary(const std::string& str)
{
    std::string binary = "";
    for (unsigned int i = 0; i < str.size(); ++i)
    {
        binary += std::bitset<8>(str.c_str()[i]).to_string();
    }
    return binary;
}

/*
std::string toString(NTL::vec_ZZ& b, int num)
{
    std::string str = "";
    for (int i = 0; i < num; i += 8)
    {
        str += static_cast<char>((1 << b[  i  ]) + (1 << b[i + 1]) +
                                 (1 << b[i + 2]) + (1 << b[i + 3]) +
                                 (1 << b[i + 4]) + (1 << b[i + 5]) +
                                 (1 << b[i + 6]) + (1 << b[i + 7]));
    }
    return str;
}*/

bool SWHS::encrypt(NTL::vec_ZZ& c, std::string text) const
{
    std::string btext = toBinary(text);
    unsigned int btext_length = btext.length();
    unsigned int b[btext_length];

    for (int i = 0; i < btext_length; i++)
    {
        b[i] = static_cast<int>(btext.c_str()[i]);
    }
    return encrypt(c, b, btext_length);
}

bool SWHS::encrypt(NTL::vec_ZZ& c, unsigned int b[], int num) const
{
    int i;
    NTL::ZZ tmp;
    unsigned long n = 1UL<<(lattice_log_dim);
    double p = (double) noise / n;
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

bool SWHS::encrypt(NTL::vec_ZZ& c, unsigned int m) const
{
    NTL::vec_ZZ vc(NTL::INIT_SIZE, 1);
    vc = c;
    if (!encrypt(vc, &m, 1)) return false;
    c = vc;
    return true;
}

unsigned int SWHS::decrypt(const NTL::ZZ& c) const
{
    NTL::ZZ e;
    MulMod(e, c, private_key, lattice_det);
    unsigned int lsb = bit(e, 0);
    if (NTL::IsOdd(lattice_det))
    {
        e <<= 1;
        if (e > lattice_det)
        {
            lsb ^= 1; // toggle the bit
        }
    }
    return lsb;
}

int SWHS::decrypt(unsigned int m[], const NTL::vec_ZZ& c) const
{
    int symbol = 0;
    for (symbol = 0; symbol < c.length(); symbol++)
    {
        if ((m[symbol] = decrypt(c[symbol])) == -1)
            break;
    }
    return symbol;
}

#endif //GENTRYCRYPTO_HOMOMORPHIC_SWHS_HPP
