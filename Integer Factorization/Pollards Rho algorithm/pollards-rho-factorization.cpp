#include <iostream>
#include <string>
#include <algorithm>
#include <map>
#include <random>
#include <chrono>
#include <cctype>
#include <cstdint>

using u128 = unsigned __int128;
using u64  = std::uint64_t;

// ---------- Input / Output for u128 ----------
std::istream& operator>>(std::istream& is, u128& val) {
    std::string s;
    is >> s;
    val = 0;
    for (char c : s) {
        if (std::isdigit((unsigned char)c)) val = val * 10 + (u128)(c - '0');
    }
    return is;
}

std::ostream& operator<<(std::ostream& os, u128 val) {
    if (val == 0) return os << "0";
    std::string s;
    while (val) {
        int digit = (int)(val % 10);
        s.push_back(char('0' + digit));
        val /= 10;
    }
    std::reverse(s.begin(), s.end());
    return os << s;
}

// ---------- small helpers ----------
static inline u128 gcd_u128(u128 a, u128 b) {
    while (b) {
        u128 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

// add modulo n without overflow
static inline u128 add_mod(u128 a, u128 b, u128 n) {
    // return (a + b) % n without overflowing u128
    if (a >= n - b) return a + b - n;
    return a + b;
}

// sub modulo n (result in [0,n))
static inline u128 sub_mod(u128 a, u128 b, u128 n) {
    return (a >= b) ? (a - b) : (n - (b - a));
}

// compute 2^k mod n via doubling (no overflow)
static inline u128 pow2_mod(u128 n, int k) {
    u128 x = 1 % n;
    for (int i = 0; i < k; ++i) {
        x = add_mod(x, x, n);
    }
    return x;
}

// ---------- RNG ----------
static std::mt19937_64 rng((u64)std::chrono::high_resolution_clock::now().time_since_epoch().count());

static inline u128 rand_u128() {
    u128 hi = (u128)rng();
    u128 lo = (u128)rng();
    return (hi << 64) | lo;
}
static inline u128 rand_range(u128 lo, u128 hi) { // inclusive
    u128 span = hi - lo + 1;
    return lo + (rand_u128() % span);
}

// ---------- Montgomery arithmetic for modulus < 2^128 (2 limbs of 64-bit) ----------
struct MontCtx {
    u128 n;
    u64 n0, n1;      // limbs of n: n = n0 + 2^64*n1
    u64 n0prime;     // -n0^{-1} mod 2^64
    u128 r2;         // R^2 mod n, where R = 2^128
    u128 oneM;       // 1 in Montgomery domain (== R mod n)
    u128 nm1M;       // (n-1) in Montgomery domain

    explicit MontCtx(u128 mod) : n(mod) {
        n0 = (u64)n;
        n1 = (u64)(n >> 64);

        // n must be odd for Montgomery inverse modulo 2^64 to exist
        // inv via Newton iteration: inv ≡ n0^{-1} (mod 2^64)
        u64 inv = 1;
        for (int i = 0; i < 6; ++i) {
            inv *= (u64)(2 - (u64)((u128)n0 * inv));
        }
        n0prime = (u64)(0 - inv);

        r2 = pow2_mod(n, 256);     // 2^256 mod n

        oneM = mont_mul(1 % n, r2);
        nm1M = mont_mul((n - 1) % n, r2);
    }

    // Montgomery multiplication: returns a*b*R^{-1} mod n, with R = 2^128
    inline u128 mont_mul(u128 a, u128 b) const {
        u64 a0 = (u64)a, a1 = (u64)(a >> 64);
        u64 b0 = (u64)b, b1 = (u64)(b >> 64);

        // 256-bit product t = a*b in 4x64-bit limbs (little-endian)
        u64 t0, t1, t2, t3;

        u128 p00 = (u128)a0 * b0;
        t0 = (u64)p00;
        u64 carry = (u64)(p00 >> 64);

        u128 p01 = (u128)a0 * b1;
        u128 p10 = (u128)a1 * b0;

        u128 s1 = (u128)(u64)p01 + (u128)(u64)p10 + carry;
        t1 = (u64)s1;
        carry = (u64)(s1 >> 64);

        u128 p11 = (u128)a1 * b1;

        u128 s2 = (u128)(u64)(p01 >> 64) + (u128)(u64)(p10 >> 64) + (u128)(u64)p11 + carry;
        t2 = (u64)s2;
        carry = (u64)(s2 >> 64);

        u128 s3 = (u128)(u64)(p11 >> 64) + carry;
        t3 = (u64)s3;

        // Montgomery reduction for k=2 limbs
        // Step 0
        u64 m0 = (u64)((u128)t0 * n0prime); // mod 2^64 automatically

        u128 u = (u128)m0 * n0 + t0;
        t0 = (u64)u;
        carry = (u64)(u >> 64);

        u = (u128)m0 * n1 + t1 + carry;
        t1 = (u64)u;
        carry = (u64)(u >> 64);

        u = (u128)t2 + carry;
        t2 = (u64)u;
        carry = (u64)(u >> 64);

        t3 += carry;
        
        // shift right by 64 (divide by base)
        t0 = t1; t1 = t2; t2 = t3; t3 = 0;

        // Step 1
        u64 m1 = (u64)((u128)t0 * n0prime);

        u = (u128)m1 * n0 + t0;
        t0 = (u64)u;
        carry = (u64)(u >> 64);

        u = (u128)m1 * n1 + t1 + carry;
        t1 = (u64)u;
        carry = (u64)(u >> 64);

        u = (u128)t2 + carry;
        t2 = (u64)u;
        carry = (u64)(u >> 64);

        // shift right by 64
        t0 = t1; t1 = t2;

        u128 res = ((u128)t1 << 64) | (u128)t0;

        if (res >= n) res -= n;
        return res;
    }

    inline u128 to_mont(u128 x) const {
        x %= n;
        return mont_mul(x, r2); // x * R^2 * R^{-1} = x*R (Montgomery)
    }

    inline u128 from_mont(u128 xM) const {
        return mont_mul(xM, 1);  // xM * 1 * R^{-1} = x (normal)
    }
};

static inline u128 mont_pow(u128 a, u128 e, const MontCtx& ctx) {
    u128 aM = ctx.to_mont(a);
    u128 rM = ctx.oneM;
    while (e) {
        if (e & 1) rM = ctx.mont_mul(rM, aM);
        aM = ctx.mont_mul(aM, aM);
        e >>= 1;
    }
    return rM; // still Montgomery domain
}

// ---------- Miller–Rabin prime test (probabilistic) ----------
static bool isPrime(u128 n) {
    if (n < 2) return false;

    static const u64 smallPrimes[] = {2ULL,3ULL,5ULL,7ULL,11ULL,13ULL,17ULL,19ULL,23ULL,29ULL,31ULL,37ULL};
    for (u64 p : smallPrimes) {
        if (n == (u128)p) return true;
        if (n % (u128)p == 0) return false;
    }
    if ((n & 1) == 0) return false;

    // n-1 = d * 2^s
    u128 d = n - 1;
    unsigned s = 0;
    while ((d & 1) == 0) { d >>= 1; ++s; }

    MontCtx ctx(n);

    auto trial = [&](u128 a) -> bool {
        a %= n;
        if (a == 0) return true;

        u128 xM = mont_pow(a, d, ctx);
        if (xM == ctx.oneM || xM == ctx.nm1M) return true;

        for (unsigned i = 1; i < s; ++i) {
            xM = ctx.mont_mul(xM, xM);
            if (xM == ctx.nm1M) return true;
        }
        return false;
    };

    for (u64 a : smallPrimes) if (!trial((u128)a)) return false;

    // extra random rounds (tune)
    for (int i = 0; i < 6; ++i) {
        u128 a = rand_range(2, n - 2);
        if (!trial(a)) return false;
    }
    return true;
}

// ---------- Pollard Rho (Brent version), Montgomery domain ----------
static u128 pollard_rho(u128 n) {
    if ((n & 1) == 0) return 2;
    if (n % 3 == 0) return 3;

    MontCtx ctx(n);
    const u128 m = 1024; // batch size

    while (true) {
        // choose random seed + c
        u128 y = rand_range(1, n - 1);
        u128 c = rand_range(1, n - 1);

        u128 yM = ctx.to_mont(y);
        u128 cM = ctx.to_mont(c);

        auto fM = [&](u128 xM) -> u128 {
            // x^2 + c in Montgomery domain
            u128 xxM = ctx.mont_mul(xM, xM);
            return add_mod(xxM, cM, n);
        };

        u128 g = 1, r = 1;
        u128 qM = ctx.oneM;
        u128 xM = 0, ysM = 0;

        while (g == 1) {
            xM = yM;
            for (u128 i = 0; i < r; ++i) yM = fM(yM);

            u128 k = 0;
            while (k < r && g == 1) {
                ysM = yM;
                u128 lim = (m < (r - k)) ? m : (r - k);

                for (u128 i = 0; i < lim; ++i) {
                    yM = fM(yM);
                    u128 diffM = sub_mod(xM, yM, n); // (x-y)R mod n
                    qM = ctx.mont_mul(qM, diffM);
                }

                g = gcd_u128(qM, n);
                k += lim;
            }
            r <<= 1;
        }

        if (g == n) {
            // fallback: single-step gcds
            do {
                ysM = fM(ysM);
                u128 diffM = sub_mod(xM, ysM, n);
                g = gcd_u128(diffM, n);
            } while (g == 1);
        }

        if (g != n) return g; // success; otherwise restart
    }
}

// ---------- factoring ----------
static void factor(u128 n, std::map<u128,int>& out) {
    if (n == 1) return;
    if (isPrime(n)) { out[n]++; return; }

    // small trial division
    for (u64 p = 2; p <= 10000; ++p) {
        if (n % (u128)p == 0) {
            int cnt = 0;
            while (n % (u128)p == 0) { n /= (u128)p; ++cnt; }
            out[(u128)p] += cnt;
            factor(n, out);
            return;
        }
    }

    u128 d = pollard_rho(n);
    factor(d, out);
    factor(n / d, out);
}

static void printFactors(const std::map<u128,int>& f) {
    bool first = true;
    for (auto& kv : f) {
        if (!first) std::cout << " * ";
        first = false;
        std::cout << kv.first;
        if (kv.second > 1) std::cout << "^" << kv.second;
    }
    std::cout << "\n";
}

int main() {
    u128 n;
    std::cin >> n;

    auto start = std::chrono::high_resolution_clock::now();

    std::map<u128,int> factors;
    factor(n, factors);

    auto end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - start).count();

    printFactors(factors);
    std::cout << "Time: " << ms << " ms\n";
    return 0;
}
