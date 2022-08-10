/*
 * Description: Modulo integer class, usually only fastmod is useful
 */
#include <bits/stdc++.h>

using namespace std;
using ll = int64_t;

struct Z {
    uint32_t v;
    using u64 = uint64_t;
    Z() = default;
    Z(const Z &) = default;
    static inline u64 p, ip;
    static void setMod(int m) { p = m, ip = -1ull / p; }
    // Useful if and only if P is constant but not known at compile time
    static u64 mod(u64 a) {
        u64 r = a - (u64)(((__uint128_t)ip * a) >> 64) * p;
        return r >= p ? r - p : r;
    }
    template <class T = int> Z(T v = 0) : v((v %= int(p)) < 0 ? p + v : v) {}
    template <class T = int> explicit operator T() const { return T(v); }
    Z inv() const {
        int x = 1, y = 0, a = v, b = p, q, tx, ta;
        for (; b; a = b, b = ta - b * q)
            q = a / b, tx = x, ta = a, x = y, y = tx - y * q;
        return x;
    }
    Z &operator-() { return v && (v = p - v), *this; }
    Z &operator-=(const Z &b) { return v -= b.v, v += (p & -(v >> 31)), *this; }
    Z &operator+=(const Z &b) { return operator-=(p - b.v); }
    Z &operator*=(const Z &b) { return v = mod((u64) v * b.v), *this; }
    Z &operator/=(const Z &b) { return operator*=(b.inv()); }
    template <class T> bool operator==(const T &b) const { return T(v) == b; }
    template <class T> bool operator!=(const T &b) const { return T(v) != b; }
    friend Z operator+(const Z &a, const Z &b) { return Z(a) += b; }
    friend Z operator-(const Z &a, const Z &b) { return Z(a) -= b; }
    friend Z operator*(const Z &a, const Z &b) { return Z(a) *= b; }
    friend Z operator/(const Z &a, const Z &b) { return Z(a) /= b; }
};


// correctness test
const int P = 1e9 + 7;
int POW(ll a, int b = P - 2, ll x = 1) {
    for (; b >= 1; b >>= 1, a = a * a % P)
        if (b & 1) x = x * a % P;
    return x;
}
mt19937 rd(time(nullptr));
int range(int a, int b) { return uniform_int_distribution<int>(a, b)(rd); }
int main() {
    Z::setMod(P);
    int n = P - 1, m = 1e5;
    while (m--) {
        ll i = range(-n, n), j = range(-n, n), y;
        Z a = i, b = j;
        y = ((i + j) % P + P) % P;
        assert(a + b == y);
        y = ((i - j) % P + P) % P;
        assert(a - b == y);
        y = ((i * j) % P + P) % P;
        assert(a * b == y);
        y = ((i * POW(j % P + P)) % P + P) % P;
        assert(a / b == y);
    }
    return 0;
}