#include <bits/stdc++.h>

using namespace std;
const int N = 1 << 18, P = 998244353;
using ll = int64_t;

#define inc(a, b) (((a) += (b)) >= P ? (a) -= P : 0)
#define dec(a, b) (((a) -= (b)) < 0 ? (a) += P : 0)
#define mul(a, b) (ll(a) * (b) % P)
int POW(ll a, int b = P - 2, ll x = 1) {
    for (; b; b >>= 1, a = a * a % P)
        if (b & 1) x = x * a % P;
    return x;
}

int inv[N], fac[N], ifac[N], _ = [] {
    fac[0] = fac[1] = ifac[0] = ifac[1] = inv[1] = 1;
    for (ll i = 2; i < N; ++i) {
        fac[i] = fac[i - 1] * i % P;
        inv[i] = (P - P / i) * inv[P % i] % P;
        ifac[i] = ifac[i - 1] * inv[i] % P;
    }
    return 0;
}();

namespace NTT {
const int G = 3, L = 1 << 18;
int W[L], _ = [] {
    W[L / 2] = 1;
    for (int i = L / 2 + 1, wn = POW(G, P / L); i < L; ++i) W[i] = mul(W[i - 1], wn);
    for (int i = L / 2 - 1; ~i; --i) W[i] = W[i << 1];
    return 0;
}();
void dft(int *a, int n) {
    for (int k = n >> 1; k; k >>= 1)
        for (int i = 0; i < n; i += k << 1)
            for (int j = 0; j < k; ++j) {
                int &x = a[i + j], y = a[i + j + k];
                a[i + j + k] = mul(x - y + P, W[k + j]);
                inc(x, y);
            }
}
void idft(int *a, int n) {
    for (int k = 1; k < n; k <<= 1)
        for (int i = 0; i < n; i += k << 1)
            for (int j = 0; j < k; ++j) {
                int x = a[i + j], y = mul(a[i + j + k], W[k + j]);
                a[i + j + k] = x < y ? x - y + P : x - y;
                inc(a[i + j], y);
            }
    for (int i = 0, in = P - (P - 1) / n; i < n; ++i)
        a[i] = mul(a[i], in);
    reverse(a + 1, a + n);
}
} // namespace NTT

int norm(int n) { return 1 << (__lg(n - 1) + 1); }

struct Poly : public vector<int> {
#define T (*this)
    using vector<int>::vector;
    int deg() const { return size(); }
    Poly rev() const { return Poly(rbegin(), rend()); }
    void append(const Poly &a) { insert(end(), a.begin(), a.end()); }

    Poly operator-() const {
        Poly a(T);
        for (auto &x : a) x = x ? P - x : 0;
        return a;
    }
    Poly &operator+=(const Poly &a) {
        if (a.deg() > deg()) resize(a.deg());
        for (int i = 0; i < a.deg(); ++i) inc(T[i], a[i]);
        return T;
    }
    Poly &operator-=(const Poly &a) {
        if (a.deg() > deg()) resize(a.deg());
        for (int i = 0; i < a.deg(); ++i) dec(T[i], a[i]);
        return T;
    }
    Poly &operator^=(const Poly &a) {
        if (a.deg() < deg()) resize(a.deg());
        for (int i = 0; i < deg(); ++i) T[i] = mul(T[i], a[i]);
        return T;
    }
    Poly &operator*=(int a) {
        for (auto &x : T) x = mul(x, a);
        return T;
    }

    Poly operator+(const Poly &a) const { return Poly(T) += a; }
    Poly operator-(const Poly &a) const { return Poly(T) -= a; }
    Poly operator^(const Poly &a) const { return Poly(T) ^= a; }
    Poly operator*(int a) const { return Poly(T) *= a; }

    Poly &operator<<=(int k) { return insert(begin(), k, 0), T; }
    Poly operator<<(int r) const { return Poly(T) <<= r; }
    Poly operator>>(int r) const { return r >= deg() ? Poly() : Poly(begin() + r, end()); }
    Poly &operator>>=(int r) { return T = T >> r; }

    Poly pre(int k) const { return k < deg() ? Poly(begin(), begin() + k) : T; }
    friend void dft(Poly &a) { NTT::dft(a.data(), a.size()); }
    friend void idft(Poly &a) { NTT::idft(a.data(), a.size()); }
    friend Poly conv(Poly a, Poly b, int n) {
        a.resize(n), dft(a);
        b.resize(n), dft(b);
        return idft(a ^= b), a;
    }
    Poly operator*(const Poly &a) const {
        int n = deg() + a.deg() - 1;
        return conv(T, a, norm(n)).pre(n);
    }
    // f_i = sum a_{i + j} * T_j
    Poly mulT(const Poly& a) const { return T * a.rev() >> (a.deg() - 1); }

    // solve f_1, ..., f_{n - 1} with f_i = sum f_{i - j} * T_j, f_0 is known
    // O(n log^2 n)
    template <typename func>
    Poly semiConv(Poly f, int n, func calc) {
        f.resize(n);
        vector<Poly> va(n); // storage dft result, faster
        for (int m = 1; m < n; m++) {
            int k = m & -m, l = m - k, r = min(m + k, n);
            Poly &p = va[k], q(f.begin() + l, f.begin() + m);
            if (p.empty()) p = pre(r - l - 1), p.resize(k * 2), dft(p);
            q.resize(k << 1), dft(q), idft(q ^= p);
            // Poly q = conv(Poly(f.begin() + l, f.begin() + m), pre(r - l - 1), k << 1);
            for (int i = m; i < r; i++) inc(f[i], q[i - l - 1]);
            calc(f, m);
        }
        return f;
    }

    Poly inv(int n) const {
        Poly x{POW(T[0])};
        for (int k = 1; k < n; k <<= 1)
            x.append(-((conv(pre(k << 1), x, k << 1) >> k) * x).pre(k));
        return x.pre(n);
    }
    Poly deriv() const {
        if (empty()) return {};
        Poly a(deg() - 1);
        for (ll i = 1; i < deg(); ++i) a[i - 1] = i * T[i] % P;
        return a;
    }
    Poly integ() const {
        if (empty()) return {};
        Poly a(deg() + 1);
        for (int i = 1; i <= deg(); ++i) a[i] = mul(::inv[i], T[i - 1]);
        return a;
    }
    Poly log(int n) const { return (deriv() * inv(n)).integ().pre(n); }
    Poly exp(int n) const {
        Poly x{1};
        for (int k = 1; k < n; k <<= 1)
            x.append((x * ((pre(k << 1) - x.log(k << 1)) >> k)).pre(k));
        return x.pre(n);
    }
    Poly exp2(int n) const { return deriv().semiConv({1}, n, [&](Poly &f, int m) { f[m] = mul(f[m], ::inv[m]); }); }
    Poly sqrt(int n) const {
        Poly x{1}, y{1}; // x[0] = sqrt(T[0]), default T[0] = 1
        for (int k = 1; k < n; k <<= 1) {
            x.append((((pre(k << 1) - x * x) >> k) * y).pre(k) * ((P + 1) >> 1));
            k << 1 < n ? y.append(-((conv(x.pre(k << 1), y, k << 1) >> k) * y).pre(k)) : void();
        }
        return x.pre(n);
    }
    Poly pow(int k, int kp, int n) const { // k = K % P, kp = K % phi(P)
        int i = 0;
        while (i < deg() && !T[i]) i++;
        if (1ll * i * k >= n) return Poly(n);
        int v = T[i], m = n - i * k;
        return ((((T >> i) * POW(v)).log(m) * k).exp(m) << (i * k)) * POW(v, kp);
    }
    vector<Poly> operator/(const Poly &a) const {
        int k = deg() - a.deg() + 1;
        if (k < 0) return {{0}, T};
        Poly q = (rev().pre(k) * (a.rev().inv(k))).pre(k).rev(), r = T - a * q;
        return {q, r.pre(a.deg() - 1)};
    }

    // [x ^ k](f / g)
    int divAt(Poly f, Poly g, ll k) {
        int n = max(f.deg(), g.deg()), m = norm(n);
        for (; k; k >>= 1) {
            f.resize(m * 2), dft(f);
            g.resize(m * 2), dft(g);
            for (int i = 0; i < 2 * m; ++i) f[i] = mul(f[i], g[i ^ 1]);
            for (int i = 0; i < m; ++i) g[i] = mul(g[2 * i], g[2 * i + 1]);
            g.resize(m), idft(f), idft(g);
            for (int i = 0, j = k & 1; i < n; i++, j += 2) f[i] = f[j];
            f.resize(n), g.resize(n);
        }
        return f[0];
    }
    // T[k] where T[n] = sum c[i] * T[n - i]
    int recur(Poly c, ll k) { return c[0] = P - 1, divAt((T * c).pre(c.deg() - 1), c, k); }

    Poly eval(Poly x) const {
        if (empty()) return Poly(x.deg());
        const int m = x.deg(), n = norm(m);
        vector<Poly> q(2 * n, {1});
        Poly ans(m), temp, d(2 * n);
        for (int i = n; i < n + m; ++i) q[i] = Poly{1, P - x[i - n]};
        for (int i = n - 1; i; --i) q[i] = q[i << 1] * q[i << 1 | 1];
        q[1] = mulT(q[1].inv(n));
        for (int i = 1, l = 2, r = 3; i < n; ++i, l += 2, r += 2) {
            temp = q[l], d[l] = d[r] = d[i] + 1;
            q[l] = q[i].mulT(q[r]).pre(n >> d[l]);
            q[r] = q[i].mulT(temp).pre(n >> d[r]);
        }
        for (int i = n; i < n + m; ++i) ans[i - n] = q[i][0];
        return ans;
    }

#undef T
};

// useless algorithm
// G(F(x))
Poly compound(Poly f, Poly g) { // (5 * sqrt(n) + 3) * M(n) + n ^ 2
    int n = f.size(), k = norm(2 * n), L = sqrt(n) + 1;
    vector<Poly> G(L + 1);
    Poly H, h(k), t; // H = g ^ (iL)
    auto dft = [&](Poly &a) { a.resize(k), NTT::dft(a.data(), k); };
    auto idft = [&](Poly &a) { NTT::idft(a.data(), k), a.resize(n); };
    G[0] = H = {1}, dft(g);
    for (int i = 1; i <= L; ++i) dft(G[i] = G[i - 1]), idft(G[i] ^= g);
    dft(g = G[L]);
    for (int i = 0; i < L; ++i, t = {}, idft(H)) {
        for (int j = 0; j < min(L, n - i * L); ++j) t += G[j] * f[i * L + j];
        dft(t), dft(H);
        for (int j = 0; j < k; ++j) h[j] = (h[j] + (ll)t[j] * H[j]) % P, H[j] = mul(H[j], g[j]);
    }
    return idft(h), h;
}
// get F(x) for F(G(x)) = x
Poly compoundInv(Poly g) { // (4 * sqrt(n) + 2) * M(n) + n ^ 2
    int n = g.size(), k = norm(2 * n), L = sqrt(n) + 1;
    vector<Poly> G(L + 1);
    Poly H, f(n);
    auto dft = [&](Poly &a) { a.resize(k), NTT::dft(a.data(), k); };
    auto idft = [&](Poly &a) { NTT::idft(a.data(), k), a.resize(n); };
    G[0] = H = {1}, H.resize(n), dft(g = (g >> 1).inv(n));
    for (int i = 1; i <= L; ++i) dft(G[i] = G[i - 1]), idft(G[i] ^= g);
    dft(g = G[L]);
    for (int i = 0, t = 1; i < L; ++i) {
        for (int j = 1; j <= L && t < n; ++j, ++t)
            for (int r = 0; r < t; ++r) f[t] = (f[t] + (ll)G[j][r] * H[t - 1 - r]) % P;
        dft(H), idft(H ^= g);
    }
    return f ^ Poly(inv, inv + N);
}

int main() { return 0; }