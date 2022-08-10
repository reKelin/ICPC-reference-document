#include <bits/stdc++.h>

using namespace std;
using ll = int64_t;
const int P = 1e9 + 7;

// O(n ^ {3 / 4} / log n)
template <class PF, PF &F>
struct Min_25 {
    vector<ll> w, G, pr;
    ll n; int m, sqrtn;
    int id(ll x) { return x <= sqrtn ? x : m - n / x; }
    void init(ll _n) {
        n = _n, m = 1, sqrtn = sqrt(n) + 1;
        w.resize(sqrtn * 2), G = w, pr = {};
        for (ll i = 1; i <= n; ++i, ++m)
            w[m] = i = n / (n / i), G[m] = F.G(i); // copy G
        for (int p = 2; p <= sqrtn; ++p)
            if (G[p] != G[p - 1]) {
                double ip = 1. / p; pr.push_back(p);
                for (int i = m - 1, j = id((ll)p * p); i >= j; --i)
                    G[i] = (G[i] + F.g(p) * (G[p - 1] - G[id(w[i] * ip + 1e-9)] + P)) % P; // copy G
            }
        pr.push_back(sqrtn + 1); // prevent array access out of bounds
    }
    int S(ll x, int y) {
        int res = (G[id(x)] - G[pr[y] - 1] + P) % P; // res = sum a_i * (G_i(x) - G_i(pr_y - 1))
        for (int i = y; pr[i] * pr[i] <= x; ++i)
            for (ll k = 1, p = pr[i], q = p; q <= x / p; q *= p, ++k)
                res = (res + F.f(p, k + 1, q * p) + F.f(p, k, q) * S(x / q, i + 1)) % P;
        return res;
    }
    int S(ll x) { return init(x), S(x, 0) + 1; }
};

// modify  function
struct fun {
    static ll g(int p) { return p; }                                // f(p) = poly(p) = sum a_i * g_i(p), g_i(p) = p ^ i
    static ll G(ll n) { return n %= P, (n * (n + 1) / 2 - 1) % P; } // G(n) = sum g - g(1)
    static ll f(int p, int k, ll pk) { return k == 1 ? p : 0; }     // f(p ^ k) can be calculated quickly
} Fun;

Min_25<fun, Fun> m25;

void Gym103428I() { // Correctness test
    ll n;
    scanf("%lld", &n), m25.init(n);
    int m = sqrt(n) + 0.5, ans = 0;
    auto calc = [&](ll x) { return x % P * ((n - x) % P) % P; };
    for (auto p : m25.pr) {
        ll s = 0, q = p * p;
        for (; q <= n; q *= p) s += calc(n / q);
        ans = (ans + s % P * p) % P;
    }
    for (ll l = 1, r; l <= n; l = r + 1) {
        r = n / (n / l);
        ans = (ans + calc(n / l) * (m25.G[m25.id(r)] - m25.G[m25.id(l - 1)] + P)) % P;
    }
    printf("%d\n", ans * 2 % P);
}

int main() {
    Gym103428I();
    return 0;
}