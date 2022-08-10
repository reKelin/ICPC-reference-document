#include <bits/stdc++.h>

using namespace std;
using ll = int64_t;
const int N = 1e7, P = 1e9 + 7;

#define inc(a, b) (((a) += (b)) >= P ? (a) -= P : 0)
int POW(ll a, int b = P - 2, ll x = 1) {
    for (; b; b >>= 1, a = a * a % P)
        if (b & 1) x = x * a % P;
    return x;
}

vector<int> pr;
void sieves() {
    bitset<N> vis;
    for (int i = 2; i < N; ++i) {
        if (!vis[i]) pr.push_back(i);
        for (auto p: pr) {
            if (i * p > N) break;
            vis[i * p] = 1;
            if (i % p == 0) break;
        }
    }
    pr.push_back(N + 1); // prevent array access out of bounds
}

// The complexity is O(sqrt n * H) if and only if f(p) = p ^ k, \
    where H is the complexity of calculating h(p ^ k), otherwise Min_25
template <class PF, PF &F>
struct PN {
    int S(ll x, int y) {
        int res = F.G(x);
        for (int i = y; (ll)pr[i] * pr[i] <= x; ++i)
            for (ll k = 2, p = pr[i], q = p; q <= x / p; ++k)
                q *= p, res = (res + F.h(p, k, q) * S(x / q, i + 1)) % P;
        return res;
    }
};

vector<vector<int>> f;

struct fun {
    static ll g(int p) { return 1; }                                              // g(p) = f(p) = p ^ c
    static ll G(ll n) { return n % P; }                                           // G(n) = sum g
    static ll h(int p, int k, ll pk) { return (f[p][k] - f[p][k - 1] + P) % P; }  // f(p ^ k) = sum g(p ^ {k - i}) * h(p ^ i) -> h(p ^ k)
} F;

PN<fun, F> pn;

void Gym103306F() { // Correctness test
    f.resize(N);
    const ll M = 1e14;
    for (int i = 0; i + 1 < pr.size(); ++i) {
        ll p = pr[i], q = p;
        f[p] = {0, 1};
        for (int j = 2; q <= M / p; ++j)
            f[p].push_back(POW(j, p)), q *= p;
    }
    ll n;
    scanf("%lld", &n);
    printf("%d\n", pn.S(n, 0));
}
int main() {
    sieves();
    Gym103306F();
    return 0;
}