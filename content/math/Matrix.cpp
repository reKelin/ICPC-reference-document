#include <bits/stdc++.h>
#define fp(i, a, b) for(int i = a, i##_ = (b) + 1; i < i##_; ++i)

using namespace std;
using ll = int64_t;
using db = double;
// const int P = 998244353;
int P;
int POW(ll a, int b = P - 2, ll x = 1) {
    for (; b; b >>= 1, a = a * a % P)
        if (b & 1) x = x * a % P;
    return x;
}

struct Z {
    using u64 = uint64_t;
    static inline u64 p, ip;
    static void setMod(int m) { p = m, ip = -1ull / p; }
    static u64 mod(u64 a) {
        u64 r = a - (u64)(((__uint128_t)ip * a) >> 64) * p;
        return r >= p ? r - p : r;
    }
};

// Defualt Time: O(n^3)
struct Mat : vector<vector<int>> {
#define T (*this)
    Mat(int n = 0, int m = 0) { T.assign(n, vector<int>(m)); }
    static Mat I(int n) {
        Mat a(n, n);
        fp(i, 0, n - 1) a[i][i] = 1;
        return a;
    }

    Mat operator*(const Mat &a) const {
        const int n = size(), o = T[0].size(), m = a[0].size();
        Mat c(n, m);
        fp(k, 0, o - 1) fp(i, 0, n - 1) if (T[i][k]) fp(j, 0, m - 1)
            c[i][j] = (c[i][j] + 1ll * T[i][k] * a[k][j]) % P;
        return c;
    }

    Mat pow(ll b) const {
        Mat a = T, x = I(size());
        for (; b; b >>= 1, a = a * a)
            if (b & 1) x = x * a;
        return x;
    }

#ifndef prime
    // P no limit
    int det() const {
        const int n = size();
        Mat a = T; int d = 1;
        fp(i, 0, n - 1) {
            auto x = a[0].end(), y = x;
            fp(j, i, n - 1) if (a[j][i] && (x == y || a[j][i] < x[i])) x = a[j].begin();
            if (x == y) return 0;
            fp(j, i, n - 1)
                if (a[j].begin() != x && a[j][i])
                    for (y = a[j].begin();; std::swap(x, y)) {
                        int w = P - y[i] / x[i];
                        fp(k, i, n - 1) y[k] = (y[k] + (ll) x[k] * w) % P;
                        if (!y[i]) break;
                    }
            if (x != a[i].begin()){
                d ^= 1;
                fp(k, i, n - 1) std::swap(x[k], a[i][k]);
            }
        }
        if (!d) d = P - 1;
        fp(i, 0, n - 1) d = ((ll)d * a[i][i]) % P;
        return d;
    }
#else
    // P is a prime number
    int det() const {
        const int n = size();
        Mat a = T; int d = 1, ia, w;
        fp(i, 0, n - 1) {
            if (!a[i][i])
                fp(j, i, n - 1)
                    if (a[j][i]) {
                        std::swap(a[i], a[j]), d ^= 1;
                        break;
                    }
            if (!a[i][i]) return 0;
            ia = P - POW(a[i][i]);
            fp(j, i + 1, n - 1) {
                w = (ll)a[j][i] * ia % P;
                fp(k, i, n - 1) a[j][k] = (a[j][k] + (ll)a[i][k] * w) % P;
            }
        }
        if (!d) d = P - 1;
        fp(i, 0, n - 1) d = (ll) d * a[i][i] % P;
        return d;
    }
#endif

    // P is a prime number
    Mat inv() const {
        const int n = size();
        Mat a = T;
        vector<int> c(n);
        iota(c.begin(), c.end(), 0);
        for(int i = 0, p; i < n; ++i) {
            for (p = i; p < n && !a[p][i]; ++p);
            if (p == n) return {};
            std::swap(a[i], a[p]), std::swap(c[i], c[p]);
            ll w = POW(a[i][i]); a[i][i] = 1;
            for (auto &x : a[i]) x = x * w % P;
            for (auto &r : a) if (&r != &a[i] && r[i]) {
                w = P - r[i], r[i] = 0;
                fp(j, 0, n - 1) r[j] = (r[j] + a[i][j] * w) % P;
            }
        }
        fp(i, 0, n - 1)
            for (int j = c[i]; j != i; std::swap(c[i], c[j]), j = c[i])
                fp(k, 0, n - 1) std::swap(a[k][i], a[k][j]);
        return a;
    }

    // Turn into upper hessenberg matrix
    Mat hessenberg() const {
        const int n = size();
        Mat a = T;
        fp(i, 0, n - 3) {
            if (!a[i + 1][i])
                fp(j, i + 2, n - 1) if (a[j][i]) {
                    fp(k, i, n - 1) std::swap(a[i + 1][k], a[j][k]);
                    fp(k, 0, n - 1) std::swap(a[k][i + 1], a[k][j]);
                    break;
                }
            if (!a[i + 1][i]) continue;
            int inv = POW(a[i + 1][i]);
            fp(j, i + 2, n - 1) {
                int w = (ll) inv * a[j][i] % P, nw = P - w;
                fp(k, i, n - 1) a[j][k] = (a[j][k] + (ll) a[i + 1][k] * nw) % P;
                fp(k, 0, n - 1) a[k][i + 1] = (a[k][i + 1] + (ll) a[k][j] * w) % P;
            }
        }
        return a;
    }

    // p_A(x) = det(xI_n - A)
    vector<int> characteristic() const {
        const int n = size();
        Mat a = T.hessenberg();
        Mat g(n + 1); g[0] = {1};
        for (auto &r: a) for (auto &x: r) x = P - x;
        fp(i, 1, n) {
            g[i] = g[i - 1], g[i].insert(g[i].begin(), 0);
            for (int j = i - 1, w = 1, wa; ~j; --j) {
                wa = (ll) w * a[j][i - 1] % P;
                fp(k, 0, j) g[i][k] = (g[i][k] + (ll) g[j][k] * wa) % P;
                if (j) w = (ll) w * (P - a[j][j - 1]) % P;
            }
        }
        return g[n];
    }

    // solve Ax = b
    vector<db> solveLinear(vector<int> b) {
        const int n = size();
        vector<vector<db>> a(n, vector<db>(n + 1));
        fp(i, 0, n - 1) a[i][n] = b[i];
        fp(i, 0, n - 1) fp(j, 0, n - 1) a[i][j] = T[i][j];
        for (int i = 0, p = 0; i < n; ++i, p = i) {
            fp(k, i, n - 1) if (abs(a[k][i]) > abs(a[p][i])) p = k;
            if (abs(a[p][i]) < 1e-8) return {};
            std::swap(a[i], a[p]);
            db w = 1 / a[i][i]; 
            fp(j, i, n) a[i][j] *= w;
            fp(k, 0, n - 1) if (k != i) {
                w = a[k][i];
                fp(j, i, n) a[k][j] -= a[i][j] * w;
            }
        }
        vector<db> ans(n);
        fp(i, 0, n - 1) ans[i] = a[i][n];
        return ans;
    }
#undef T
};

// correctness test
int n;
void P3389() {
    scanf("%d", &n);
    Mat a(n, n);
    vector<int> b(n);
    fp(i, 0, n - 1) {
        for (auto &x : a[i]) scanf("%d", &x);
        scanf("%d", &b[i]);
    }
    auto r = a.solveLinear(b);
    if (r.empty()) puts("No Solution");
    else for (auto x : r) printf("%.2lf\n", x);
}
void P4783() { // P =  1e9 + 7
    scanf("%d", &n);
    Mat a(n, n);
    for (auto &r : a) for (auto &x: r) scanf("%d", &x);
    a = a.inv();
    if (a.empty()) return puts("No Solution"), void();
    for (auto &r : a) fp(j, 0, n - 1) printf("%d%c", r[j], " \n"[j == n - 1]);
}
void P7112(){ // P no limit
    scanf("%d%d", &n, &P);
    // Z::setMod(P);
    Mat a(n, n);
    for (auto &r : a) for (auto &x : r) scanf("%d", &x), x %= P;
    printf("%d\n", a.det());
}
void P7776(){ // P = 998244353
    scanf("%d", &n);
    Mat a(n, n);
    for (auto &r : a) for (auto &x : r) scanf("%d", &x);
    auto f = a.characteristic();
    fp(i, 0, n) printf("%d%c", f[i], " \n"[i == n]);
}
int main() {
    P3389();
    return 0;
}