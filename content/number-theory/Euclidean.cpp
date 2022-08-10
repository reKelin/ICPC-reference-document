#include <bits/stdc++.h>
#define fp(i, a, b) for (int i = (a), i##_ = (b) + 1; i < i##_; ++i)
#define fd(i, a, b) for (int i = (a), i##_ = (b) - 1; i > i##_; --i)

using namespace std;
using ll = int64_t;
const int P = 998244353;
int POW(ll a, int b = P - 2, ll x = 1) {
    for (; b; b >>= 1, a = a * a % P)
        if (b & 1) x = x * a % P;
    return x;
}

template<const int N = 11>
struct Eucildean {
    int f[90][N][N], B[N][N], C[N][N];
    Eucildean() {
        C[0][0] = 1;
        fp(i, 1, N - 1) {
            C[i][0] = 1;
            fp(j, 1, i) C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]) % P;
        }
        A[0][1] = A[0][0] = 1;
        fp(i, 1, N - 2) {
            fp(j, 0, i + 1) A[i][j] = C[i + 1][j];
            fp(j, 0, i - 1) fp(k, 0, j + 1)
                A[i][k] = (A[i][k] + ll(P - C[i + 1][j]) * A[j][k]) % P;
            ll iv = POW(i + 1);
            fp(j, 0, i + 1) A[i][j] = A[i][j] * iv % P;
        }
    }
    int S(ll n, int m) {
        int s = 0, pw = 1; n %= P;
        fp(i, 0, m + 1) s = (s + (ll)A[m][i] * pw) % P, pw = pw * n % P;
        return s;
    }
    int calc(int k1, int k2, int a, int b, int c, ll n, int d) {
        if (f[d][k1][k2] != -1) return f[d][k1][k2];
        int &F = f[d][k1][k2] = 0;
        if (a == 0) return F = POW(b / c, k2, S(n, k1)) % P;
        if (a >= c || b >= c) {
            fp(i, 0, k2) fp(j, 0, k2 - i)
                F = (F + (ll)C[k2][i] * C[k2 - i][j] % P * POW(a / c, i) % P * POW(b / c, j) % P *
                     calc(k1 + i, k2 - i - j, a % c, b % c, c, n, d + 1)) % P;
            return F;
        }
        ll m = (a * n + b) / c - 1;
        F = POW((m + 1) % P, k2, S(n, k1));
        fp(i, 0, k2 - 1) fp(j, 0, k1 + 1)
            F = (F + ll(P - C[k2][i]) * A[k1][j] % P * calc(i, j, c, c - b - 1, a, m, d + 1)) % P;
        return F;
    }
    // sum x ^ k1 * [(ax + b) / c] ^ k2, 0 <= x <= n
    int solve(int k1, int k2, ll a, ll b, ll c, ll n) {
        return memset(f, -1, sizeof f), calc(k1, k2, a, b, c, n, 0);
    }
};

// Correctness test
int n, a, b, c;
void LOJ138() {
    Eucildean<> E{};
    scanf("%*d");
    for (int k1, k2; ~scanf("%d%d%d%d%d%d", &n, &a, &b, &c, &k1, &k2);)
        printf("%d\n", E.solve(k1, k2, a, b, c, n));
}
void P5170() {
    Eucildean<5> E{};
    scanf("%*d");
    while (~scanf("%d%d%d%d", &n, &a, &b, &c)) {
        E.solve(0, 1, a, b, c, n), E.calc(0, 2, a, b, c, n, 0), E.calc(1, 1, a, b, c, n, 0);
        printf("%d %d %d\n", E.f[0][0][1], E.f[0][0][2], E.f[0][1][1]);
    }
}
int main() {
    P5170();
    return 0;
}