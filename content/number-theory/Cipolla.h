#pragma once

#include <bits/stdc++.h>

using namespace std;
using ll = int64_t;

// O(log P)
struct Cipolla {
    int P, I2;
    using pll = pair<ll, ll>;
#define X first
#define Y second
    ll mul(ll a, ll b) const { return a * b % P; }
    pll mul(pll a, pll b) const {
        return {(a.X * b.X + I2 * a.Y % P * b.Y) % P,
                (a.X * b.Y + a.Y * b.X) % P};
    }
    template <class T> T POW(T a, int b, T x) {
        for (; b; b >>= 1, a = mul(a, a))
            if (b & 1) x = mul(x, a);
        return x;
    }

    Cipolla(int p = 0) : P(p) {}
    pair<int, int> sqrt(int n) {
        int a = rand(), x;
        if (!(n %= P)) return {0, 0};
        if (POW(n, (P - 1) >> 1, 1) == P - 1) return {-1, -1};
        while (POW(I2 = ((ll)a * a - n + P) % P, (P - 1) >> 1, 1) == 1) a = rand();
        x = (int)POW(pll{a, 1}, (P + 1) >> 1, {1, 0}).X;
        if (2 * x > P) x = P - x;
        return {x, P - x};
    }
#undef X
#undef Y
};
