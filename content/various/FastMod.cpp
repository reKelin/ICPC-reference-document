/*
 * Description:  Compute a % b about 5 times faster than usual,
 * where b is constant but not known at compile time.
 */
#include <bits/stdc++.h>

using namespace std;
using ll = int64_t;

struct Z {
    static inline u64 p, ip;
    static void setMod(int m) { p = m, ip = -1ull / p; }
    static u64 mod(u64 a) {
        u64 r = a - (u64)(((__uint128_t)ip * a) >> 64) * p;
        return r >= p ? r - p : r;
    }
};

int main() { return 0; }