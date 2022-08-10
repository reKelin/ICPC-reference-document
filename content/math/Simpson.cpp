/*
 * Description: calculate $\int_l^r f(x) dx$
 * Time: O(?)
 */
#include <bits/stdc++.h>

using namespace std;
using db = double;
const db eps = 1e-10;
struct Simpson {
    function<db(db)> f;
    db simp(db l, db r) const { return (f(l) + 4 * f((l + r) * 0.5) + f(r)) * (r - l) / 6; }
    db calc(db l, db r, db S) {
        db m = (l + r) * 0.5, sl = simp(l, m), sr = simp(m, r);
        if (abs(S - sl - sr) < 15 * eps) return sl + sr + (sl + sr - S) / 15;
        return calc(l, m, sl) + calc(m, r, sr);
    }
    template <class fun>
    db simpson(fun F, db l, db r) { return f = F, calc(l, r, simp(l, r)); }
};

// correctness test
void P4526 () {
    db a;
    scanf("%lf", &a);
    auto f = [&](db x) { return pow(x, (a / x) - x); };
    if (a < 0) puts("orz");
    else printf("%.5lf", (new Simpson)->simpson(f, eps, 15));
}
int main() {
    P4526();
    return 0;
}