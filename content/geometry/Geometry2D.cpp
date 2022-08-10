#include <bits/stdc++.h>

using namespace std;

using T = double;
using db = double;

const T eps = 1e-8, pi = acos(-1);

int sgn(T x) { return (x > eps) - (x < -eps); }

struct Vec {
    T x, y;
    bool operator<(Vec p) const { return tie(x, y) < tie(p.x, p.y); }
    bool operator==(Vec p) const { return tie(x, y) == tie(p.x, p.y); }
    Vec operator+(Vec p) const { return {x + p.x, y + p.y}; }
    Vec operator-(Vec p) const { return {x - p.x, y - p.y}; }
    Vec operator*(T d) const { return {x * d, y * d}; }
    Vec operator/(T d) const { return {x / d, y / d}; }
    T operator*(Vec p) const { return x * p.x + y * p.y; }
    T cross(Vec p) const { return x * p.y - y * p.x; }
    T cross(Vec a, Vec b) const { return (a - *this).cross(b - *this); }
    db len() const { return sqrt(x * x + y * y); }
    int onLeft(Vec p) const { return sgn(cross(p)); }
    int half() const { return y > 0 || (y == 0 && x > 0) ? 1 : -1; }
    friend bool argcmp(Vec a, Vec b) { return a.half() != b.half() ? a.half() == 1 : a.cross(b) > eps; }
    // rarely used
    T len2() const { return x * x + y * y; }
    db angle() const { return atan2(y, x); } // angle to x-axis in [-pi, pi]
    Vec unit() const { return *this / len(); }
    Vec perp() const { return {-y, x}; } // rotates 90 degrees ccw
    Vec norm() const { return perp().unit(); }
    Vec rot(db a) const { return {x * cos(a) - y * sin(a), x * sin(a) + y * cos(a)}; } // rotates a radians ccw
};

using Poly = vector<Vec>; // polygon, points, vectors

struct Line {
    Vec p, v; // point on line, direction vector
    int onLeft(Vec a) const { return v.onLeft(a - p); }
    Vec inter(const Line &a) const { return p + v * (a.v.cross(p - a.p) / (v.cross(a.v))); }
    db dis(Vec a) const { return v.cross(a - p) / v.len(); }
    Vec proj(Vec a) const { return p + v * (v * (a - p) / v.len2()); }
};

struct Seg {
    Vec a, b; // endpoint of line segment
    int on(Vec p) const { return p == a || p == b || (p.cross(a, b) == 0 && (p - a) * (p - b) < -eps); }
    Poly inter(const Line &l) const { 
        if (l.onLeft(a) * l.onLeft(b) == -1)
            return {{l.inter({a, b - a})}};
        set<Vec> t; // l.v != {0, 0}
        if (!l.onLeft(a)) t.insert(a);
        if (!l.onLeft(b)) t.insert(b);
        return Poly(t.begin(), t.end());
    }
    Poly inter(const Seg &s) const {
        const Line u{a, b - a}, v{s.a, s.b - s.a};
        if (u.onLeft(s.a) * u.onLeft(s.b) == -1 && v.onLeft(a) * v.onLeft(b) == -1)
            return {u.inter(v)};
        set<Vec> t;
        if (on(s.a)) t.insert(s.a);
        if (on(s.b)) t.insert(s.b);
        if (s.on(a)) t.insert(a);
        if (s.on(b)) t.insert(b);
        return Poly(t.begin(), t.end());
    }
    db dis(Vec p) const {
        if ((p - a) * (b - a) < -eps || (p - b) * (a - b) < -eps)
            return min((p - a).len(), (p - b).len());
        return Line{a, b - a}.dis(p);
    }
};

// polygon area: O(n)
T area(const Poly &p) {
    if (p.empty()) return 0;
    T a = p.back().cross(p[0]);
    for (unsigned i = 1; i < p.size(); ++i)
        a += p[i - 1].cross(p[i]);
    return a / 2;
}

// point inside polygon: O(n)
int inPolygon(const Poly &p, Vec a, int strict = 1) {
    int t = 0, n = p.size();
    for (int i = 0; i < n; ++i) {
        Vec q = p[(i + 1) % n];
        if (Seg{p[i], q}.on(a)) return !strict;
        t ^= ((a.y < p[i].y) - (a.y < q.y)) * a.cross(p[i], q) > 0;
    }
    return t;
}

// point inside convex hull: O(logn)
int inHull(const Poly &h, Vec p, int strict = 1) {
    if (h.empty()) return 0;
    int a = 1, b = h.size() - 1, c, r = !strict;
    if (b < 2) return r && Seg{h[0], h[b]}.on(p);
    auto left = [](Vec a, Vec b, Vec p) { return sgn(a.cross(b, p)); };
    if (left(h[0], h[a], h[b]) > 0) swap(a, b);
    if (left(h[0], h[a], p) >= r || left(h[0], h[b], p) <= -r) return 0;
    while (abs(a - b) > 1) c = (a + b) >> 1, (left(h[0], h[c], p) > 0 ? b : a) = c;
    return left(h[a], h[b], p) < r;
}

// convex hull: O(nlogn)
Poly convexHull(Poly p) {
    if (p.size() <= 1) return p;
    sort(p.begin(), p.end());
    Poly h(p.size() + 1);
    int s = 0, t = 0, i = 2;
    for (; i--; s = --t, reverse(p.begin(), p.end()))
        for (auto a : p) {
            while (t >= s + 2 && h[t - 2].cross(h[t - 1], a) <= 0) t--;
            h[t++] = a;
        }
    return {h.begin(), h.begin() + t - (t == 2 && h[0] == h[1])};
}

// halfplane intersection: O(nlogn)
Poly hpInter(vector<Line> L, const T Inf = 1e9) {
    static vector<Line> b = {{{Inf, -Inf}, {0, 1}}, {{Inf, Inf}, {-1, 0}}, {{-Inf, Inf}, {0, -1}}, {{-Inf, -Inf}, {1, 0}}};
    L.insert(L.end(), b.begin(), b.end());
    sort(L.begin(), L.end(), [](Line a, Line b) { return argcmp(a.v, b.v); });
    deque<Vec> ps; deque<Line> ls;
    for (auto l : L) {
        if (ls.empty()) {
            ls.push_back(l);
            continue;
        }
        while (!ps.empty() && l.onLeft(ps.back()) != 1)
            ps.pop_back(), ls.pop_back();
        while (!ps.empty() && l.onLeft(ps[0]) != 1)
            ps.pop_front(), ls.pop_front();
        if (l.v.cross(ls.back().v) == 0) {
            if (l.v * ls.back().v > 0) {
                if (l.onLeft(ls.back().p) != 1)
                    ls[0] = l;
                continue;
            }
            return {};
        }
        ps.push_back(l.inter(ls.back())), ls.push_back(l);
    }
    while (!ps.empty() && ls[0].onLeft(ps.back()) != 1)
        ps.pop_back(), ls.pop_back();
    if (ls.size() <= 2) return {};
    ps.push_back(ls[0].inter(ls.back()));
    return Poly(ps.begin(), ps.end());
}

// mincowsky sum: O(n + m)
Poly mincowskySum(Poly a, Poly b) {
    int n = a.size(), m = b.size();
    Poly c(n + m + 1), sa(n), sb(m);
    auto cmp = [](Vec a, Vec b) { return tie(a.y, a.x) < tie(b.y, b.x); };
    rotate(a.begin(), min_element(a.begin(), a.end(), cmp), a.end());
    rotate(b.begin(), min_element(b.begin(), b.end(), cmp), b.end());
    c[0] = a[0] + b[0];
    for(int i = 0; i < n; ++i) sa[i] = a[(i + 1) % n] - a[i];
    for(int i = 0; i < m; ++i) sb[i] = b[(i + 1) % m] - b[i];
    merge(sa.begin(), sa.end(), sb.begin(), sb.end(), c.begin() + 1, [](Vec a, Vec b) { return argcmp(a, b); });
    for(int i = 1; i <= n + m; ++i) c[i] = c[i] + c[i - 1];
    return c;
}

// axes of polygon: O(n)
vector<Line> polygonAxes(Poly a) {
    int n = a.size(), m = 3 * n;
    a.push_back(a[0]), a.push_back(a[1]);
    vector<pair<T, T>> s(n * 2);
    for (int i = 1, j = 0; i <= n; ++i, j += 2)
        s[j] = {(a[i] - a[i - 1]).len2(), 0},
        s[j + 1] = {a[i].cross(a[i - 1], a[i + 1]), (a[i + 1] - a[i]) * (a[i - 1] - a[i])};
    s.insert(s.end(), s.begin(), s.begin() + n);
    vector<int> f(m), res;
    for (int r = 0, p = 0, i = 0; i < m; i++) {
        f[i] = r > i ? min(r - i, f[p * 2 - i]) : 1;
        while (i >= f[i] && i + f[i] < m && s[i - f[i]] == s[i + f[i]]) ++f[i];
        if (i + f[i] > r) r = i + f[i], p = i;
        if (f[i] > n) res.push_back(i);
    }
    auto get = [&](int i) {
        int x = (i + 1) / 2;
        return i & 1 ? a[x] : (a[x] + a[x + 1]) / 2;
    };
    vector<Line> axe;
    for (auto i : res)
        axe.push_back({get(i), get(i) - get(i - n)});
    return axe;
}

int main() { return 0; };