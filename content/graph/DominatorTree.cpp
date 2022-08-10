/*
 * Description: Solving domination tree of directed graph.
 * Time: O(n alpha(n))
 */

#include <bits/stdc++.h>

using namespace std;
struct domT {
    using graph = vector<vector<int>>;
    int n, dft;
    graph G, iG, sG;
    vector<int> dfn, ord, fa, p, mn, sdom, idom;
    void init(int _n) {
        n = _n, dft = 0;
        ord.reserve(n);
        G.resize(n), iG = sG = G;
        p.resize(n), dfn = fa = idom = p;
        iota(p.begin(), p.end(), 0), sdom = mn = p;
    }
    void add(int u, int v) { G[u].push_back(v), iG[v].push_back(u); }
    void dfs(int u) {
        dfn[u] = ++dft, ord.push_back(u);
        for (int v : G[u])
            if (!dfn[v])
                fa[v] = u, dfs(v);
    }
    int GF(int u) {
        if (u == p[u]) return u;
        int &f = p[u], pu = GF(f);
        if (dfn[sdom[mn[f]]] < dfn[sdom[mn[u]]])
            mn[u] = mn[f];
        return f = pu;
    }
    void LengauerTarjan(int rt) {
        dfs(rt);
        for (int i = ord.size() - 1, u; i; --i) {
            u = ord[i];            
            for (int v : iG[u]) {
                if (!dfn[v]) continue;
                GF(v);
                if (dfn[sdom[mn[v]]] < dfn[sdom[u]])
                    sdom[u] = sdom[mn[v]];
            }
            sG[sdom[u]].push_back(u);
            u = p[u] = fa[u];
            for (auto v : sG[u])
                GF(v), idom[v] = u == sdom[mn[v]] ? u : mn[v];
            sG[u].clear();
        }
        for (auto u : ord)
            if (idom[u] != sdom[u])
                idom[u] = idom[idom[u]];
    }
    graph work(int rt) {
        LengauerTarjan(rt);
        graph res(n);
        for (int i = 0; i < n; ++i)
            if (idom[i] != i) res[idom[i]].push_back(i);
        return res;
    }
};

int main() { return 0; }