/*
 * Desciption: Dinic with current edge optimization
 * Time: O(n^2 m)
 */
#include <bits/stdc++.h>
using namespace std;
const int N = 100 + 5;

template<const int N = N, class T = int>
struct Dinic {
    const T Inf = numeric_limits<T>::max();
    struct Edge { int v, rev; T cap; };
    int n, st, ed; T Flow;
    vector<int> d, cur;
    vector<Edge> G[N];
    bool bfs() {
        d.assign(n + 1, -1), d[st] = 0;
        queue<int> q; q.push(st);
        while (!q.empty()) {
            int u = q.front(), v; q.pop();
            for (auto E : G[u])
                if (E.cap && d[v = E.v] == -1)
                    d[v] = d[u] + 1, q.push(v);
        }
        return d[ed] != -1;
    }
    T dfs(int u, T lim) {
        if (u == ed) return lim;
        T f = 0, v; auto &E = G[u];
        for (auto &i = cur[u]; i < (int)E.size(); ++i)
            if (d[v = E[i].v] == d[u] + 1 && E[i].cap) {
                int t = dfs(v, min(lim - f, E[i].cap));
                E[i].cap -= t, G[v][E[i].rev].cap += t, f += t;
                if (f == lim) return f;
            }
        return f;
    }
    void add(int u, int v, T a) {
        G[u].push_back({v, (int)G[v].size(), a});
        G[v].push_back({u, (int)G[u].size() - 1, 0});
    }
    T work(int n_, int s, int t) {
        n = n_, st = s, ed = t, Flow = 0;
        for (T f; bfs();) {
            cur.assign(n + 1, 0);
            do Flow += f = dfs(st, Inf); while (f);
        }
        return Flow;
    }
};

int main() { return 0; }