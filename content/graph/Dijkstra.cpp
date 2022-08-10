/*
 * Description: Single source shortest path without negative weight edge.
 * Time: O(m log n)
 */
#include <bits/stdc++.h>

using namespace std;
using T = int64_t;
const T Inf = numeric_limits<T>::max();
const int N = 2e5 + 5;

int n;
vector<pair<int, int>> G[N];

vector<T> dijkstra(int S) {
    bitset<N> vis;
    vector<T> dis(n + 1, Inf);
    priority_queue<pair<T, int>> q;
    q.push({dis[S] = 0, S});
    while (!q.empty()) {
        int u = q.top().second;
        q.pop();
        if (vis[u]) continue;
        vis[u] = 1;
        for (auto [v, w] : G[u])
            if (dis[v] > dis[u] + w)
                q.push({-(dis[v] = dis[u] + w), v});
    }
    return dis;
}

int main() { return 0; }