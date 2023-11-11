struct odt : map<int, ll> {
    odt::iterator find(int p) { return --upper_bound(p); }
    void split(int p) { this->insert({p, find(p)->second}); }
    template<typename F> void work(int l, int r, F fun) { // [l, r)
        split(l), split(r);
        auto i = find(l), j = find(r);
        while (i != j) fun(i++);
    }
};