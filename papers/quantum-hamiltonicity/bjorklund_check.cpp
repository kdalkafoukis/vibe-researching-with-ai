// End-to-end check of Björklund's general-case test (arXiv:1008.0541, Lemmas 6, 9 and 11) on small
// graphs: a C++ port of bjorklund_check.py, fast enough to go up to 12 vertices.
//
// Part 1 reruns the graphs of bjorklund_check.py and reproduces its 2,560 results bit for bit. It
// uses the same field and the same random numbers: Python's Mersenne Twister seeding, randrange,
// and even CPython's set iteration order, which decides the order of the lazy random draws.
// Part 2 runs bigger graphs (10 and 12 vertices) on all cores, with one random stream per case.
//
// Build: clang++ -O3 -std=c++17 -pthread bjorklund_check.cpp -o bjorklund_check   (or g++)
// Run:   ./bjorklund_check           add --dump to print every case as "name V1 m p formula q"
//
// For each graph G, balanced partition V1|V2 and cutoff m, the labelled cycle cover sum
// Lambda(D, V2 ∪ L_m, f) is computed at a random point of GF(2^16) three ways, as in the Python:
//   p : exponential-space route, Lemma 6/7
//   q : polynomial-space route, Section 4 / Lemma 11 (walks instead of paths)
//   formula : closed form of Lemma 9 (sum over oriented Hamiltonian cycles with exactly m direct arcs)
// Soundness requires Lambda == 0 exactly whenever no Hamiltonian cycle has exactly m arcs inside V1.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using u32 = uint32_t;
using Edges = std::vector<std::pair<int, int>>;

// ---------------- GF(2^16) via log tables ----------------
static const int Q = 65535;
static const u32 POLY = 0x1100B;                        // x^16+x^12+x^3+x+1
static int EXPT[2 * Q], LOGT[Q + 1];
static void gf_init() {
    u32 v = 1;
    for (int i = 0; i < Q; i++) {
        if (i > 0 && v == 1) { fprintf(stderr, "modulus not primitive\n"); exit(1); }
        EXPT[i] = v; LOGT[v] = i; v <<= 1;
        if (v & 0x10000) v ^= POLY;
    }
    if (v != 1) { fprintf(stderr, "modulus not primitive\n"); exit(1); }
    for (int i = Q; i < 2 * Q; i++) EXPT[i] = EXPT[i - Q];
}
static inline u32 mul(u32 a, u32 b) { return (a && b) ? EXPT[LOGT[a] + LOGT[b]] : 0; }
static inline u32 inv(u32 a) { return EXPT[Q - LOGT[a]]; }

static const int MAXD = 8;
typedef u32 Mat[MAXD][MAXD];
static u32 det(const Mat Min, int n) {                  // char 2: det == permanent, row swaps free
    Mat M; memcpy(M, Min, sizeof(Mat));
    u32 d = 1;
    for (int c = 0; c < n; c++) {
        int p = -1;
        for (int i = c; i < n; i++) if (M[i][c]) { p = i; break; }
        if (p < 0) return 0;
        if (p != c) for (int j = 0; j < n; j++) std::swap(M[c][j], M[p][j]);
        u32 piv = M[c][c]; d = mul(d, piv); u32 ip = inv(piv);
        for (int i = c + 1; i < n; i++)
            if (M[i][c]) {
                u32 f = mul(M[i][c], ip);
                for (int j = c; j < n; j++) if (M[c][j]) M[i][j] ^= mul(f, M[c][j]);
            }
    }
    return d;
}

// coefficient of r^t of the polynomial through (xs[i], ys[i])
static u32 interp_coeff(const std::vector<u32>& xs, const std::vector<u32>& ys, int t) {
    std::vector<u32> N = {1};
    for (u32 a : xs) {                                  // N(r) = prod (r + a)
        std::vector<u32> R(N.size() + 1, 0);
        for (size_t i = 0; i < N.size(); i++) { R[i + 1] ^= N[i]; R[i] ^= mul(N[i], a); }
        N.swap(R);
    }
    int deg = (int)N.size() - 1; std::vector<u32> Qd(deg); u32 res = 0;
    for (size_t x = 0; x < xs.size(); x++) {
        u32 xi = xs[x];
        Qd[deg - 1] = N[deg];
        for (int k = deg - 1; k > 0; k--) Qd[k - 1] = N[k] ^ mul(xi, Qd[k]);
        u32 den = 0, pw = 1;
        for (u32 c : Qd) { den ^= mul(c, pw); pw = mul(pw, xi); }
        res ^= mul(ys[x], mul(Qd[t], inv(den)));
    }
    return res;
}

// ---------------- Python's random module, reproduced ----------------
struct PyRandom {                                       // MT19937 seeded like random.Random(seed)
    u32 mt[624]; int mti = 624;
    explicit PyRandom(const std::vector<u32>& key) {    // key = seed's 32-bit words, low word first
        mt[0] = 19650218u;
        for (int i = 1; i < 624; i++) mt[i] = 1812433253u * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i;
        int i = 1, j = 0, kl = (int)key.size();
        for (int k = std::max(624, kl); k; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525u)) + key[j] + j;
            if (++i >= 624) { mt[0] = mt[623]; i = 1; }
            if (++j >= kl) j = 0;
        }
        for (int k = 623; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941u)) - i;
            if (++i >= 624) { mt[0] = mt[623]; i = 1; }
        }
        mt[0] = 0x80000000u;
    }
    u32 next() {
        if (mti >= 624) {
            for (int k = 0; k < 624; k++) {
                u32 y = (mt[k] & 0x80000000u) | (mt[(k + 1) % 624] & 0x7fffffffu);
                mt[k] = mt[(k + 397) % 624] ^ (y >> 1) ^ ((y & 1) ? 0x9908b0dfu : 0u);
            }
            mti = 0;
        }
        u32 y = mt[mti++];
        y ^= y >> 11; y ^= (y << 7) & 0x9d2c5680u; y ^= (y << 15) & 0xefc60000u; y ^= y >> 18;
        return y;
    }
    double random() { u32 a = next() >> 5, b = next() >> 6; return (a * 67108864.0 + b) / 9007199254740992.0; }
    u32 randrange65536() { u32 r; do r = next() >> 15; while (r >= 65536); return r; }  // getrandbits(17), reject
};

// Iteration order of a CPython set of small non-negative ints, given the insertion order
// (open addressing: 8 slots, then 32 after the 5th insert; linear probes only in big tables).
static std::vector<int> py_set_order(const std::vector<int>& inserted) {
    std::vector<int> t(8, -1); size_t mask = 7, used = 0;
    auto put = [](std::vector<int>& tab, size_t msk, int v) {
        size_t i = (size_t)v & msk, perturb = (size_t)v;
        for (;;) {
            size_t probes = (i + 9 <= msk) ? 9 : 0;
            for (size_t e = i; ; e++) { if (tab[e] < 0) { tab[e] = v; return; } if (probes-- == 0) break; }
            perturb >>= 5; i = (i * 5 + 1 + perturb) & msk;
        }
    };
    for (int v : inserted) {
        if (std::find(t.begin(), t.end(), v) != t.end()) continue;
        put(t, mask, v); used++;
        if (used * 5 >= mask * 3) {
            size_t size = 8; while (size <= used * 4) size <<= 1;
            std::vector<int> nt(size, -1);
            for (int x : t) if (x >= 0) put(nt, size - 1, x);
            t.swap(nt); mask = size - 1;
        }
    }
    std::vector<int> out;
    for (int x : t) if (x >= 0) out.push_back(x);
    return out;
}

// ---------------- graphs and variables ----------------
struct Graph {
    int n; std::vector<u32> nb; std::vector<std::vector<int>> adjl;   // adjl: Python's iteration order
    Graph(int n_, const Edges& edges) : n(n_), nb(n_, 0) {
        std::vector<std::vector<int>> ins(n);
        for (auto [u, v] : edges) { ins[u].push_back(v); ins[v].push_back(u); nb[u] |= 1u << v; nb[v] |= 1u << u; }
        for (int u = 0; u < n; u++) adjl.push_back(py_set_order(ins[u]));
    }
    bool adj(int u, int v) const { return nb[u] >> v & 1; }
};

struct Vars {                                           // s-oriented mirror variables, drawn lazily
    PyRandom* rng; int s, n; std::vector<int32_t> xv, dv;
    Vars(PyRandom* r, int s_, int n_) : rng(r), s(s_), n(n_), xv(n_ * n_, -1), dv(n_ * n_ * n_, -1) {}
    int key(int u, int w) const { return (u == s || w == s) ? u * n + w : std::min(u, w) * n + std::max(u, w); }
    u32 get(int32_t& slot) { if (slot < 0) slot = (int32_t)rng->randrange65536(); return (u32)slot; }
    u32 x(int u, int w) { return get(xv[key(u, w)]); }             // traverse edge u -> w
    u32 xd(int u, int v, int d) { return get(dv[key(u, v) * n + d]); }  // direct arc u -> v, extra label d
};

struct Ctx {
    const Graph& G; Vars& var; int n, d, k; std::vector<int> V1, V2; int idx[32]; std::vector<int32_t> memo;
    Ctx(const Graph& g, const std::vector<int>& v1, Vars& vr) : G(g), var(vr), n(g.n), V1(v1) {
        std::fill(idx, idx + 32, -1);
        for (int x = 0; x < n; x++) if (!std::count(V1.begin(), V1.end(), x)) { idx[x] = (int)V2.size(); V2.push_back(x); }
        d = (int)V1.size(); k = (int)V2.size();
    }
    u32 h(int w, int v, int mask) {                     // simple paths w -> ... -> v, interior (incl. w) == mask
        int32_t& slot = memo[((size_t)idx[w] * n + v) * (1u << k) + mask];
        if (slot >= 0) return (u32)slot;
        int rest = mask & ~(1 << idx[w]); u32 val = 0;
        if (rest == 0) val = G.adj(w, v) ? var.x(w, v) : 0;
        else
            for (int z : G.adjl[w])
                if (idx[z] >= 0 && (rest >> idx[z] & 1)) { u32 t = h(z, v, rest); if (t) val ^= mul(var.x(w, z), t); }
        slot = (int32_t)val; return val;
    }
};

// ---------------- the three evaluations ----------------
static u32 lam_p(Ctx& c, int m) {
    const int d = c.d, k = c.k, K = 1 << k, Ls = k + m;
    c.memo.assign((size_t)k * c.n * K, -1);
    std::vector<u32> f((size_t)d * d * K, 0);
    for (int a = 0; a < d; a++)
        for (int b = 0; b < d; b++) {
            int u = c.V1[a], v = c.V1[b]; if (u == v) continue;
            for (int mask = 1; mask < K; mask++) {
                u32 val = 0;
                for (int w : c.G.adjl[u])
                    if (c.idx[w] >= 0 && (mask >> c.idx[w] & 1)) { u32 t = c.h(w, v, mask); if (t) val ^= mul(c.var.x(u, w), t); }
                f[((size_t)a * d + b) * K + mask] = val;
            }
        }
    // direct arcs, D[Yd][a][b] = sum over j in Yd of xd(u, v, j); same draw order as the Python
    std::vector<u32> D((size_t)(1 << m) * d * d, 0);
    for (int yd = 1; yd < (1 << m); yd++)
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++) {
                int u = c.V1[a], v = c.V1[b];
                if (u == v || !c.G.adj(u, v)) continue;
                for (int j = 0; j < m; j++) if (yd >> j & 1) D[((size_t)yd * d + a) * d + b] ^= c.var.xd(u, v, j);
            }
    const int npts = d * Ls + 1;
    std::vector<u32> xs(npts), ys(npts), S((size_t)d * d * K);
    for (int ri = 0; ri < npts; ri++) {
        u32 r = xs[ri] = EXPT[ri + 1], rp[64]; rp[0] = 1;
        for (int i = 1; i <= Ls; i++) rp[i] = mul(rp[i - 1], r);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++) {
                if (a == b) continue;
                u32* arr = &S[((size_t)a * d + b) * K]; const u32* fr = &f[((size_t)a * d + b) * K];
                arr[0] = 0;
                for (int mask = 1; mask < K; mask++) arr[mask] = mul(rp[__builtin_popcount(mask)], fr[mask]);
                for (int i = 0; i < k; i++)
                    for (int mask = 0; mask < K; mask++) if (mask >> i & 1) arr[mask] ^= arr[mask ^ (1 << i)];
            }
        u32 tot = 0; Mat M;
        for (int Y = 0; Y < (1 << Ls); Y++) {
            int Yv = Y & (K - 1), Yd = Y >> k;
            for (int a = 0; a < d; a++)
                for (int b = 0; b < d; b++)
                    M[a][b] = a == b ? 0 : S[((size_t)a * d + b) * K + Yv] ^ mul(r, D[((size_t)Yd * d + a) * d + b]);
            tot ^= det(M, d);
        }
        ys[ri] = tot;
    }
    return interp_coeff(xs, ys, Ls);
}

static u32 lam_q(Ctx& c, int m) {
    // Same walk sums as the Python; the walks are computed once per Y as coefficients of r^l,
    // then evaluated at every interpolation point (Horner), instead of once per point.
    const int d = c.d, k = c.k, n = c.n, Ls = k + m, npts = d * n + 1;
    std::vector<u32> xs(npts), ys(npts, 0), A((size_t)(n + 1) * d * d), Dm((size_t)d * d);
    for (int ri = 0; ri < npts; ri++) xs[ri] = EXPT[ri + 1];
    std::vector<std::pair<int, u32>> vec, nv; int pos[32];
    for (int Y = 0; Y < (1 << Ls); Y++) {
        std::vector<int> X1;
        for (int i = 0; i < k; i++) if (Y >> i & 1) X1.push_back(c.V2[i]);
        int Yd = Y >> k;
        std::fill(A.begin(), A.end(), 0); std::fill(Dm.begin(), Dm.end(), 0);
        for (int a = 0; a < d; a++) {
            int u = c.V1[a];
            vec.clear();
            for (int w : X1) if (c.G.adj(u, w)) vec.push_back({w, c.var.x(u, w)});
            for (int l = 1; l <= n; l++) {
                for (int b = 0; b < d; b++) {
                    int v = c.V1[b]; if (v == u) continue;
                    u32 s = 0;
                    for (auto [w, cw] : vec) if (c.G.adj(w, v)) s ^= mul(cw, c.var.x(w, v));
                    A[((size_t)l * d + a) * d + b] ^= s;
                }
                if (l == n) break;
                nv.clear(); std::fill(pos, pos + 32, -1);
                for (auto [w, cw] : vec)
                    for (int z : X1)
                        if (c.G.adj(w, z)) {
                            u32 val = mul(cw, c.var.x(w, z));
                            if (pos[z] < 0) { pos[z] = (int)nv.size(); nv.push_back({z, val}); }
                            else nv[pos[z]].second ^= val;
                        }
                vec.swap(nv);
            }
            for (int b = 0; b < d; b++) {
                int v = c.V1[b];
                if (v != u && c.G.adj(u, v))
                    for (int j = 0; j < m; j++) if (Yd >> j & 1) Dm[(size_t)a * d + b] ^= c.var.xd(u, v, j);
            }
        }
        Mat M;
        for (int ri = 0; ri < npts; ri++) {
            u32 r = xs[ri];
            for (int a = 0; a < d; a++)
                for (int b = 0; b < d; b++) {
                    u32 acc = A[((size_t)n * d + a) * d + b];
                    for (int l = n - 1; l >= 1; l--) acc = mul(acc, r) ^ A[((size_t)l * d + a) * d + b];
                    M[a][b] = mul(r, acc ^ Dm[(size_t)a * d + b]);
                }
            ys[ri] ^= det(M, d);
        }
    }
    return interp_coeff(xs, ys, Ls);
}

typedef std::vector<std::vector<int>> Cycles;
static Cycles oriented_hcs(const Graph& G) {
    int n = G.n; Cycles out; std::vector<int> p(n - 1);
    std::iota(p.begin(), p.end(), 1);
    do {
        if (p[0] > p[n - 2]) continue;
        bool ok = G.adj(0, p[0]) && G.adj(p[n - 2], 0);
        for (int i = 0; ok && i + 1 < n - 1; i++) ok = G.adj(p[i], p[i + 1]);
        if (!ok) continue;
        std::vector<int> cyc = {0}; cyc.insert(cyc.end(), p.begin(), p.end()); out.push_back(cyc);
        cyc = {0}; cyc.insert(cyc.end(), p.rbegin(), p.rend()); out.push_back(cyc);
    } while (std::next_permutation(p.begin(), p.end()));
    return out;
}

static u32 lam_formula(const Cycles& hcs, u32 in1, int m, Vars& var) {
    u32 tot = 0;
    for (const auto& H : hcs) {
        int n = (int)H.size(), cnt = 0;
        for (int i = 0; i < n; i++) cnt += (in1 >> H[i] & 1) && (in1 >> H[(i + 1) % n] & 1);
        if (cnt != m) continue;
        u32 base = 1; std::vector<std::pair<int, int>> U;
        for (int i = 0; i < n; i++) {
            int u = H[i], v = H[(i + 1) % n];
            if ((in1 >> u & 1) && (in1 >> v & 1)) U.push_back({u, v}); else base = mul(base, var.x(u, v));
        }
        u32 perm = 1;
        if (m > 0) {
            Mat P;
            for (int a = 0; a < m; a++) for (int j = 0; j < m; j++) P[a][j] = var.xd(U[a].first, U[a].second, j);
            perm = det(P, m);
        }
        tot ^= mul(base, perm);
    }
    return tot;
}

// ---------------- driver ----------------
struct Stats { long cases = 0, sound_zero = 0, nonzero_ok = 0, sz_miss = 0, mismatch = 0, false_pos = 0; };
struct Out { u32 p, f, q; };
static bool DUMP = false;
static Stats TOTAL[2];

// parallel == false: one random stream for the whole run, exactly as bjorklund_check.py.
// parallel == true : case (partition i, cutoff m) gets its own stream, seeded with the words (seed, i, m).
static std::map<int, int> run(const std::string& name, int n, const Edges& edges, bool use_q, u32 seed, bool parallel) {
    auto t0 = std::chrono::steady_clock::now();
    Graph G(n, edges); Cycles hcs = oriented_hcs(G);
    std::vector<std::vector<int>> parts; std::vector<int> sel(n, 0);
    std::fill(sel.begin(), sel.begin() + n / 2, 1);
    do {                                                // itertools.combinations order
        std::vector<int> V1; for (int i = 0; i < n; i++) if (sel[i]) V1.push_back(i);
        parts.push_back(V1);
    } while (std::prev_permutation(sel.begin(), sel.end()));
    const int M = n / 2;
    std::vector<std::vector<Out>> res(parts.size(), std::vector<Out>(M));
    auto one = [&](size_t i, int m, PyRandom& rng) {
        Vars var(&rng, parts[i][0], n); Ctx c(G, parts[i], var);
        u32 in1 = 0; for (int u : parts[i]) in1 |= 1u << u;
        u32 p = lam_p(c, m), f = lam_formula(hcs, in1, m, var), q = use_q ? lam_q(c, m) : p;
        res[i][m] = {p, f, q};
    };
    if (!parallel) {
        PyRandom rng({seed});
        for (size_t i = 0; i < parts.size(); i++) for (int m = 0; m < M; m++) one(i, m, rng);
    } else {
        std::atomic<size_t> next(0); std::vector<std::thread> pool;
        unsigned T = std::max(1u, std::thread::hardware_concurrency());
        for (unsigned t = 0; t < T; t++)
            pool.emplace_back([&] {
                for (size_t i; (i = next++) < parts.size();)
                    for (int m = 0; m < M; m++) { PyRandom rng({seed, (u32)i, (u32)m}); one(i, m, rng); }
            });
        for (auto& th : pool) th.join();
    }
    Stats st; std::map<int, int> hist;
    for (size_t i = 0; i < parts.size(); i++) {
        u32 in1 = 0; for (int u : parts[i]) in1 |= 1u << u;
        std::set<int> have;
        for (const auto& H : hcs) {
            int cnt = 0;
            for (int j = 0; j < n; j++) cnt += (in1 >> H[j] & 1) && (in1 >> H[(j + 1) % n] & 1);
            have.insert(cnt);
        }
        for (int m = 0; m < M; m++) {
            const Out& o = res[i][m]; st.cases++;
            if (!(o.p == o.f && o.f == o.q)) st.mismatch++;
            if (have.count(m)) { if (o.p) { st.nonzero_ok++; hist[m]++; } else st.sz_miss++; }
            else { if (o.p) st.false_pos++; else st.sound_zero++; }
            if (DUMP) {
                std::string v; for (int u : parts[i]) v += std::to_string(u);
                printf("%s %s %d %u %u %u\n", name.c_str(), v.c_str(), m, o.p, o.f, o.q);
            }
        }
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("%s: n=%d, %zu undirected Hamiltonian cycles, %zu partitions, %s, %.1fs\n", name.c_str(), n,
           hcs.size() / 2, parts.size(), use_q ? "p+q+formula" : "p+formula", secs);
    printf("   cases=%ld, sound_zero=%ld, nonzero_ok=%ld, sz_miss=%ld, mismatch=%ld, false_pos=%ld\n",
           st.cases, st.sound_zero, st.nonzero_ok, st.sz_miss, st.mismatch, st.false_pos);
    fflush(stdout);
    Stats& T = TOTAL[parallel];
    T.cases += st.cases; T.sound_zero += st.sound_zero; T.nonzero_ok += st.nonzero_ok;
    T.sz_miss += st.sz_miss; T.mismatch += st.mismatch; T.false_pos += st.false_pos;
    return hist;
}

static Edges randgraph(int n, double p, u32 seed) {     // same graphs as the Python's randgraph
    PyRandom r({seed}); Edges e;
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (r.random() < p) e.push_back({u, v});
    return e;
}
static Edges cycle(int n) { Edges e; for (int i = 0; i < n; i++) e.push_back({i, (i + 1) % n}); return e; }
static long comb(int a, int b) { long r = 1; for (int i = 1; i <= b; i++) r = r * (a - b + i) / i; return r; }
static void cycle_run(int n, bool use_q, bool parallel) {
    // A single Hamiltonian cycle: the test must fire on exactly the partitions counted by the
    // paper's exact direct-arc formula (Appendix A).
    auto hist = run("cycle C_" + std::to_string(n), n, cycle(n), use_q, n, parallel);
    int k = n / 2; long tot = comb(n, k);
    printf("   m  fired on    exact formula (n/r)C(k-1,r-1)^2\n");
    for (int m = 0; m < k; m++) {
        int r = k - m;
        printf("   %d  %4d/%ld    %4ld/%ld\n", m, hist[m], tot, n * comb(k - 1, r - 1) * comb(k - 1, r - 1) / r, tot);
    }
}
static void summary(const char* what, const Stats& T) {
    printf("%s: %ld cases, sound_zero=%ld, nonzero_ok=%ld, sz_miss=%ld, mismatch=%ld, false_pos=%ld\n\n",
           what, T.cases, T.sound_zero, T.nonzero_ok, T.sz_miss, T.mismatch, T.false_pos);
}

int main(int argc, char** argv) {
    for (int i = 1; i < argc; i++) if (!strcmp(argv[i], "--dump")) DUMP = true;
    gf_init();
    printf("Part 1: the graphs of bjorklund_check.py, same random draws\n");
    const struct { int n; double p; u32 seed; const char* name; } R[] = {
        {6, 0.6, 1, "random graph G(n=6, p=0.6)"}, {6, 1.0, 2, "random graph G(n=6, p=1.0)"},
        {8, 0.5, 3, "random graph G(n=8, p=0.5)"}, {8, 0.7, 4, "random graph G(n=8, p=0.7)"}};
    for (auto& g : R) run(g.name, g.n, randgraph(g.n, g.p, g.seed), true, g.seed, false);
    run("hand-made 6-vertex graph", 6, {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {2, 3}, {0, 4}}, true, 9, false);
    run("hand-made 8-vertex graph", 8, {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {2, 6}},
        true, 10, false);
    cycle_run(8, false, false);
    cycle_run(10, false, false);
    summary("Part 1 total", TOTAL[0]);

    printf("Part 2: bigger graphs, all partitions, all three evaluations\n");
    // Petersen graph: the classic non-Hamiltonian graph, and hypohamiltonian (deleting any vertex
    // leaves a Hamiltonian graph), so it is as close to Hamiltonian as a graph can be.
    Edges pet;
    for (int i = 0; i < 5; i++) { pet.push_back({i, (i + 1) % 5}); pet.push_back({i, i + 5}); pet.push_back({5 + i, 5 + (i + 2) % 5}); }
    run("Petersen graph", 10, pet, true, 11, true);
    run("random graph G(n=10, p=0.5)", 10, randgraph(10, 0.5, 12), true, 12, true);
    run("random graph G(n=10, p=0.7)", 10, randgraph(10, 0.7, 13), true, 13, true);
    cycle_run(10, true, true);
    cycle_run(12, true, true);
    run("random graph G(n=12, p=0.5)", 12, randgraph(12, 0.5, 14), true, 14, true);
    summary("Part 2 total", TOTAL[1]);
    return 0;
}
