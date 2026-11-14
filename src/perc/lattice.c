#include "hp/perc/lattice.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "hp/core/memory.h"
#include "hp/perc/uf.h"

static size_t n_sites(const hp_lattice *lat) {
    return (size_t)lat->L * (size_t)lat->L;
}

hp_lattice *hp_lattice_new(int L, hp_perc_kind kind, int periodic) {
    if (L <= 0) return NULL;
    if (kind != HP_PERC_SITE && kind != HP_PERC_BOND) return NULL;
    size_t n = (size_t)L * (size_t)L;
    hp_lattice *lat = (hp_lattice *)hp_xcalloc(1, sizeof(*lat));
    lat->L = L;
    lat->periodic = periodic ? 1 : 0;
    lat->kind = kind;
    lat->site = hp_u8s(n);
    lat->h = hp_u8s(n);
    lat->v = hp_u8s(n);
    return lat;
}

void hp_lattice_free(hp_lattice *lat) {
    if (!lat) return;
    hp_free(lat->site);
    hp_free(lat->h);
    hp_free(lat->v);
    hp_free(lat);
}

void hp_lattice_clear(hp_lattice *lat) {
    if (!lat) return;
    size_t n = n_sites(lat);
    memset(lat->site, 0, n);
    memset(lat->h, 0, n);
    memset(lat->v, 0, n);
}

hp_status hp_lattice_bernoulli(hp_lattice *lat, double p, hp_rng *rng) {
    if (!lat || !rng || !(p >= 0.0 && p <= 1.0)) return HP_ERR_INVALID;
    size_t n = n_sites(lat);
    if (lat->kind == HP_PERC_BOND) {
        memset(lat->site, 1, n);
        memset(lat->h, 0, n);
        memset(lat->v, 0, n);
        for (size_t i = 0; i < n; ++i) {
            lat->h[i] = hp_rng_uniform(rng) < p ? 1 : 0;
            lat->v[i] = hp_rng_uniform(rng) < p ? 1 : 0;
        }
    } else {
        memset(lat->h, 0, n);
        memset(lat->v, 0, n);
        for (size_t i = 0; i < n; ++i)
            lat->site[i] = hp_rng_uniform(rng) < p ? 1 : 0;
    }
    return HP_OK;
}

int hp_lattice_site_open(const hp_lattice *lat, int i, int j) {
    if (!lat || i < 0 || j < 0 || i >= lat->L || j >= lat->L) return 0;
    if (lat->kind == HP_PERC_BOND) return 1;
    return lat->site[(size_t)i * (size_t)lat->L + (size_t)j] ? 1 : 0;
}

/* Builds the occupied cluster forest.  active[] is filled for every site. */
static void lattice_union(const hp_lattice *lat, hp_uf *uf, uint8_t *active) {
    const int L = lat->L;
    const size_t N = (size_t)L * (size_t)L;
    if (lat->kind == HP_PERC_BOND) {
        memset(active, 1, N);
    } else {
        memcpy(active, lat->site, N);
    }
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            size_t u = (size_t)i * (size_t)L + (size_t)j;
            if (!active[u]) continue;
            /* right */
            int jr = j + 1;
            if (jr == L) {
                if (lat->periodic) jr = 0;
                else jr = -1;
            }
            if (jr >= 0) {
                size_t v = (size_t)i * (size_t)L + (size_t)jr;
                int edge = lat->kind == HP_PERC_SITE ? 1 : lat->h[u];
                if (active[v] && edge) hp_uf_union(uf, u, v);
            }
            /* down */
            int id = i + 1;
            if (id == L) {
                if (lat->periodic) id = 0;
                else id = -1;
            }
            if (id >= 0) {
                size_t v = (size_t)id * (size_t)L + (size_t)j;
                int edge = lat->kind == HP_PERC_SITE ? 1 : lat->v[u];
                if (active[v] && edge) hp_uf_union(uf, u, v);
            }
        }
    }
}

hp_status hp_lattice_analyze(const hp_lattice *lat, hp_lattice_stats *out) {
    if (!lat || !out) return HP_ERR_INVALID;
    memset(out, 0, sizeof(*out));
    const int L = lat->L;
    const size_t N = n_sites(lat);
    hp_uf *uf = hp_uf_new(N);
    uint8_t *active = hp_u8s(N);
    lattice_union(lat, uf, active);

    size_t *cnt = hp_sizes(N);
    size_t occupied = 0;
    for (size_t i = 0; i < N; ++i) {
        if (!active[i]) continue;
        ++occupied;
        cnt[hp_uf_find(uf, i)]++;
    }

    uint8_t *tl = hp_u8s(N), *tr = hp_u8s(N), *tt = hp_u8s(N), *tb = hp_u8s(N);
    for (int k = 0; k < L; ++k) {
        size_t left = (size_t)k * (size_t)L;
        size_t right = left + (size_t)(L - 1);
        size_t top = (size_t)k;
        size_t bottom = (size_t)(L - 1) * (size_t)L + (size_t)k;
        if (active[left]) tl[hp_uf_find(uf, left)] = 1;
        if (active[right]) tr[hp_uf_find(uf, right)] = 1;
        if (active[top]) tt[hp_uf_find(uf, top)] = 1;
        if (active[bottom]) tb[hp_uf_find(uf, bottom)] = 1;
    }

    size_t max_size = 0, n_clusters = 0;
    double sumsq = 0.0;
    for (size_t i = 0; i < N; ++i) {
        if (cnt[i] == 0) continue;
        ++n_clusters;
        if (cnt[i] > max_size) max_size = cnt[i];
        sumsq += (double)cnt[i] * (double)cnt[i];
        if (tl[i] && tr[i]) out->spanning_x = 1;
        if (tt[i] && tb[i]) out->spanning_y = 1;
    }
    out->spanning = out->spanning_x || out->spanning_y;
    out->n_sites = N;
    out->n_occupied = occupied;
    out->n_clusters = n_clusters;
    out->max_cluster = max_size;
    out->s_max = (double)max_size / (double)N;
    out->occupancy = (double)occupied / (double)N;
    out->susceptibility = (sumsq - (double)max_size * (double)max_size) / (double)N;

    hp_free(cnt);
    hp_free(tl);
    hp_free(tr);
    hp_free(tt);
    hp_free(tb);
    hp_free(active);
    hp_uf_free(uf);
    return HP_OK;
}

static int cmp_size_desc(const void *a, const void *b) {
    size_t x = *(const size_t *)a, y = *(const size_t *)b;
    return x > y ? -1 : (x < y ? 1 : 0);
}

size_t *hp_lattice_cluster_sizes(const hp_lattice *lat, size_t *n_out) {
    if (n_out) *n_out = 0;
    if (!lat) return NULL;
    const size_t N = n_sites(lat);
    hp_uf *uf = hp_uf_new(N);
    uint8_t *active = hp_u8s(N);
    lattice_union(lat, uf, active);
    size_t *cnt = hp_sizes(N);
    for (size_t i = 0; i < N; ++i)
        if (active[i]) cnt[hp_uf_find(uf, i)]++;
    size_t nc = 0;
    for (size_t i = 0; i < N; ++i)
        if (cnt[i] > 0) ++nc;
    size_t *sizes = hp_sizes(nc ? nc : 1);
    size_t k = 0;
    for (size_t i = 0; i < N; ++i)
        if (cnt[i] > 0) sizes[k++] = cnt[i];
    qsort(sizes, nc, sizeof(size_t), cmp_size_desc);
    if (n_out) *n_out = nc;
    hp_free(cnt);
    hp_free(active);
    hp_uf_free(uf);
    return sizes;
}

const char *hp_perc_kind_name(hp_perc_kind kind) {
    return kind == HP_PERC_BOND ? "bond" : "site";
}

static double spanning_prob_seed(const hp_perc_exp *e, double p) {
    hp_rng rng;
    uint64_t stream = (uint64_t)(unsigned)e->L * 1000003u +
                      (uint64_t)e->kind * 17u + (uint64_t)e->periodic;
    hp_rng_seed_stream(&rng, e->seed, stream);
    hp_lattice *lat = hp_lattice_new(e->L, e->kind, e->periodic);
    size_t hits = 0;
    for (size_t r = 0; r < e->reps; ++r) {
        hp_lattice_bernoulli(lat, p, &rng);
        hp_lattice_stats st;
        hp_lattice_analyze(lat, &st);
        if (st.spanning) ++hits;
    }
    hp_lattice_free(lat);
    return e->reps ? (double)hits / (double)e->reps : NAN;
}

double hp_perc_spanning_prob(const hp_perc_exp *exp, double p) {
    if (!exp || exp->L <= 0 || exp->reps == 0 || !(p >= 0.0 && p <= 1.0))
        return NAN;
    return spanning_prob_seed(exp, p);
}

hp_status hp_perc_crossing(const hp_perc_exp *exp, double p_lo, double p_hi,
                           unsigned iters, double *p_cross,
                           double *span_prob) {
    if (!exp || exp->L <= 0 || exp->reps == 0 || !(p_lo >= 0.0) ||
        !(p_hi <= 1.0) || !(p_lo < p_hi))
        return HP_ERR_INVALID;
    double lo = p_lo, hi = p_hi;
    double last = NAN;
    for (unsigned i = 0; i < iters; ++i) {
        double mid = 0.5 * (lo + hi);
        last = spanning_prob_seed(exp, mid);
        if (last >= 0.5)
            hi = mid;
        else
            lo = mid;
    }
    if (p_cross) *p_cross = 0.5 * (lo + hi);
    if (span_prob) *span_prob = last;
    return HP_OK;
}
