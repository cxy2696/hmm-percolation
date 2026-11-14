#include "hp/core/rng.h"

#include <math.h>
#include <string.h>

#include "hp/core/memory.h"
#include "hp/core/stats.h"

static uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

void hp_rng_seed(hp_rng *rng, uint64_t seed) {
    uint64_t sm = seed;
    for (int i = 0; i < 4; ++i) rng->s[i] = splitmix64(&sm);
}

void hp_rng_seed_stream(hp_rng *rng, uint64_t seed, uint64_t stream) {
    /* Decorrelate streams: two consecutive streams of the same seed must not
     * share the splitmix sequence. */
    uint64_t mixed = seed ^ (stream * 0xD1B54A32D192ED03ULL +
                            0x9E3779B97F4A7C15ULL);
    hp_rng_seed(rng, mixed);
}

uint64_t hp_rng_u64(hp_rng *rng) {
    const uint64_t result = rotl64(rng->s[0] + rng->s[3], 23) + rng->s[0];
    const uint64_t t = rng->s[1] << 17;

    rng->s[2] ^= rng->s[0];
    rng->s[3] ^= rng->s[1];
    rng->s[1] ^= rng->s[2];
    rng->s[0] ^= rng->s[3];
    rng->s[2] ^= t;
    rng->s[3] = rotl64(rng->s[3], 45);
    return result;
}

uint32_t hp_rng_u32(hp_rng *rng) { return (uint32_t)(hp_rng_u64(rng) >> 32); }

double hp_rng_uniform(hp_rng *rng) {
    return (double)(hp_rng_u64(rng) >> 11) * 0x1.0p-53;
}

double hp_rng_uniform_range(hp_rng *rng, double lo, double hi) {
    return lo + (hi - lo) * hp_rng_uniform(rng);
}

double hp_rng_normal(hp_rng *rng) {
    /* Marsaglia polar method: no cached spare, so the stream position is a
     * pure function of the number of calls. */
    for (;;) {
        double u = 2.0 * hp_rng_uniform(rng) - 1.0;
        double v = 2.0 * hp_rng_uniform(rng) - 1.0;
        double s = u * u + v * v;
        if (s > 0.0 && s < 1.0) return u * sqrt(-2.0 * log(s) / s);
    }
}

double hp_rng_exponential(hp_rng *rng) {
    double u = hp_rng_uniform(rng);
    while (u <= 0.0) u = hp_rng_uniform(rng);
    return -log(u);
}

double hp_rng_gamma(hp_rng *rng, double shape) {
    if (!(shape > 0.0) || !isfinite(shape)) return NAN;
    if (shape < 1.0) {
        double u = hp_rng_uniform(rng);
        while (u <= 0.0) u = hp_rng_uniform(rng);
        return hp_rng_gamma(rng, shape + 1.0) * pow(u, 1.0 / shape);
    }
    const double d = shape - 1.0 / 3.0;
    const double c = 1.0 / sqrt(9.0 * d);
    for (;;) {
        double x, v;
        do {
            x = hp_rng_normal(rng);
            v = 1.0 + c * x;
        } while (v <= 0.0);
        v = v * v * v;
        double u = hp_rng_uniform(rng);
        if (u <= 0.0) u = 0x1.0p-53;
        if (log(u) < 0.5 * x * x + d - d * v + d * log(v)) return d * v;
    }
}

double hp_rng_student_t(hp_rng *rng, double df) {
    if (!(df > 0.0)) return NAN;
    double z = hp_rng_normal(rng);
    double chi2 = 2.0 * hp_rng_gamma(rng, 0.5 * df);
    return z / sqrt(chi2 / df);
}

size_t hp_rng_below(hp_rng *rng, size_t n) {
    if (n <= 1) return 0;
    uint64_t bound = (uint64_t)n;
    uint64_t limit = UINT64_MAX - (UINT64_MAX % bound);
    for (;;) {
        uint64_t r = hp_rng_u64(rng);
        if (r < limit) return (size_t)(r % bound);
    }
}

void hp_rng_dirichlet(hp_rng *rng, const double *alpha, size_t n, double *out) {
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double a = alpha ? alpha[i] : 1.0;
        out[i] = hp_rng_gamma(rng, a);
        if (!(out[i] >= 0.0) || !isfinite(out[i])) out[i] = 0.0;
        sum += out[i];
    }
    if (!(sum > 0.0)) {
        for (size_t i = 0; i < n; ++i) out[i] = 1.0 / (double)n;
        return;
    }
    for (size_t i = 0; i < n; ++i) out[i] /= sum;
}

void hp_rng_shuffle(hp_rng *rng, void *base, size_t n, size_t size) {
    unsigned char *p = (unsigned char *)base;
    for (size_t i = n; i > 1; --i) {
        size_t j = hp_rng_below(rng, i);
        unsigned char *a = p + (i - 1) * size;
        unsigned char *b = p + j * size;
        for (size_t k = 0; k < size; ++k) {
            unsigned char t = a[k];
            a[k] = b[k];
            b[k] = t;
        }
    }
}
