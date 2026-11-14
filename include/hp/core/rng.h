/* HP · deterministic random numbers.
 *
 * Generator: xoshiro256++ (Blackman & Vigna, public domain).  Seeding goes
 * through splitmix64.  Two runs with the same (seed, stream) produce exactly
 * the same sequence on every platform: all draws are built only from 64-bit
 * integer arithmetic and IEEE-754 double operations.
 *
 * Every stochastic experiment in the project takes a seed and derives
 * independent streams with hp_rng_seed_stream(): never share one generator
 * between two logically independent experiments, or adding a draw in one
 * place will silently change the other.
 */
#ifndef HP_CORE_RNG_H
#define HP_CORE_RNG_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint64_t s[4];
} hp_rng;

void hp_rng_seed(hp_rng *rng, uint64_t seed);
/* Independent stream: seed' = splitmix64(seed ^ mix(stream)). */
void hp_rng_seed_stream(hp_rng *rng, uint64_t seed, uint64_t stream);

uint64_t hp_rng_u64(hp_rng *rng);
uint32_t hp_rng_u32(hp_rng *rng);
/* Uniform in [0, 1) with 53 random bits. */
double hp_rng_uniform(hp_rng *rng);
double hp_rng_uniform_range(hp_rng *rng, double lo, double hi);
double hp_rng_normal(hp_rng *rng);      /* N(0,1), Marsaglia polar method */
double hp_rng_exponential(hp_rng *rng); /* rate 1 */
double hp_rng_gamma(hp_rng *rng, double shape); /* scale 1, shape > 0 */
double hp_rng_student_t(hp_rng *rng, double df);
/* Uniform integer in [0, n), unbiased by rejection. */
size_t hp_rng_below(hp_rng *rng, size_t n);

/* Draw a Dirichlet(alpha[0..n-1]) vector; out must have n entries. */
void hp_rng_dirichlet(hp_rng *rng, const double *alpha, size_t n, double *out);
/* In-place Fisher-Yates over an array of n elements of `size` bytes. */
void hp_rng_shuffle(hp_rng *rng, void *base, size_t n, size_t size);

#ifdef __cplusplus
}
#endif

#endif /* HP_CORE_RNG_H */
