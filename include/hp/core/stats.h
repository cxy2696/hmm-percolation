/* HP · small statistical toolkit (no external dependencies). */
#ifndef HP_CORE_STATS_H
#define HP_CORE_STATS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define HP_PI 3.14159265358979323846
#define HP_LOG_2PI 1.83787706640934548356

int hp_is_finite(double x);

double hp_sum(const double *x, size_t n);
double hp_mean(const double *x, size_t n);
/* Unbiased sample variance (n-1).  Returns NAN for n < 2. */
double hp_var(const double *x, size_t n);
double hp_sd(const double *x, size_t n);
double hp_min(const double *x, size_t n);
double hp_max(const double *x, size_t n);
double hp_dot(const double *a, const double *b, size_t n);

/* Sample covariance / correlation; pairs where either value is non-finite
 * are skipped.  Return NAN when fewer than 2 complete pairs exist. */
double hp_cov_complete(const double *x, const double *y, size_t n);
double hp_corr_complete(const double *x, const double *y, size_t n);

/* Pairwise-complete correlation matrix of a T x N column-major-free layout
 * (row t, column j at x[t*N + j]).  corr is N*N row-major; diagonal = 1 when
 * the column has >= 2 finite values, otherwise diagonal = 0. */
void hp_corr_matrix(const double *x, size_t T, size_t N, double *corr);

/* Quantile with linear interpolation; `q` in [0,1].  hp_quantile makes a
 * sorted copy, hp_quantile_sorted consumes an already sorted array. */
double hp_quantile(const double *x, size_t n, double q);
double hp_quantile_sorted(const double *sorted, size_t n, double q);
void hp_sort(double *x, size_t n);

/* log(sum exp(v[i])) in a numerically stable way.  Returns -INFINITY when
 * every entry is -INFINITY, and NAN if any entry is NAN. */
double hp_logsumexp(const double *v, size_t n);
/* log(sum exp(v[i])) where all entries are finite. */
double hp_logsumexp_finite(const double *v, size_t n);

double hp_norm_logpdf(double x, double mean, double var);
double hp_norm_pdf(double x, double mean, double var);
double hp_clamp(double x, double lo, double hi);

/* Ordinary least squares y = a + b x.  Returns R^2 (0 if undefined). */
double hp_linreg(const double *x, const double *y, size_t n, double *slope,
                 double *intercept);

/* Expanding-window percentile helpers used by the app layer.  `sorted` must
 * hold the last `n` values in ascending order; insertion helpers keep it
 * sorted without re-sorting the whole window. */
void hp_sorted_insert(double *sorted, size_t n, double value);
void hp_sorted_remove(double *sorted, size_t n, double value);

#ifdef __cplusplus
}
#endif

#endif /* HP_CORE_STATS_H */
