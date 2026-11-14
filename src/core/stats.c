#include "hp/core/stats.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "hp/core/memory.h"

int hp_is_finite(double x) { return isfinite(x) ? 1 : 0; }

double hp_sum(const double *x, size_t n) {
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) s += x[i];
    return s;
}

double hp_mean(const double *x, size_t n) {
    if (n == 0) return NAN;
    return hp_sum(x, n) / (double)n;
}

double hp_var(const double *x, size_t n) {
    if (n < 2) return NAN;
    double m = hp_mean(x, n);
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double d = x[i] - m;
        s += d * d;
    }
    return s / (double)(n - 1);
}

double hp_sd(const double *x, size_t n) { return sqrt(hp_var(x, n)); }

double hp_min(const double *x, size_t n) {
    if (n == 0) return NAN;
    double v = x[0];
    for (size_t i = 1; i < n; ++i)
        if (x[i] < v) v = x[i];
    return v;
}

double hp_max(const double *x, size_t n) {
    if (n == 0) return NAN;
    double v = x[0];
    for (size_t i = 1; i < n; ++i)
        if (x[i] > v) v = x[i];
    return v;
}

double hp_dot(const double *a, const double *b, size_t n) {
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

double hp_cov_complete(const double *x, const double *y, size_t n) {
    size_t cnt = 0;
    double mx = 0.0, my = 0.0;
    for (size_t i = 0; i < n; ++i) {
        if (!isfinite(x[i]) || !isfinite(y[i])) continue;
        ++cnt;
        mx += x[i];
        my += y[i];
    }
    if (cnt < 2) return NAN;
    mx /= (double)cnt;
    my /= (double)cnt;
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) {
        if (!isfinite(x[i]) || !isfinite(y[i])) continue;
        s += (x[i] - mx) * (y[i] - my);
    }
    return s / (double)(cnt - 1);
}

double hp_corr_complete(const double *x, const double *y, size_t n) {
    size_t cnt = 0;
    double mx = 0.0, my = 0.0;
    for (size_t i = 0; i < n; ++i) {
        if (!isfinite(x[i]) || !isfinite(y[i])) continue;
        ++cnt;
        mx += x[i];
        my += y[i];
    }
    if (cnt < 2) return NAN;
    mx /= (double)cnt;
    my /= (double)cnt;
    double sxy = 0.0, sxx = 0.0, syy = 0.0;
    for (size_t i = 0; i < n; ++i) {
        if (!isfinite(x[i]) || !isfinite(y[i])) continue;
        double dx = x[i] - mx;
        double dy = y[i] - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    if (!(sxx > 0.0) || !(syy > 0.0)) return 0.0; /* constant column */
    double r = sxy / sqrt(sxx * syy);
    return hp_clamp(r, -1.0, 1.0);
}

void hp_corr_matrix(const double *x, size_t T, size_t N, double *corr) {
    for (size_t i = 0; i < N; ++i) {
        size_t cnt = 0;
        for (size_t t = 0; t < T; ++t)
            if (isfinite(x[t * N + i])) ++cnt;
        corr[i * N + i] = cnt >= 2 ? 1.0 : 0.0;
        for (size_t j = i + 1; j < N; ++j) {
            size_t c = 0;
            double mx = 0.0, my = 0.0;
            for (size_t t = 0; t < T; ++t) {
                double a = x[t * N + i], b = x[t * N + j];
                if (!isfinite(a) || !isfinite(b)) continue;
                ++c;
                mx += a;
                my += b;
            }
            double r = 0.0;
            if (c >= 2) {
                mx /= (double)c;
                my /= (double)c;
                double sxy = 0.0, sxx = 0.0, syy = 0.0;
                for (size_t t = 0; t < T; ++t) {
                    double a = x[t * N + i], b = x[t * N + j];
                    if (!isfinite(a) || !isfinite(b)) continue;
                    double da = a - mx, db = b - my;
                    sxy += da * db;
                    sxx += da * da;
                    syy += db * db;
                }
                r = (sxx > 0.0 && syy > 0.0) ? sxy / sqrt(sxx * syy) : 0.0;
                r = hp_clamp(r, -1.0, 1.0);
            }
            corr[i * N + j] = r;
            corr[j * N + i] = r;
        }
    }
}

static int cmp_double(const void *pa, const void *pb) {
    double a = *(const double *)pa, b = *(const double *)pb;
    return (a > b) - (a < b);
}

void hp_sort(double *x, size_t n) {
    if (n > 1) qsort(x, n, sizeof(double), cmp_double);
}

double hp_quantile_sorted(const double *sorted, size_t n, double q) {
    if (n == 0) return NAN;
    if (q <= 0.0) return sorted[0];
    if (q >= 1.0) return sorted[n - 1];
    double pos = q * (double)(n - 1);
    size_t lo = (size_t)floor(pos);
    size_t hi = lo + 1 < n ? lo + 1 : lo;
    double frac = pos - (double)lo;
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

double hp_quantile(const double *x, size_t n, double q) {
    double *tmp = hp_doubles(n);
    size_t m = 0;
    for (size_t i = 0; i < n; ++i)
        if (isfinite(x[i])) tmp[m++] = x[i];
    double out = NAN;
    if (m) {
        hp_sort(tmp, m);
        out = hp_quantile_sorted(tmp, m, q);
    }
    hp_free(tmp);
    return out;
}

double hp_logsumexp(const double *v, size_t n) {
    if (n == 0) return -INFINITY;
    double mx = -INFINITY;
    for (size_t i = 0; i < n; ++i) {
        if (isnan(v[i])) return NAN;
        if (v[i] > mx) mx = v[i];
    }
    if (mx == -INFINITY) return -INFINITY;
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) s += exp(v[i] - mx);
    return mx + log(s);
}

double hp_logsumexp_finite(const double *v, size_t n) {
    if (n == 0) return -INFINITY;
    double mx = v[0];
    for (size_t i = 1; i < n; ++i)
        if (v[i] > mx) mx = v[i];
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) s += exp(v[i] - mx);
    return mx + log(s);
}

double hp_norm_logpdf(double x, double mean, double var) {
    if (!(var > 0.0) || !isfinite(var) || !isfinite(x) || !isfinite(mean))
        return -INFINITY;
    double d = x - mean;
    return -0.5 * (HP_LOG_2PI + log(var) + d * d / var);
}

double hp_norm_pdf(double x, double mean, double var) {
    double lp = hp_norm_logpdf(x, mean, var);
    return lp == -INFINITY ? 0.0 : exp(lp);
}

double hp_clamp(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

double hp_linreg(const double *x, const double *y, size_t n, double *slope,
                 double *intercept) {
    if (n < 2) {
        if (slope) *slope = NAN;
        if (intercept) *intercept = NAN;
        return NAN;
    }
    double mx = hp_mean(x, n), my = hp_mean(y, n);
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double dx = x[i] - mx, dy = y[i] - my;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    if (!(sxx > 0.0)) {
        if (slope) *slope = 0.0;
        if (intercept) *intercept = my;
        return 0.0;
    }
    double b = sxy / sxx;
    if (slope) *slope = b;
    if (intercept) *intercept = my - b * mx;
    if (!(syy > 0.0)) return 0.0;
    double r2 = (sxy * sxy) / (sxx * syy);
    return hp_clamp(r2, 0.0, 1.0);
}

void hp_sorted_insert(double *sorted, size_t n, double value) {
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (sorted[mid] < value)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo < n) memmove(sorted + lo + 1, sorted + lo, (n - lo) * sizeof(double));
    sorted[lo] = value;
}

void hp_sorted_remove(double *sorted, size_t n, double value) {
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (sorted[mid] < value)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo < n && sorted[lo] == value) {
        if (lo + 1 < n)
            memmove(sorted + lo, sorted + lo + 1, (n - lo - 1) * sizeof(double));
    }
}
