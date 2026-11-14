#include <math.h>
#include <stdio.h>
#include <string.h>

#include "hp/core/csv.h"
#include "hp/core/io.h"
#include "hp/core/memory.h"
#include "hp/core/rng.h"
#include "hp/core/stats.h"
#include "hp_test.h"

HP_TEST(test_rng_determinism) {
    hp_rng a, b, c;
    hp_rng_seed(&a, 12345);
    hp_rng_seed(&b, 12345);
    hp_rng_seed_stream(&c, 12345, 1);
    for (int i = 0; i < 1000; ++i) {
        uint64_t x = hp_rng_u64(&a);
        uint64_t y = hp_rng_u64(&b);
        HP_CHECK(x == y);
    }
    hp_rng_seed(&a, 12345);
    hp_rng_seed_stream(&b, 12345, 1);
    int different = 0;
    for (int i = 0; i < 100; ++i)
        if (hp_rng_u64(&a) != hp_rng_u64(&b)) different = 1;
    HP_CHECK(different);
    return 0;
}

HP_TEST(test_rng_moments) {
    hp_rng rng;
    hp_rng_seed(&rng, 7);
    const int n = 200000;
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        double u = hp_rng_uniform(&rng);
        HP_CHECK(u >= 0.0 && u < 1.0);
        s += u;
    }
    HP_CHECK_NEAR(s / n, 0.5, 0.01);

    double mn = 0.0, mn2 = 0.0;
    for (int i = 0; i < n; ++i) {
        double z = hp_rng_normal(&rng);
        mn += z;
        mn2 += z * z;
    }
    mn /= n;
    mn2 = mn2 / n - mn * mn;
    HP_CHECK_NEAR(mn, 0.0, 0.01);
    HP_CHECK_NEAR(mn2, 1.0, 0.01);

    double mg = 0.0;
    for (int i = 0; i < n; ++i) mg += hp_rng_gamma(&rng, 2.0);
    HP_CHECK_NEAR(mg / n, 2.0, 0.02);

    double alpha[3] = {1.0, 2.0, 3.0}, w[3];
    hp_rng_dirichlet(&rng, alpha, 3, w);
    HP_CHECK_NEAR(w[0] + w[1] + w[2], 1.0, 1e-12);
    HP_CHECK(w[0] > 0.0 && w[1] > 0.0 && w[2] > 0.0);

    double st = 0.0, st2 = 0.0;
    for (int i = 0; i < n; ++i) {
        double t = hp_rng_student_t(&rng, 5.0);
        st += t;
        st2 += t * t;
    }
    st /= n;
    st2 = st2 / n - st * st;
    HP_CHECK_NEAR(st2, 5.0 / 3.0, 0.06);
    return 0;
}

HP_TEST(test_stats_basics) {
    double x[4] = {1, 2, 3, 4};
    HP_CHECK_NEAR(hp_mean(x, 4), 2.5, 1e-15);
    HP_CHECK_NEAR(hp_var(x, 4), 5.0 / 3.0, 1e-15);
    HP_CHECK_NEAR(hp_quantile(x, 4, 0.5), 2.5, 1e-15);
    HP_CHECK_NEAR(hp_quantile(x, 4, 0.25), 1.75, 1e-15);

    double y[4] = {2, 4, 6, 8};
    HP_CHECK_NEAR(hp_corr_complete(x, y, 4), 1.0, 1e-15);
    double z[4] = {4, 3, 2, 1};
    HP_CHECK_NEAR(hp_corr_complete(x, z, 4), -1.0, 1e-15);
    double c[4] = {2, 2, 2, 2};
    HP_CHECK_NEAR(hp_corr_complete(x, c, 4), 0.0, 1e-15);

    double data[8] = {1, 2, 2, 4, 3, 6, 4, 8}; /* 4 rows x 2 cols */
    double corr[4];
    hp_corr_matrix(data, 4, 2, corr);
    HP_CHECK_NEAR(corr[0], 1.0, 1e-15);
    HP_CHECK_NEAR(corr[1], 1.0, 1e-15);
    HP_CHECK_NEAR(corr[2], 1.0, 1e-15);
    HP_CHECK_NEAR(corr[3], 1.0, 1e-15);

    double v[3] = {1000.0, 1001.0, 1002.0};
    HP_CHECK_NEAR(hp_logsumexp(v, 3),
                  log(exp(v[0]) + exp(v[1]) + exp(v[2])), 1e-12);
    double ninf[2] = {-INFINITY, -INFINITY};
    HP_CHECK(hp_logsumexp(ninf, 2) == -INFINITY);

    double lx[3] = {1, 2, 3}, ly[3] = {5, 7, 9};
    double slope = 0, intercept = 0;
    double r2 = hp_linreg(lx, ly, 3, &slope, &intercept);
    HP_CHECK_NEAR(slope, 2.0, 1e-12);
    HP_CHECK_NEAR(intercept, 3.0, 1e-12);
    HP_CHECK_NEAR(r2, 1.0, 1e-12);

    double sorted[6] = {1, 3, 5, 7, 9, 0};
    hp_sorted_insert(sorted, 5, 6);
    HP_CHECK_NEAR(sorted[3], 6.0, 1e-15);
    HP_CHECK_NEAR(sorted[4], 7.0, 1e-15);
    HP_CHECK_NEAR(sorted[5], 9.0, 1e-15);
    hp_sorted_remove(sorted, 6, 6);
    HP_CHECK_NEAR(sorted[3], 7.0, 1e-15);
    HP_CHECK_NEAR(sorted[5], 9.0, 1e-15);
    return 0;
}

HP_TEST(test_csv_roundtrip) {
    const char *path = "build/test_tmp/core_roundtrip.csv";
    const char *names[3] = {"a", "b", "c"};
    double data[9] = {1.0, NAN, 3.0, 4.0, 5.5, 6.0, -7.5, 0.0, 1e300};
    HP_CHECK(hp_csv_write_numeric(path, names, 3, data, 3) == HP_OK);

    hp_csv csv;
    HP_CHECK(hp_csv_read(path, &csv) == HP_OK);
    HP_CHECK(csv.n_rows == 3);
    HP_CHECK(csv.n_cols == 3);
    HP_CHECK_STREQ(csv.header[1], "b");
    double v = 0;
    HP_CHECK(hp_csv_parse_double(hp_csv_cell(&csv, 0, 0), &v) == HP_OK);
    HP_CHECK_NEAR(v, 1.0, 1e-15);
    HP_CHECK(hp_csv_parse_double(hp_csv_cell(&csv, 0, 1), &v) == HP_OK);
    HP_CHECK(isnan(v));
    HP_CHECK(hp_csv_parse_double(hp_csv_cell(&csv, 2, 2), &v) == HP_OK);
    HP_CHECK_NEAR(v, 1e300, 1e-15);
    HP_CHECK(hp_csv_col_index(&csv, "C") == 2);
    HP_CHECK(hp_csv_col_index(&csv, "zz") == -1);
    hp_csv_free(&csv);
    return 0;
}

HP_TEST(test_csv_label_and_errors) {
    const char *path = "build/test_tmp/core_label.csv";
    const char *names[3] = {"date", "x", "y"};
    hp_csv_writer *w = hp_csv_writer_open(path, names, 3, 1);
    HP_CHECK(w != NULL);
    double row1[2] = {1.5, 2.5};
    double row2[2] = {-1.0, 0.0};
    HP_CHECK(hp_csv_writer_row_label(w, "2020-01-01", row1, 2) == HP_OK);
    HP_CHECK(hp_csv_writer_row_label(w, "2020-01-02", row2, 2) == HP_OK);
    HP_CHECK(hp_csv_writer_close(w) == HP_OK);

    hp_csv csv;
    HP_CHECK(hp_csv_read(path, &csv) == HP_OK);
    HP_CHECK(csv.n_rows == 2 && csv.n_cols == 3);
    HP_CHECK_STREQ(hp_csv_cell(&csv, 0, 0), "2020-01-01");
    double v = 0;
    HP_CHECK(hp_csv_parse_double(hp_csv_cell(&csv, 1, 1), &v) == HP_OK);
    HP_CHECK_NEAR(v, -1.0, 1e-15);
    hp_csv_free(&csv);

    const char *bad = "build/test_tmp/core_bad.csv";
    HP_CHECK(hp_write_text(bad, "a,b\n1,2,3\n") == HP_OK);
    hp_csv bad_csv;
    HP_CHECK(hp_csv_read(bad, &bad_csv) == HP_ERR_FORMAT);

    HP_CHECK(hp_write_text(bad, "a,b\n1\n") == HP_OK);
    HP_CHECK(hp_csv_read(bad, &bad_csv) == HP_OK);
    HP_CHECK(bad_csv.n_rows == 1 && bad_csv.n_cols == 2);
    HP_CHECK(hp_csv_parse_double(hp_csv_cell(&bad_csv, 0, 1), &v) == HP_OK);
    HP_CHECK(isnan(v));
    hp_csv_free(&bad_csv);
    return 0;
}

HP_TEST(test_io_helpers) {
    char s[32];
    strcpy(s, "  hello \t");
    HP_CHECK_STREQ(hp_trim(s), "hello");
    HP_CHECK(hp_str_ieq("NA", "na"));
    HP_CHECK(!hp_str_ieq("na", "nan"));
    HP_CHECK_STREQ(hp_status_name(HP_ERR_FORMAT), "HP_ERR_FORMAT");
    double v = 0;
    HP_CHECK(hp_csv_parse_double("1x", &v) == HP_ERR_FORMAT);
    HP_CHECK(hp_csv_parse_double("", &v) == HP_OK && isnan(v));
    return 0;
}

void hp_test_core_register(void) {
    HP_REG(test_rng_determinism);
    HP_REG(test_rng_moments);
    HP_REG(test_stats_basics);
    HP_REG(test_csv_roundtrip);
    HP_REG(test_csv_label_and_errors);
    HP_REG(test_io_helpers);
}
