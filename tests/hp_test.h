/* Minimal test harness: no external dependency, deterministic order. */
#ifndef HP_TEST_H
#define HP_TEST_H

#include <math.h>
#include <stdio.h>
#include <string.h>

typedef int (*hp_test_fn)(void);

void hp_test_register(const char *name, hp_test_fn fn);
int hp_test_run_all(const char *filter);

#define HP_TEST(name) static int name(void)
#define HP_REG(fn) hp_test_register(#fn, fn)

#define HP_CHECK(cond)                                                      \
    do {                                                                    \
        if (!(cond)) {                                                      \
            fprintf(stderr, "  CHECK failed at %s:%d: %s\n", __FILE__,      \
                    __LINE__, #cond);                                       \
            return 1;                                                       \
        }                                                                   \
    } while (0)

#define HP_CHECK_NEAR(a, b, eps)                                            \
    do {                                                                    \
        double hp__a = (a), hp__b = (b);                                    \
        if (!(fabs(hp__a - hp__b) <= (eps) * (1.0 + fabs(hp__b)))) {        \
            fprintf(stderr,                                                 \
                    "  NEAR failed at %s:%d: %s=%.17g vs %s=%.17g (tol %g)\n", \
                    __FILE__, __LINE__, #a, hp__a, #b, hp__b, (double)(eps)); \
            return 1;                                                       \
        }                                                                   \
    } while (0)

#define HP_CHECK_STREQ(a, b)                                                \
    do {                                                                    \
        const char *hp__a = (a), *hp__b = (b);                              \
        if (!hp__a || !hp__b || strcmp(hp__a, hp__b) != 0) {                \
            fprintf(stderr, "  STREQ failed at %s:%d: \"%s\" vs \"%s\"\n",  \
                    __FILE__, __LINE__, hp__a ? hp__a : "(null)",           \
                    hp__b ? hp__b : "(null)");                              \
            return 1;                                                       \
        }                                                                   \
    } while (0)

#endif /* HP_TEST_H */
