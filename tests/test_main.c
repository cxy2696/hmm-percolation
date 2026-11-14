#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "hp_test.h"

#define MAX_TESTS 256

typedef struct {
    const char *name;
    hp_test_fn fn;
} test_entry;

static test_entry g_tests[MAX_TESTS];
static int g_n = 0;

void hp_test_register(const char *name, hp_test_fn fn) {
    if (g_n >= MAX_TESTS || !name || !fn) return;
    g_tests[g_n].name = name;
    g_tests[g_n].fn = fn;
    ++g_n;
}

int hp_test_run_all(const char *filter) {
    int passed = 0, failed = 0;
    for (int i = 0; i < g_n; ++i) {
        if (filter && *filter && !strstr(g_tests[i].name, filter)) continue;
        printf("[ RUN  ] %s\n", g_tests[i].name);
        fflush(stdout);
        int rc = g_tests[i].fn();
        if (rc == 0) {
            printf("[  OK  ] %s\n", g_tests[i].name);
            ++passed;
        } else {
            printf("[ FAIL ] %s\n", g_tests[i].name);
            ++failed;
        }
    }
    printf("\n%d passed, %d failed\n", passed, failed);
    return failed == 0 ? 0 : 1;
}

void hp_test_core_register(void);
void hp_test_hmm_register(void);
void hp_test_ghmm_register(void);
void hp_test_perc_register(void);
void hp_test_market_register(void);
void hp_test_strategy_register(void);

int main(int argc, char **argv) {
    mkdir("build", 0755);
    mkdir("build/test_tmp", 0755);
    hp_test_core_register();
    hp_test_hmm_register();
    hp_test_ghmm_register();
    hp_test_perc_register();
    hp_test_market_register();
    hp_test_strategy_register();
    return hp_test_run_all(argc > 1 ? argv[1] : NULL);
}
