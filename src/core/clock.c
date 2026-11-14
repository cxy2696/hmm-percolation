#include "hp/core/clock.h"

#include <time.h>

double hp_seconds(void) {
#if defined(CLOCK_MONOTONIC)
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0)
        return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
#endif
    return (double)clock() / (double)CLOCKS_PER_SEC;
}
