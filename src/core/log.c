#include "hp/core/log.h"

#include <stdarg.h>
#include <stdio.h>

static hp_log_level g_level = HP_LOG_WARN;

void hp_log_set_level(hp_log_level level) { g_level = level; }
hp_log_level hp_log_get_level(void) { return g_level; }

void hp_log_message(hp_log_level level, const char *file, int line,
                    const char *fmt, ...) {
    if (level < g_level || g_level == HP_LOG_OFF) return;
    static const char *tag[] = {"debug", "info", "warn", "error"};
    const char *name = (level >= 0 && level <= HP_LOG_ERROR) ? tag[level] : "?";
    fprintf(stderr, "[hp %s] ", name);
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    if (file && level >= HP_LOG_WARN) fprintf(stderr, " (%s:%d)", file, line);
    fputc('\n', stderr);
}
