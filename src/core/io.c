#include "hp/core/io.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hp/core/memory.h"

hp_status hp_read_file(const char *path, char **out, size_t *out_len) {
    if (!path || !out) return HP_ERR_INVALID;
    *out = NULL;
    if (out_len) *out_len = 0;

    FILE *f = fopen(path, "rb");
    if (!f) return HP_ERR_IO;
    if (fseek(f, 0, SEEK_END) != 0) {
        fclose(f);
        return HP_ERR_IO;
    }
    long end = ftell(f);
    if (end < 0) {
        fclose(f);
        return HP_ERR_IO;
    }
    if (fseek(f, 0, SEEK_SET) != 0) {
        fclose(f);
        return HP_ERR_IO;
    }
    size_t len = (size_t)end;
    char *buf = (char *)hp_xmalloc(len + 1);
    size_t got = fread(buf, 1, len, f);
    fclose(f);
    if (got != len) {
        hp_free(buf);
        return HP_ERR_IO;
    }
    buf[len] = '\0';
    *out = buf;
    if (out_len) *out_len = len;
    return HP_OK;
}

hp_status hp_write_file(const char *path, const char *data, size_t len) {
    if (!path || (!data && len)) return HP_ERR_INVALID;
    FILE *f = fopen(path, "wb");
    if (!f) return HP_ERR_IO;
    size_t wrote = len ? fwrite(data, 1, len, f) : 0;
    int bad = (wrote != len);
    if (fclose(f) != 0) bad = 1;
    return bad ? HP_ERR_IO : HP_OK;
}

hp_status hp_write_text(const char *path, const char *text) {
    return hp_write_file(path, text, text ? strlen(text) : 0);
}

char *hp_strdup(const char *s) {
    if (!s) return NULL;
    size_t n = strlen(s);
    char *p = (char *)hp_xmalloc(n + 1);
    memcpy(p, s, n + 1);
    return p;
}

char *hp_trim(char *s) {
    if (!s) return NULL;
    char *p = s;
    while (*p && isspace((unsigned char)*p)) ++p;
    size_t n = strlen(p);
    while (n > 0 && isspace((unsigned char)p[n - 1])) --n;
    p[n] = '\0';
    if (p != s) memmove(s, p, n + 1);
    return s;
}

int hp_str_ieq(const char *a, const char *b) {
    if (!a || !b) return a == b;
    while (*a && *b) {
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) return 0;
        ++a;
        ++b;
    }
    return *a == '\0' && *b == '\0';
}
