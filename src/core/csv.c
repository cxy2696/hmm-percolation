#include "hp/core/csv.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hp/core/io.h"
#include "hp/core/memory.h"

static hp_status split_fields(char *line, size_t n_cols, char **out) {
    char *p = line;
    for (size_t c = 0; c < n_cols; ++c) {
        char *q = strchr(p, ',');
        if (q) *q = '\0';
        out[c] = hp_strdup(hp_trim(p));
        if (!q) {
            for (size_t k = c + 1; k < n_cols; ++k) out[k] = hp_strdup("");
            return HP_OK;
        }
        p = q + 1;
    }
    hp_trim(p);
    if (*p) return HP_ERR_FORMAT; /* too many columns */
    return HP_OK;
}

static void free_fields(char **fields, size_t n) {
    if (!fields) return;
    for (size_t i = 0; i < n; ++i) hp_free(fields[i]);
    hp_free(fields);
}

hp_status hp_csv_read(const char *path, hp_csv *out) {
    if (!path || !out) return HP_ERR_INVALID;
    memset(out, 0, sizeof(*out));

    char *buf = NULL;
    size_t len = 0;
    HP_TRY(hp_read_file(path, &buf, &len));
    if (len >= 3 && (unsigned char)buf[0] == 0xEF &&
        (unsigned char)buf[1] == 0xBB && (unsigned char)buf[2] == 0xBF) {
        memmove(buf, buf + 3, len - 2);
    }

    /* Normalise CRLF/CR and split into non-empty lines. */
    char *w = buf;
    for (char *p = buf; *p;) {
        if (*p == '\r' && p[1] == '\n') {
            ++p;
        } else if (*p == '\r') {
            *w++ = '\n';
            ++p;
        } else {
            *w++ = *p++;
        }
    }
    *w = '\0';

    char **lines = NULL;
    size_t n_lines = 0, cap = 0;
    for (char *p = buf; p && *p;) {
        char *nl = strchr(p, '\n');
        if (nl) *nl = '\0';
        hp_trim(p);
        if (*p) {
            if (n_lines == cap) {
                cap = cap ? cap * 2 : 64;
                lines = (char **)hp_xrealloc(lines, cap * sizeof(char *));
            }
            lines[n_lines++] = p;
        }
        p = nl ? nl + 1 : NULL;
    }
    if (n_lines < 1) {
        hp_free(lines);
        hp_free(buf);
        return HP_ERR_FORMAT;
    }

    size_t n_cols = 1;
    for (const char *p = lines[0]; *p; ++p)
        if (*p == ',') ++n_cols;
    size_t n_rows = n_lines - 1;

    char **header = (char **)hp_xcalloc(n_cols, sizeof(char *));
    char **cells = (char **)hp_xcalloc(n_rows ? n_rows * n_cols : 1, sizeof(char *));

    hp_status st = split_fields(lines[0], n_cols, header);
    if (st != HP_OK) {
        free_fields(header, n_cols);
        hp_free(cells);
        hp_free(lines);
        hp_free(buf);
        return st;
    }

    for (size_t r = 0; r < n_rows; ++r) {
        char **row = cells + r * n_cols;
        st = split_fields(lines[r + 1], n_cols, row);
        if (st != HP_OK) {
            free_fields(cells, n_rows * n_cols);
            free_fields(header, n_cols);
            hp_free(lines);
            hp_free(buf);
            return st;
        }
    }

    out->n_rows = n_rows;
    out->n_cols = n_cols;
    out->header = header;
    out->cells = cells;

    hp_free(lines);
    hp_free(buf);
    return HP_OK;
}

void hp_csv_free(hp_csv *csv) {
    if (!csv) return;
    if (csv->header) free_fields(csv->header, csv->n_cols);
    if (csv->cells) free_fields(csv->cells, csv->n_rows * csv->n_cols);
    memset(csv, 0, sizeof(*csv));
}

const char *hp_csv_cell(const hp_csv *csv, size_t row, size_t col) {
    if (!csv || row >= csv->n_rows || col >= csv->n_cols) return NULL;
    return csv->cells[row * csv->n_cols + col];
}

int hp_csv_col_index(const hp_csv *csv, const char *name) {
    if (!csv || !name) return -1;
    for (size_t i = 0; i < csv->n_cols; ++i)
        if (hp_str_ieq(csv->header[i], name)) return (int)i;
    return -1;
}

hp_status hp_csv_parse_double(const char *text, double *out) {
    if (!text || !out) return HP_ERR_INVALID;
    while (*text && isspace((unsigned char)*text)) ++text;
    if (*text == '\0') {
        *out = NAN;
        return HP_OK;
    }
    if (hp_str_ieq(text, "NA") || hp_str_ieq(text, "NAN") ||
        hp_str_ieq(text, "NULL")) {
        *out = NAN;
        return HP_OK;
    }
    char *end = NULL;
    double v = strtod(text, &end);
    if (end == text) return HP_ERR_FORMAT;
    while (end && *end && isspace((unsigned char)*end)) ++end;
    if (end && *end != '\0') return HP_ERR_FORMAT;
    *out = v;
    return HP_OK;
}

struct hp_csv_writer {
    FILE *f;
    size_t n_names;
    int with_label;
    hp_status err;
    char *path;
};

static void write_number(FILE *f, double v) {
    if (isfinite(v))
        fprintf(f, "%.17g", v);
    /* non-finite -> empty cell, parsed back as NAN */
}

hp_csv_writer *hp_csv_writer_open(const char *path, const char *const *names,
                                  size_t n_names, int with_label) {
    if (!path || !names || n_names == 0) return NULL;
    hp_csv_writer *w = (hp_csv_writer *)hp_xcalloc(1, sizeof(*w));
    w->f = fopen(path, "wb");
    if (!w->f) {
        hp_free(w);
        return NULL;
    }
    w->n_names = n_names;
    w->with_label = with_label != 0;
    w->path = hp_strdup(path);
    for (size_t i = 0; i < n_names; ++i) {
        if (i) fputc(',', w->f);
        fputs(names[i] ? names[i] : "", w->f);
    }
    fputc('\n', w->f);
    if (ferror(w->f)) {
        fclose(w->f);
        remove(path);
        hp_free(w->path);
        hp_free(w);
        return NULL;
    }
    return w;
}

static hp_status writer_prefix(hp_csv_writer *w, const char *label) {
    if (w->with_label) {
        if (!label) return HP_ERR_INVALID;
        fputs(label, w->f);
        fputc(',', w->f);
    }
    return HP_OK;
}

hp_status hp_csv_writer_row(hp_csv_writer *w, const double *values,
                            size_t n_values) {
    if (!w || w->err != HP_OK) return w ? w->err : HP_ERR_INVALID;
    if (w->with_label || n_values != w->n_names) {
        w->err = HP_ERR_INVALID;
        return w->err;
    }
    for (size_t i = 0; i < n_values; ++i) {
        if (i) fputc(',', w->f);
        write_number(w->f, values[i]);
    }
    fputc('\n', w->f);
    if (ferror(w->f)) w->err = HP_ERR_IO;
    return w->err;
}

hp_status hp_csv_writer_row_label(hp_csv_writer *w, const char *label,
                                  const double *values, size_t n_values) {
    if (!w || w->err != HP_OK) return w ? w->err : HP_ERR_INVALID;
    if (!w->with_label || n_values + 1 != w->n_names) {
        w->err = HP_ERR_INVALID;
        return w->err;
    }
    HP_TRY(writer_prefix(w, label));
    for (size_t i = 0; i < n_values; ++i) {
        if (i) fputc(',', w->f);
        write_number(w->f, values[i]);
    }
    fputc('\n', w->f);
    if (ferror(w->f)) w->err = HP_ERR_IO;
    return w->err;
}

hp_status hp_csv_writer_close(hp_csv_writer *w) {
    if (!w) return HP_ERR_INVALID;
    hp_status st = w->err;
    if (fclose(w->f) != 0) st = HP_ERR_IO;
    hp_free(w->path);
    hp_free(w);
    return st;
}

void hp_csv_writer_abort(hp_csv_writer *w) {
    if (!w) return;
    fclose(w->f);
    if (w->path) {
        remove(w->path);
        hp_free(w->path);
    }
    hp_free(w);
}

hp_status hp_csv_write_numeric(const char *path, const char *const *names,
                               size_t n_cols, const double *data,
                               size_t n_rows) {
    if (!path || !names || n_cols == 0 || (!data && n_rows)) return HP_ERR_INVALID;
    hp_csv_writer *w = hp_csv_writer_open(path, names, n_cols, 0);
    if (!w) return HP_ERR_IO;
    for (size_t r = 0; r < n_rows; ++r) {
        hp_status st = hp_csv_writer_row(w, data + r * n_cols, n_cols);
        if (st != HP_OK) {
            hp_csv_writer_abort(w);
            return st;
        }
    }
    return hp_csv_writer_close(w);
}
