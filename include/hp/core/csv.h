/* HP · minimal CSV reader/writer.
 *
 * Scope on purpose: comma-separated, optional UTF-8 BOM, LF/CRLF, no embedded
 * newlines or quotes.  Data panels and experiment tables fit that shape.
 * Empty cells and NA/NaN/nan become NAN when parsed as numbers.
 */
#ifndef HP_CORE_CSV_H
#define HP_CORE_CSV_H

#include <stddef.h>
#include "hp/core/status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t n_rows;
    size_t n_cols;
    char **header;  /* n_cols strings */
    char **cells;   /* n_rows * n_cols strings, row-major, never NULL */
} hp_csv;

/* Reads the whole file.  On failure *out is zeroed and nothing leaks. */
hp_status hp_csv_read(const char *path, hp_csv *out);
void hp_csv_free(hp_csv *csv);

const char *hp_csv_cell(const hp_csv *csv, size_t row, size_t col);
/* -1 when absent. */
int hp_csv_col_index(const hp_csv *csv, const char *name);
/* "": NAN, "NA": NAN, "NaN": NAN; otherwise full-string strtod. */
hp_status hp_csv_parse_double(const char *text, double *out);

/* Streaming writer.  Names[0..n_names-1] form the header.  When
 * `with_label` is 1 the first name is a text label column and rows are
 * written with hp_csv_writer_row_label(). */
typedef struct hp_csv_writer hp_csv_writer;

hp_csv_writer *hp_csv_writer_open(const char *path, const char *const *names,
                                  size_t n_names, int with_label);
hp_status hp_csv_writer_row(hp_csv_writer *w, const double *values,
                            size_t n_values);
hp_status hp_csv_writer_row_label(hp_csv_writer *w, const char *label,
                                  const double *values, size_t n_values);
hp_status hp_csv_writer_close(hp_csv_writer *w);
/* Close and delete the partial file. */
void hp_csv_writer_abort(hp_csv_writer *w);

/* One-shot numeric table: header names[n_cols], data row-major. */
hp_status hp_csv_write_numeric(const char *path, const char *const *names,
                               size_t n_cols, const double *data,
                               size_t n_rows);

#ifdef __cplusplus
}
#endif

#endif /* HP_CORE_CSV_H */
