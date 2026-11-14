#include "hp/core/status.h"

const char *hp_status_name(hp_status st) {
    switch (st) {
        case HP_OK: return "HP_OK";
        case HP_ERR_INVALID: return "HP_ERR_INVALID";
        case HP_ERR_NOMEM: return "HP_ERR_NOMEM";
        case HP_ERR_IO: return "HP_ERR_IO";
        case HP_ERR_FORMAT: return "HP_ERR_FORMAT";
        case HP_ERR_DOMAIN: return "HP_ERR_DOMAIN";
        case HP_ERR_CONVERGE: return "HP_ERR_CONVERGE";
        case HP_ERR_UNSUPPORTED: return "HP_ERR_UNSUPPORTED";
        case HP_ERR_RANGE: return "HP_ERR_RANGE";
        case HP_ERR_STATE: return "HP_ERR_STATE";
        default: return "HP_ERR_UNKNOWN";
    }
}

const char *hp_status_message(hp_status st) {
    switch (st) {
        case HP_OK: return "success";
        case HP_ERR_INVALID: return "invalid argument";
        case HP_ERR_NOMEM: return "out of memory";
        case HP_ERR_IO: return "I/O error";
        case HP_ERR_FORMAT: return "malformed input";
        case HP_ERR_DOMAIN: return "numerical domain error";
        case HP_ERR_CONVERGE: return "iteration limit reached without convergence";
        case HP_ERR_UNSUPPORTED: return "unsupported operation or layout";
        case HP_ERR_RANGE: return "index or size out of range";
        case HP_ERR_STATE: return "object is in the wrong state";
        default: return "unknown error";
    }
}
