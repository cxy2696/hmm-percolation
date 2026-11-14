/* HP · logging.  Quiet by default in library code; CLIs opt in. */
#ifndef HP_CORE_LOG_H
#define HP_CORE_LOG_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    HP_LOG_DEBUG = 0,
    HP_LOG_INFO = 1,
    HP_LOG_WARN = 2,
    HP_LOG_ERROR = 3,
    HP_LOG_OFF = 4
} hp_log_level;

void hp_log_set_level(hp_log_level level);
hp_log_level hp_log_get_level(void);
void hp_log_message(hp_log_level level, const char *file, int line,
                    const char *fmt, ...);

#define HP_LOGD(...) hp_log_message(HP_LOG_DEBUG, __FILE__, __LINE__, __VA_ARGS__)
#define HP_LOGI(...) hp_log_message(HP_LOG_INFO, __FILE__, __LINE__, __VA_ARGS__)
#define HP_LOGW(...) hp_log_message(HP_LOG_WARN, __FILE__, __LINE__, __VA_ARGS__)
#define HP_LOGE(...) hp_log_message(HP_LOG_ERROR, __FILE__, __LINE__, __VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif /* HP_CORE_LOG_H */
