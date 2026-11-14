/* HP · library version.
 *
 * Keep this in sync with CMakeLists.txt and README.md.  Bindings read the
 * runtime string instead of duplicating the number.
 */
#ifndef HP_CORE_VERSION_H
#define HP_CORE_VERSION_H

#ifdef __cplusplus
extern "C" {
#endif

#define HP_VERSION_MAJOR 0
#define HP_VERSION_MINOR 1
#define HP_VERSION_PATCH 4
#define HP_VERSION_STRING "0.1.4"

const char *hp_version(void);

#ifdef __cplusplus
}
#endif

#endif /* HP_CORE_VERSION_H */
