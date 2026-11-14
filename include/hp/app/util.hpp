/* HP app · small RAII and filesystem helpers. */
#ifndef HP_APP_UTIL_HPP
#define HP_APP_UTIL_HPP

#include <memory>
#include <string>

#include "hp/hmm/ghmm.h"
#include "hp/market/panel.h"

namespace hp::app {

using PanelPtr = std::unique_ptr<hp_panel, void (*)(hp_panel *)>;
using GhmmPtr = std::unique_ptr<hp_ghmm, void (*)(hp_ghmm *)>;

PanelPtr load_panel(const std::string &path);
GhmmPtr make_ghmm(size_t n_states, size_t dim, size_t n_mix);

std::string path_join(const std::string &dir, const std::string &name);
void ensure_dir(const std::string &path);
void ensure_parent(const std::string &path);
std::string fmt(double v, int precision = 6);
std::string now_utc();

}  // namespace hp::app

#endif
