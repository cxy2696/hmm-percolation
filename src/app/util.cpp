#include "hp/app/util.hpp"

#include <chrono>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <stdexcept>

#include "hp/core/status.h"

namespace hp::app {

PanelPtr load_panel(const std::string &path) {
    hp_status st = HP_OK;
    hp_panel *p = hp_panel_load_csv(path.c_str(), &st);
    if (!p) throw std::runtime_error(std::string("cannot load panel ") + path + ": " + hp_status_message(st));
    return PanelPtr(p, hp_panel_free);
}

GhmmPtr make_ghmm(size_t n_states, size_t dim, size_t n_mix) {
    hp_ghmm *m = hp_ghmm_new(n_states, dim, n_mix);
    if (!m) throw std::runtime_error("hp_ghmm_new failed");
    return GhmmPtr(m, hp_ghmm_free);
}

std::string path_join(const std::string &dir, const std::string &name) {
    if (dir.empty()) return name;
    if (dir.back() == '/') return dir + name;
    return dir + "/" + name;
}

void ensure_dir(const std::string &path) {
    if (path.empty()) return;
    std::error_code ec;
    std::filesystem::create_directories(path, ec);
    if (ec) throw std::runtime_error("cannot create directory " + path + ": " + ec.message());
}

void ensure_parent(const std::string &path) {
    std::filesystem::path p(path);
    auto parent = p.parent_path();
    if (!parent.empty()) ensure_dir(parent.string());
}

std::string fmt(double v, int precision) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", precision, v);
    return std::string(buf);
}

std::string now_utc() {
    using namespace std::chrono;
    auto t = system_clock::to_time_t(system_clock::now());
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return std::string(buf);
}

}  // namespace hp::app
