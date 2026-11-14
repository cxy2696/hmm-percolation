#include "hp/app/args.hpp"

#include <cstdlib>
#include <stdexcept>

namespace hp::app {

Args::Args(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a.rfind("--", 0) == 0) {
            std::string key = a.substr(2);
            std::string value;
            auto eq = key.find('=');
            if (eq != std::string::npos) {
                value = key.substr(eq + 1);
                key = key.substr(0, eq);
            } else if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) {
                value = argv[++i];
            } else {
                flags_.insert(key);
                value = "1";
            }
            values_[key] = value;
        } else {
            positional_.push_back(a);
        }
    }
}

bool Args::has(const std::string &key) const { return values_.count(key) != 0; }

std::string Args::get(const std::string &key, const std::string &fallback) const {
    auto it = values_.find(key);
    return it == values_.end() ? fallback : it->second;
}

std::string Args::require(const std::string &key) const {
    auto it = values_.find(key);
    if (it == values_.end()) throw std::runtime_error("missing required argument --" + key);
    return it->second;
}

double Args::get_double(const std::string &key, double fallback) const {
    auto it = values_.find(key);
    if (it == values_.end()) return fallback;
    char *end = nullptr;
    double v = std::strtod(it->second.c_str(), &end);
    if (!end || *end != '\0') throw std::runtime_error("invalid number for --" + key + ": " + it->second);
    return v;
}

long Args::get_int(const std::string &key, long fallback) const {
    auto it = values_.find(key);
    if (it == values_.end()) return fallback;
    char *end = nullptr;
    long v = std::strtol(it->second.c_str(), &end, 10);
    if (!end || *end != '\0') throw std::runtime_error("invalid integer for --" + key + ": " + it->second);
    return v;
}

}  // namespace hp::app
