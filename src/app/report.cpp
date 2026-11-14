#include "hp/app/report.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>

namespace hp::app {

void Json::escape(const std::string &s) {
    out_ += '"';
    for (unsigned char c : s) {
        switch (c) {
            case '"': out_ += "\\\""; break;
            case '\\': out_ += "\\\\"; break;
            case '\n': out_ += "\\n"; break;
            case '\r': out_ += "\\r"; break;
            case '\t': out_ += "\\t"; break;
            default:
                if (c < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    out_ += buf;
                } else {
                    out_ += static_cast<char>(c);
                }
        }
    }
    out_ += '"';
}

void Json::prefix_value() {
    if (stack_.empty()) return;
    Frame &f = stack_.back();
    if (f.is_obj) {
        if (!f.after_key) throw std::runtime_error("JSON object value without key");
    } else if (!f.first) {
        out_ += ',';
    }
    f.first = false;
    f.after_key = false;
}

Json &Json::obj() {
    prefix_value();
    out_ += '{';
    stack_.push_back({true, true, false});
    return *this;
}

Json &Json::arr() {
    prefix_value();
    out_ += '[';
    stack_.push_back({false, true, false});
    return *this;
}

Json &Json::end() {
    if (stack_.empty()) throw std::runtime_error("JSON end without open container");
    bool obj = stack_.back().is_obj;
    stack_.pop_back();
    out_ += obj ? '}' : ']';
    return *this;
}

Json &Json::key(const std::string &k) {
    if (stack_.empty() || !stack_.back().is_obj) throw std::runtime_error("JSON key outside object");
    Frame &f = stack_.back();
    if (!f.first) out_ += ',';
    f.first = false;
    escape(k);
    out_ += ':';
    f.after_key = true;
    return *this;
}

Json &Json::val(double v) {
    prefix_value();
    if (!std::isfinite(v)) {
        out_ += "null";
    } else {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.17g", v);
        out_ += buf;
    }
    return *this;
}

Json &Json::val(int v) {
    prefix_value();
    out_ += std::to_string(v);
    return *this;
}

Json &Json::val(unsigned int v) {
    prefix_value();
    out_ += std::to_string(v);
    return *this;
}

Json &Json::val(size_t v) {
    prefix_value();
    out_ += std::to_string(v);
    return *this;
}

Json &Json::val(unsigned long long v) {
    prefix_value();
    out_ += std::to_string(v);
    return *this;
}

Json &Json::val(bool v) {
    prefix_value();
    out_ += v ? "true" : "false";
    return *this;
}

Json &Json::val(const std::string &v) {
    prefix_value();
    escape(v);
    return *this;
}

Json &Json::val(const char *v) {
    prefix_value();
    escape(v ? v : "");
    return *this;
}

void write_text(const std::string &path, const std::string &text) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot write " + path);
    f.write(text.data(), static_cast<std::streamsize>(text.size()));
    if (!f) throw std::runtime_error("write failed " + path);
}

}  // namespace hp::app
