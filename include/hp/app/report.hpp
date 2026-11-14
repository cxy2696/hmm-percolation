/* HP app · minimal streaming JSON writer (no external dependency). */
#ifndef HP_APP_REPORT_HPP
#define HP_APP_REPORT_HPP

#include <string>
#include <vector>

namespace hp::app {

class Json {
public:
    Json &obj();
    Json &arr();
    Json &end();
    Json &key(const std::string &k);

    Json &val(double v);
    Json &val(int v);
    Json &val(size_t v);
    Json &val(unsigned int v);
    Json &val(unsigned long long v);
    Json &val(bool v);
    Json &val(const std::string &v);
    Json &val(const char *v);

    template <typename T>
    Json &kv(const std::string &k, const T &v) {
        key(k);
        return val(v);
    }

    std::string dump() const { return out_; }

private:
    struct Frame {
        bool is_obj;
        bool first;
        bool after_key;
    };
    void prefix_value();
    void escape(const std::string &s);
    std::vector<Frame> stack_;
    std::string out_;
};

void write_text(const std::string &path, const std::string &text);

}  // namespace hp::app

#endif
