/* HP app · tiny command-line parser (--key value, --key=value, --flag). */
#ifndef HP_APP_ARGS_HPP
#define HP_APP_ARGS_HPP

#include <map>
#include <set>
#include <string>
#include <vector>

namespace hp::app {

class Args {
public:
    Args(int argc, char **argv);

    bool has(const std::string &key) const;
    std::string get(const std::string &key, const std::string &fallback = "") const;
    std::string require(const std::string &key) const;
    double get_double(const std::string &key, double fallback) const;
    long get_int(const std::string &key, long fallback) const;
    const std::vector<std::string> &positional() const { return positional_; }

private:
    std::map<std::string, std::string> values_;
    std::set<std::string> flags_;
    std::vector<std::string> positional_;
};

}  // namespace hp::app

#endif
