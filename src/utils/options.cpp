#include "utils/options.hpp"

namespace pairhmm {
template <typename... Args>
auto parse(bpo::options_description &options,
           bpo::positional_options_description &positional,
           bool allow_unregistered, const Args &...args) {}
auto argparse(int argc, char **argv) {}
} // namespace pairhmm_test