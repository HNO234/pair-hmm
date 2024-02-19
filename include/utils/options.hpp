#pragma once
#include <boost/program_options.hpp>

namespace pairhmm {
namespace bpo = boost::program_options;
template <typename... Args>
auto parse(bpo::options_description &options,
           bpo::positional_options_description &positional,
           bool allow_unregistered, const Args &...args);
auto argparse(int argc, char **argv);
} // namespace pairhmm_test