#pragma once
#include <boost/program_options.hpp>

namespace pairhmm {
namespace bpo = boost::program_options;
template <typename... Args>
auto parse(bpo::options_description &options,
           bpo::positional_options_description &positional,
           const Args &...args);
bpo::variables_map argparse(int argc, char **argv);
} // namespace pairhmm