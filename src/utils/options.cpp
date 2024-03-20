#include "utils/options.hpp"
#include "utils/constant.hpp"
#include <boost/algorithm/string.hpp>
#include <iostream>

namespace pairhmm {
std::istream &operator>>(std::istream &in, PairHMMAlgorithm &algorithm) {
  std::string token;
  in >> token;
  boost::to_upper(token);
  if (token == "NAIVE") {
    algorithm = PairHMMAlgorithm::NAIVE;
  } else if (token == "SUZUKI_KASAHARA") {
    algorithm = PairHMMAlgorithm::SUZUKI_KASAHARA;
  } else {
    throw boost::program_options::validation_error(
        boost::program_options::validation_error::invalid_option_value);
  }
  return in;
}

auto options = [] { // {{{
  auto options = bpo::options_description{"Options"};
  options.add_options()("help,h", "produce help message")(
      "pairhmm-algorithm,p",
      bpo::value<pairhmm::PairHMMAlgorithm>()
          ->default_value(pairhmm::PairHMMAlgorithm::NAIVE, "NAIVE")
          ->value_name("ALGO"),
      "The PairHMM algorithm.\n"
      "Valid arguments are NAIVE and SUZUKI_KASAHARA.");
  return options;
}(); // }}}

auto [options_cmdline, options_positional] = [] { // {{{
  auto options_hidden = bpo::options_description{};
  options_hidden.add_options()(
      "haplotype-fasta", bpo::value<std::string>()->value_name("HAPLOTYPE"),
      "Haplotype fasta")("read-fasta",
                         bpo::value<std::string>()->value_name("READ"),
                         "Read fasta");

  auto options_positional = bpo::positional_options_description{};
  options_positional.add("haplotype-fasta", 1).add("read-fasta", 1);

  auto options_cmdline = bpo::options_description{};
  options_cmdline.add(options).add(options_hidden);

  return std::tuple{options_cmdline, options_positional};
}(); // }}}

auto help_options = [] {
  namespace bpo = boost::program_options;
  auto help_options = bpo::options_description(
      "pairhmm [--options...] haplotype-fasta read-fasta");
  help_options.add(options);
  return help_options;
}();

template <typename... Args>
auto parse(bpo::options_description &options,
           bpo::positional_options_description &positional,
           const Args &...args) {
  try {
    auto parsed = bpo::command_line_parser(args...).options(options).positional(
        positional);
    auto parsed_run = parsed.run();
    auto vm = bpo::variables_map{};
    bpo::store(parsed_run, vm);
    bpo::notify(vm);
    return std::tuple{parsed_run, vm};
  } catch (bpo::error &e) {
    std::cerr << help_options << std::endl;
    std::cerr << e.what() << std::endl;
    exit(1);
  }
}
bpo::variables_map argparse(int argc, char **argv) {
  namespace bpo = boost::program_options;

  auto [parsed, vm] = parse(options_cmdline, options_positional, argc, argv);

  if (vm.count("help")) {
    std::cerr << help_options << std::endl;
    exit(1);
  }
  return vm;
}
} // namespace pairhmm