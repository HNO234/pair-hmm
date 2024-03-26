#include "pairhmm/naive_pairhmm.hpp"
#include "pairhmm/nw_pairhmm.hpp"
#include "pairhmm/suzuki_pairhmm.hpp"
#include "utils/constant.hpp"
#include "utils/options.hpp"
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/utility/istring.hpp>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <variant>

int main(int argc, char **argv) {
  auto generic_vm = pairhmm::argparse(argc, argv);
  // Set haplotype
  auto haplotype_fa =
      std::ifstream{generic_vm["haplotype-fasta"].as<std::string>()};
  auto haplotype = biovoltron::istring{};
  auto haplotype_refs =
      std::ranges::istream_view<biovoltron::FastaRecord<true>>(haplotype_fa);
  for (auto &ref : haplotype_refs)
    haplotype += ref.seq;
  std::ranges::transform(haplotype, haplotype.begin(),
                         [](auto &c) { return c % 4; });
  // Set read
  auto read_fa = std::ifstream{generic_vm["read-fatsa"].as<std::string>()};
  auto read = biovoltron::istring{};
  auto read_refs =
      std::ranges::istream_view<biovoltron::FastaRecord<true>>(read_fa);
  for (auto &ref : read_refs)
    read += ref.seq;
  std::ranges::transform(read, read.begin(), [](auto &c) { return c % 4; });
  // Set GOP & GCP
  auto gop = pairhmm::standard_gop;
  auto gcp = pairhmm::standard_gcp;
  // Set pairHMM algorithm
  std::variant<pairhmm::NaivePairHMM<double>, pairhmm::NWPairHMM<double>,
    pairhmm::SuzukiPairHMM<double>> pairhmm;
  auto algo = generic_vm["pairhmm-algorithm"].as<pairhmm::PairHMMAlgorithm>();
  if (algo == pairhmm::PairHMMAlgorithm::NAIVE) {
    pairhmm = pairhmm::NaivePairHMM<double>(haplotype, read, gop, gcp);
  } else if (algo == pairhmm::PairHMMAlgorithm::NW) {
    pairhmm = pairhmm::NWPairHMM<double>(haplotype, read, gop, gcp);
  } else if (algo == pairhmm::PairHMMAlgorithm::SUZUKI_KASAHARA) {
    pairhmm = pairhmm::SuzukiPairHMM<double>(haplotype, read, gop, gcp);
  } else {
    throw std::invalid_argument("Invalid pairHMM algorithm");
  }
  // Get cigar
  auto cigar =
      std::visit([](auto &&pairhmm) { return pairhmm.get_cigar(); }, pairhmm);
  std::cout << std::string(cigar) << std::endl;
}