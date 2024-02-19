#include "pairhmm/naive_pairhmm.hpp"
#include "pairhmm/suzuki_pairhmm.hpp"
#include "utils/constant.hpp"
#include "utils/options.hpp"
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/utility/istring.hpp>
#include <boost/program_options.hpp>
#include <exception>
#include <iostream>
#include <ranges>
#include <string>
#include <variant>

int main(int argc, char **argv) {
  auto generic_vm = pairhmm_test::argparse(argc, argv);
  // Set haplotype
  auto haplotype_fa =
      std::ifstream{generic_vm["haplotype-fasta"].as<std::string>()};
  auto haplotype = biovoltron::istring{};
  auto haplotype_refs =
      std::ranges::istream_view<biovoltron::FastaRecord<true>>(haplotype);
  for (auto &ref : haplotype_refs)
    haplotype += ref.seq;
  std::ranges::transform(haplotype, haplotype.begin(),
                         [](auto &c) { return c % 4; });
  // Set read
  auto read_fa = std::ifstream{generic_vm["read-fasta"].as<std::string>()};
  auto read = biovoltron::istring{};
  auto read_refs =
      std::ranges::istream_view<biovoltron::FastaRecord<true>>(read);
  for (auto &ref : read_refs)
    read += ref.seq;
  std::ranges::transform(read, read.begin(), [](auto &c) { return c % 4; });
  // Set GOP & GCP
  auto gop = pairhmm::standard_gop;
  auto gcp = pairhmm::standard_gcp;
  // Set pairHMM algorithm
  std::variant<pairhmm::NaivePairHMM<double>, pairhmm::SuzukiPairHMM<double>>
      pairhmm;
  switch (generic_vm["pairhmm-algorithm"].as<pairhmm::SortingAlgorithm>()) {
  case pairhmm::SortingAlgorithm::NAIVE:
    pairhmm = pairhmm::NaivePairHMM<double>(haplotype, read, gop, gcp);
    break;
  case pairhmm::SortingAlgorithm::SUZUKI_KASAHARA:
    pairhmm = pairhmm::SuzukiPairHMM<double>();
    break;
  default:
    throw std::invalid_argument("Invalid pairHMM algorithm");
    break;
  }
  // Get cigar
  auto cigar =
      std::visit([](auto &&pairhmm) { return pairhmm.get_cigar(); }, pairhmm);
  std::cout << std::string(cigar) << std::endl;
}