#include "pairhmm/suzuki_pairhmm.hpp"
#include "pairhmm/fixed_point_suzuki_pairhmm.hpp"
#include "pairhmm/nw_pairhmm.hpp"
#include "utils/constant.hpp"
#include "utils/options.hpp"
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/file_io/bam.hpp>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <variant>
#include <iomanip>
#include "utils/utils.hpp"

int main(int argc, char **argv) {
  auto generic_vm = pairhmm::argparse(argc, argv);
  // Set haplotype
  auto haplotype_fa = generic_vm["haplotype-fasta"].as<std::string>();
  auto haplotype = pairhmm::get_haplotype(haplotype_fa);
  // Set read
  auto read_bam = std::filesystem::path(generic_vm["read-bam"].as<std::string>());
  auto reads = pairhmm::get_reads(read_bam);
  // Set GOP & GCP
  auto tables_file_name = generic_vm["tables"].as<std::string>();
  auto [gop, gcp] = pairhmm::get_gop_gcp(tables_file_name);
  // Set pairHMM algorithm
  std::cout << std::fixed << std::setprecision(15);
  int cnt = 0;
  for (auto &read_record : reads) {
    auto read = biovoltron::Codec::to_istring(read_record.seq);
    auto scores = pairhmm::score_to_phred(read_record.qual);
    std::cout << haplotype << std::endl;
    std::cout << read << std::endl;
    auto pairhmm_nw = pairhmm::NWPairHMM<long double>(haplotype, read, gop, gcp, scores);
    std::cout << pairhmm_nw.get_align_score() << std::endl;
    auto pairhmm_sz = pairhmm::SuzukiPairHMM<long double>(haplotype, read, gop, gcp, scores);
    std::cout << pairhmm_sz.get_align_score() << std::endl;
    auto pairhmm_fixed = pairhmm::FixedPointSuzukiPairHMM<long double>(haplotype, read, gop, gcp, scores);
    std::cout << pairhmm_fixed.get_align_score() << std::endl;
    std::cout << std::endl;
    cnt++;
    if (cnt > 100)
      break;
  }
}