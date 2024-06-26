#include "pairhmm/naive_pairhmm.hpp"
#include "pairhmm/nw_pairhmm.hpp"
#include "pairhmm/suzuki_pairhmm.hpp"
#include "utils/constant.hpp"
#include <biovoltron/utility/istring.hpp>
#include <catch_amalgamated.hpp>
#include <cmath>
#include <experimental/random>
#include <string>
#include <string_view>
#include <vector>
#include <iostream>

bool is_close(double a, double b) { return std::abs(a - b) < 1e-8; }

TEST_CASE("test repeat") {
  std::string s = "ATTTTTTCAATGTTTACACATTTCCTTCCTCCCTCCCTCCTTCCTTTCCTCCCTTCCTCC"
                  "CTTCCTCCCTTCCTTCCTGTTTGCTTTATTATTGTATTG";
  std::vector<size_t> best_period = {1, 
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
      1, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 8, 8, 1,
      1, 1, 1, 8, 8, 8, 8, 1, 1, 1, 1, 8, 8, 8, 8, 1, 1, 1, 1, 8, 8, 8, 8, 8, 1,
      1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1};
  std::vector<size_t> best_repeat = {1, 
      6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 1, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2,
      2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2,
      2, 2, 2, 3, 3, 3, 3, 1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1};
  auto haplotype = biovoltron::Codec::to_istring(std::string_view(s));
  auto read = biovoltron::Codec::to_istring(std::string_view(s));
  auto gop = pairhmm::standard_gop;
  auto gcp = pairhmm::standard_gcp;
  auto pairhmm_instance = pairhmm::NaivePairHMM<double>(haplotype, read, gop, gcp);
  // auto [calculated_period_haplotype, calculated_repeat_haplotype] =
  //     pairhmm_instance.get_haplotype_best_repeat();
  // REQUIRE(calculated_period_haplotype == best_period);
  // REQUIRE(calculated_repeat_haplotype == best_repeat);
  auto [calculated_period_read, calculated_repeat_read] =
      pairhmm_instance.get_read_best_repeat();
  REQUIRE(calculated_period_read.size() == best_period.size());
  REQUIRE(calculated_period_read == best_period);
  REQUIRE(calculated_repeat_read == best_repeat);
}

TEST_CASE("test cigar") {
  biovoltron::istring haplotype(20, 0), read(4, 0);
  for (size_t i = size_t{}; i < 20; i++)
    haplotype[i] = std::experimental::randint(0, 3);
  for (size_t i = size_t{}; i < 8; i++)
    read[i] = rand() % 4;
  auto gop = pairhmm::standard_gop;
  auto gcp = pairhmm::standard_gcp;
  auto naive_pairhmm = pairhmm::NaivePairHMM<double>(haplotype, read, gop, gcp);
  // auto nw_pairhmm = pairhmm::NWPairHMM<double>(haplotype, read, gop, gcp);
  std::vector<int> score;
  auto suzuki_pairhmm = pairhmm::SuzukiPairHMM<double>(haplotype, read, gop, gcp, score);
  naive_pairhmm.run_alignment();
  // nw_pairhmm.run_alignment();
  suzuki_pairhmm.run_alignment();
  // REQUIRE(naive_pairhmm.get_cigar() == nw_pairhmm.get_cigar());
  REQUIRE(naive_pairhmm.get_cigar() == suzuki_pairhmm.get_cigar());
}

TEST_CASE("test tables") {
  // data preparation
  biovoltron::istring haplotype(64, 0), read(16, 0);
  for (size_t i = size_t{}; i < 64; i++)
    haplotype[i] = std::experimental::randint(0, 3);
  for (size_t i = size_t{}; i < 16; i++)
    read[i] = std::experimental::randint(0, 3);
  auto gop = pairhmm::standard_gop;
  auto gcp = pairhmm::standard_gcp;
  // run alignment
  auto naive_pairhmm = pairhmm::NaivePairHMM<double>(haplotype, read, gop, gcp);
  // auto nw_pairhmm = pairhmm::NWPairHMM<double>(haplotype, read, gop, gcp);
  std::vector<int> score;
  auto suzuki_pairhmm = pairhmm::SuzukiPairHMM<double>(haplotype, read, gop, gcp, score);
  naive_pairhmm.run_alignment();
  // nw_pairhmm.run_alignment();
  suzuki_pairhmm.run_alignment();
  // get tables
  auto M = naive_pairhmm.get_M();
  auto I = naive_pairhmm.get_I();
  auto D = naive_pairhmm.get_D();
  auto S_calculated = M * (-10);
  auto E_calculated = D * (-10);
  auto F_calculated = I * (-10);
  // auto S = nw_pairhmm.get_S();
  // auto S_down = S.shifted_down();
  // auto S_right = S.shifted_right();
  // auto E = nw_pairhmm.get_E();
  // auto F = nw_pairhmm.get_F();
  // auto DeltaV_calculated = S - S_right;
  // auto DeltaH_calculated = S - S_down;
  // auto DeltaE_calculated = E - S;
  // auto DeltaF_calculated = F - S;
  // auto DeltaEp_calculated = DeltaE_calculated + DeltaV_calculated;
  // auto DeltaFp_calculated = DeltaF_calculated + DeltaH_calculated;

  auto DeltaV = suzuki_pairhmm.get_DeltaV();
  auto DeltaH = suzuki_pairhmm.get_DeltaH();
  auto DeltaEp = suzuki_pairhmm.get_DeltaEp();
  auto DeltaFp = suzuki_pairhmm.get_DeltaFp();
  // check table consistency
  // REQUIRE(table::is_close(S_calculated, S, 1e-8));
  // REQUIRE(table::is_close(E_calculated, E, 1e-8));
  // REQUIRE(table::is_close(F_calculated, F, 1e-8));
  // REQUIRE(table::is_close(DeltaV_calculated, DeltaV, 1e-8));
  // REQUIRE(table::is_close(DeltaH_calculated, DeltaH, 1e-8));
  // REQUIRE(table::is_close(DeltaEp_calculated, DeltaEp, 1e-8));
  // REQUIRE(table::is_close(DeltaFp_calculated, DeltaFp, 1e-8));
}