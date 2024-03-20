#include "pairhmm/naive_pairhmm.hpp"
#include "pairhmm/suzuki_pairhmm.hpp"
#include "utils/constant.hpp"
#include <biovoltron/utility/istring.hpp>
#include <catch_amalgamated.hpp>
#include <cmath>
#include <random>
#include <string>
#include <string_view>
#include <vector>

bool is_close(double a, double b) { return std::abs(a - b) < 1e-8; }

TEST_CASE("test repeat") {
  std::string s = "ATTTTTTCAATGTTTACACATTTCCTTCCTCCCTCCCTCCTTCCTTTCCTCCCTTCCTCC"
                  "CTTCCTCCCTTCCTTCCTGTTTGCTTTATTATTGTATTG";
  std::vector<size_t> best_period = {
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
      1, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 8, 8, 1,
      1, 1, 1, 8, 8, 8, 8, 1, 1, 1, 1, 8, 8, 8, 8, 1, 1, 1, 1, 8, 8, 8, 8, 8, 1,
      1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 1, 1, 1};
  std::vector<size_t> best_repeat = {
      6, 6, 6, 6, 6, 6, 6, 2, 2, 2, 1, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2,
      2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2,
      2, 2, 2, 3, 3, 3, 3, 1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  auto haplotype = biovoltron::Codec::to_istring(std::string_view(s));
  auto read = biovoltron::Codec::to_istring(std::string_view(s));
  auto gop = pairhmm::standard_gop;
  auto gcp = pairhmm::standard_gcp;
  auto pairhmm_instance = pairhmm::NaivePairHMM<double>(haplotype, read, gop, gcp);
  auto [calculated_period_haplotype, calculated_repeat_haplotype] =
      pairhmm_instance.get_haplotype_best_repeat();
  REQUIRE(calculated_period_haplotype == best_period);
  REQUIRE(calculated_repeat_haplotype == best_repeat);
  auto [calculated_period_read, calculated_repeat_read] =
      pairhmm_instance.get_read_best_repeat();
  REQUIRE(calculated_period_read == best_period);
  REQUIRE(calculated_repeat_read == best_repeat);
}

// TEST_CASE("test cigar") {}

TEST_CASE("test tables") {
  biovoltron::istring haplotype(20, 0), read(4, 0);
  for (size_t i = size_t{}; i < 20; i++)
    haplotype[i] = rand() % 4;
  for (size_t i = size_t{}; i < 8; i++)
    read[i] = rand() % 4;
  auto gop = pairhmm::standard_gop;
  auto gcp = pairhmm::standard_gcp;
  auto naive_pairhmm = pairhmm::NaivePairHMM<double>(haplotype, read, gop, gcp);
  auto suzuki_pairhmm = pairhmm::SuzukiPairHMM<double>(haplotype, read, gop, gcp);
  auto S = naive_pairhmm.get_S();
  auto S_down = S.shifted_down();
  auto S_right = S.shifted_right();
  auto E = naive_pairhmm.get_E();
  auto F = naive_pairhmm.get_F();
  auto DeltaV_calculated = S - S_down;
  auto DeltaH_calculated = S - S_right;
  auto DeltaE_calculated = E - S;
  auto DeltaF_calculated = F - S;
  auto DeltaEp_calculated = DeltaE_calculated + DeltaV_calculated;
  auto DeltaFp_calculated = DeltaF_calculated + DeltaH_calculated;

  auto DeltaV = suzuki_pairhmm.get_DeltaV();
  auto DeltaH = suzuki_pairhmm.get_DeltaH();
  auto DeltaEp = suzuki_pairhmm.get_DeltaEp();
  auto DeltaFp = suzuki_pairhmm.get_DeltaFp();

  REQUIRE(table::is_close(DeltaV_calculated, DeltaV, 1e-8));
  REQUIRE(table::is_close(DeltaH_calculated, DeltaH, 1e-8));
  REQUIRE(table::is_close(DeltaEp_calculated, DeltaEp, 1e-8));
  REQUIRE(table::is_close(DeltaFp_calculated, DeltaFp, 1e-8));
}