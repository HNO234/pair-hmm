#pragma once

#include "utils/fixed_point.hpp"
#include "pairhmm/pairhmm.hpp"
#include "table/fixed_point_table.hpp"
#include "table/probability_table.hpp"
#include "utils/constant.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>


namespace pairhmm {
template <typename T, size_t P0 = 16, size_t P1 = 16, size_t P2 = 16, size_t P3 = 16, 
          size_t P4 = 16, size_t P5 = 16, size_t P6 = 16, size_t P7 = 16>
          class FixedPointSuzukiPairHMM : public PairHMM<T> {
private:
  table::FixedPointTable<FixedPoint<P5>> A;
  table::FixedPointTable<FixedPoint<P6>> DeltaH, DeltaV;
  table::FixedPointTable<FixedPoint<P7>> DeltaEp, DeltaFp;
  std::vector<size_t> period_read, repeat_read;
  std::vector<int> score;

  T get_gop_phred (size_t j) {
    if (j == 0)
      return 0;
    else if (j > this->read.size()) {
      return std::numeric_limits<double>::max();
    }
    else if (j == this->read.size()) {
      return 45;
    } else {
      auto best_period_j = period_read[j];
      auto best_repeat_j = repeat_read[j];
      best_period_j = std::min(best_period_j, this->gop.get_rows());
      best_repeat_j = std::min(best_repeat_j, this->gcp.get_columns());
      auto delta_j = this->gop.get_cell(best_period_j -  1, best_repeat_j - 1);
      return std::min((T)40, std::round(delta_j));
    }
  }

  T get_gcp_phred (size_t j) {
    if (j == 0)
      return 0;
    else if (j > this->read.size()) {
      return std::numeric_limits<double>::max();
    }
    else if (j == this->read.size()) {
      return 10;
    } else {
      auto best_period_j = period_read[j];
      auto best_repeat_j = repeat_read[j];
      best_period_j = std::min(best_period_j, this->gop.get_rows());
      best_repeat_j = std::min(best_repeat_j, this->gcp.get_columns());
      auto epislon_j = this->gcp.get_cell(best_period_j -  1, best_repeat_j - 1);
      return std::round(epislon_j);
    }
  }

  FixedPoint<P0> s(size_t i, size_t j) {
    auto match_prob = match(i, j);
    auto phred_1_minus_2_times_delta_j_minus_1 = one_minus_2_delta_j(j);
    return FixedPoint<P0>(phred_1_minus_2_times_delta_j_minus_1 + match_prob);
  }

  FixedPoint<P1> oH(size_t j) {
    auto delta_j = get_gop_phred(j);
    return FixedPoint<P1>(delta_j);
  }

  FixedPoint<P1> oV(size_t j) {
    auto delta_j = get_gop_phred(j + 1);
    return FixedPoint<P1>(delta_j);
  }

  FixedPoint<P1> one_minus_2_delta_j(size_t j) {
    if (j <= 0)
      return FixedPoint<P1>(0.0);
    auto delta_j_real = std::pow(10, -0.1 * get_gop_phred(j));
    auto phred_1_minus_2_times_delta_j = -10 * std::log10(1 - 2 * delta_j_real);
    return FixedPoint<P1>(phred_1_minus_2_times_delta_j);
  }

  FixedPoint<P2> misalign() {
    return FixedPoint<P2>(0.0);
  }

  FixedPoint<P2> match(size_t i, size_t j) {
    T base_quality = score[j - 1];
    T match_prob;
    if (this->haplotype[i - 1] == this->read[j - 1]) {
      match_prob = -10 * std::log10(1 - std::pow(10, -0.1 * base_quality));
    } else {
      match_prob = -10 * std::log10(std::pow(10, -0.1 * base_quality) / 3);
    }
    return FixedPoint<P2>(match_prob);
  }

  FixedPoint<P3> eH(size_t j) {
    return FixedPoint<P3>(epislon(j));
  }

  FixedPoint<P3> eV(size_t j) {
    return FixedPoint<P3>(epislon(j + 1));
  }

  FixedPoint<P3> cH(size_t j) {
    return FixedPoint<P3>(one_minus_epsilon_j(j + 1));
  }

  FixedPoint<P3> cV(size_t j) {
    return FixedPoint<P3>(one_minus_epsilon_j(j + 1));
  }

  FixedPoint<P4> epislon(size_t j) {
    auto epsilon_j = get_gcp_phred(j);
    return FixedPoint<P4>(epsilon_j);
  }

  FixedPoint<P4> one_minus_epsilon_j(size_t j) {
    auto epsilon_j = get_gcp_phred(j);
    auto phred_1_minus_epislon_j = -10 * std::log10(1 - std::pow(10, -0.1 * epsilon_j));
    return FixedPoint<P4>(phred_1_minus_epislon_j);
  }

  T initial_value() {
  size_t j = 1;
  auto haplotype_size = this->haplotype.size();
  auto initial_value = -10 * std::log10(1.0 / haplotype_size);
  return initial_value;
}

public:
  FixedPointSuzukiPairHMM() {}
  FixedPointSuzukiPairHMM(
    biovoltron::istring haplotype_,
    biovoltron::istring read_,
    table::STRTable<T> gop_,
    table::STRTable<T> gcp_,
    std::vector<int> score_):
    PairHMM<T>(haplotype_, read_, gop_, gcp_) {
      score = score_;
      auto haplotype_size = this->haplotype.size();
      auto read_size = this->read.size();
      // std::cerr << haplotype_size << ' ' << read_size << std::endl;
      A = table::FixedPointTable<FixedPoint<P5>>(haplotype_size + 1, read_size + 1);
      DeltaH = table::FixedPointTable<FixedPoint<P6>>(haplotype_size + 1, read_size + 1);
      DeltaV = table::FixedPointTable<FixedPoint<P6>>(haplotype_size + 1, read_size + 1);
      DeltaEp = table::FixedPointTable<FixedPoint<P7>>(haplotype_size + 1, read_size + 1);
      DeltaFp = table::FixedPointTable<FixedPoint<P7>>(haplotype_size + 1, read_size + 1);
      auto [period_read_, repeat_read_] = this->get_read_best_repeat();
      period_read = period_read_;
      repeat_read = repeat_read_;
    }
  biovoltron::Cigar get_cigar() {
    // throw std::logic_error("Function not yet implemented");
    return biovoltron::Cigar{};
  }
  void run_alignment() {
    auto haplotype_size = this->haplotype.size();
    auto read_size = this->read.size();

    // calculate values for j = 0
    for (auto i = size_t{1}; i <= haplotype_size; i++) {
      DeltaH.set_cell(i, 0, 0);
      DeltaFp.set_cell(i, 0, 80);
    }

    // calculate values for i = 1
    auto last_gap_penalty_i_0 = FixedPoint<P6>(0.0);
    auto accumulated_gap_penalty_i_0 = FixedPoint<P6>(0.0);
    for (auto j = size_t{1}; j <= read_size; j++) {
      size_t i = 1;
      // Calculate DeltaV
      if (j == 1)
        accumulated_gap_penalty_i_0 = s(i, j);
      else if (j == 2)
        accumulated_gap_penalty_i_0 = accumulated_gap_penalty_i_0 + oV(j - 1) + cV(j);
      else
        accumulated_gap_penalty_i_0 = accumulated_gap_penalty_i_0 + (-cV(j - 1) + eV(j - 1) + cV(j));
      DeltaV.set_cell(i, j, accumulated_gap_penalty_i_0 - last_gap_penalty_i_0);
      last_gap_penalty_i_0 = accumulated_gap_penalty_i_0;
      // Calculate DeltaEp
      DeltaEp.set_cell(i, j, FixedPoint<P7>(DeltaV.get_cell(i, j) + oH(j)));
    }

    // calculate values for other cells
    for (auto i = size_t{2}; i <= haplotype_size; i++)
      for (auto j = size_t{1}; j <= read_size; j++) {
        A.set_cell(i, j, 
          FixedPoint<P5>(min(s(i, j), 
          min(DeltaEp.get_cell(i - 1, j) + cH(j),
              DeltaFp.get_cell(i, j - 1) + cV(j)))));
        DeltaH.set_cell(i, j, FixedPoint<P6>(A.get_cell(i, j) - DeltaV.get_cell(i - 1, j)));
        DeltaV.set_cell(i, j, FixedPoint<P6>(A.get_cell(i, j) - DeltaH.get_cell(i, j - 1)));
        DeltaEp.set_cell(i, j, FixedPoint<P7>(min(A.get_cell(i, j) + oH(j),
          DeltaEp.get_cell(i - 1, j) + eH(j)) - DeltaH.get_cell(i, j - 1)));
        DeltaFp.set_cell(i, j, FixedPoint<P7>(min(A.get_cell(i, j) + oV(j),
          DeltaFp.get_cell(i, j - 1) + eV(j)) - DeltaV.get_cell(i - 1, j)));
    }
  }

  T get_align_score() {
    auto haplotype_size = this->haplotype.size();
    auto read_size = this->read.size();
    run_alignment();
    auto S = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
    for (auto i = size_t{0}; i <= haplotype_size; i++) {
      S.set_cell(i, 0, initial_value());
    }
    for (auto i = size_t{1}; i <= haplotype_size; i++)
      for (auto j = size_t{1}; j <= read_size; j++)
        S.set_cell(i, j, S.get_cell(i, j - 1) + DeltaV.get_cell(i, j).to_float());
    auto infty = std::numeric_limits<double>::max();
    T ans = infty;
    for (int i = 1; i <= haplotype_size; i++) {
      ans = std::min(ans, S.get_cell(i, read_size));
    }
    return ans;
  }

  auto get_resource_adder() {
    return std::max(P1, P2)
      + 4 * std::max(P2, P4)
      + 2 * std::max(P3, P7)
      + 2 * std::max(P5, P6)
      + 2 * std::max(P2, P5)
      + 2 * std::max(P3, P7)
      + 2 * std::max({P2, P3, P5, P6, P7});
  }

  auto get_resource_comp() {
    return std::max({P0, P3, P7})
      + 2 * std::max({P2, P3, P5, P7});
  }
};
} // namespace pairhmm