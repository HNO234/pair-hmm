#pragma once

#include "pairhmm/pairhmm.hpp"
#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include "pairhmm/suzuki_pairhmm.hpp"
#include "table/probability_table.hpp"
#include "utils/constant.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>

namespace pairhmm {
template <typename T> class SuzukiPairHMM : public PairHMM<T> {
private:
  table::ProbabilityTable<T> A, DeltaH, DeltaV, DeltaEp, DeltaFp;
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

  T match(size_t i, size_t j) {
    T base_quality = score[j - 1];
    if (this->haplotype[i - 1] == this->read[j - 1]) {
      return -10 * std::log10(1 - std::pow(10, -0.1 * base_quality));
    } else {
      return -10 * std::log10(std::pow(10, -0.1 * base_quality) / 3);
    }
  }

  T misalign() {
    return 0;
  }

  T s(size_t i, size_t j) {
  auto match_prob = match(i, j);
  auto gop_phred = get_gop_phred(j);
  auto delta_j_real = (j ? std::pow(10, -0.1 * gop_phred) : 0);
  auto phred_1_minus_2_times_delta_j = -10 * std::log10(1 - 2 * delta_j_real);
  return phred_1_minus_2_times_delta_j + match_prob;
}
  T oH(size_t j) {
    return get_gop_phred(j);
}
  T oV(size_t j) {
    return get_gop_phred(j + 1);
}
  T eH(size_t j) {
  return get_gcp_phred(j);
}
  T eV(size_t j) {
  return get_gcp_phred(j + 1);
}
  T cH(size_t j) {
  auto epsilon_j = get_gcp_phred(j + 1);
  auto phred_1_minus_epislon_j = -10 * std::log10(1 - std::pow(10, -0.1 * epsilon_j));
  return phred_1_minus_epislon_j;
}
  T cV(size_t j) {
  auto epislon_j_minus_1 = get_gcp_phred(j + 1);
  auto phred_1_minus_epislon_j_minus_1 = -10 * std::log10(1 - std::pow(10, -0.1 * epislon_j_minus_1));
  return phred_1_minus_epislon_j_minus_1;
}

  T initial_value() {
  size_t j = 1;
  auto haplotype_size = this->haplotype.size();
  auto initial_value = -10 * std::log10(1.0 / haplotype_size);
  return initial_value;
}

public:
  // using PairHMM<T>::PairHMM;
  SuzukiPairHMM() {}
  SuzukiPairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
                    table::STRTable<T> gop_, table::STRTable<T> gcp_,
                    std::vector<int> score_)
    : PairHMM<T>(haplotype_, read_, gop_, gcp_) {
  score = score_;
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  A = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  DeltaH = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  DeltaV = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  DeltaEp = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  DeltaFp = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  auto [period_read_, repeat_read_] = this->get_read_best_repeat();
  period_read = period_read_;
  repeat_read = repeat_read_;
}
  biovoltron::Cigar get_cigar() {
  throw std::logic_error("Function not yet implemented");
  return biovoltron::Cigar{};
}

  void run_alignment() {
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  // calculate values for j = 0
  for (auto i = size_t{1}; i <= haplotype_size; i++) {
    DeltaH.set_cell(i, 0, 0);
    DeltaFp.set_cell(i, 0, 80); // max of match score + open gap score
  }
  // calculate values for i = 1
  auto last_gap_penalty_i_0 = T{};
  auto accumulated_gap_penalty_i_0 = T{};
  for (auto j = size_t{1}; j <= read_size; j++) {
    size_t i = 1;
    // Calculate DeltaV
    if (j == 1)
      accumulated_gap_penalty_i_0 = s(i, j);
    else if (j == 2)
      accumulated_gap_penalty_i_0 += oV(j - 1) + cV(j);
    else
      accumulated_gap_penalty_i_0 += (-cV(j - 1) + eV(j - 1) + cV(j));
    DeltaV.set_cell(i, j, accumulated_gap_penalty_i_0 - last_gap_penalty_i_0);
    last_gap_penalty_i_0 = accumulated_gap_penalty_i_0;
    // Calculate DeltaEp
    DeltaEp.set_cell(i, j, DeltaV.get_cell(i, j) + oH(j));

 auto gop_phred = get_gop_phred(j);
    auto delta_j_real = (j ? std::pow(10, -0.1 * gop_phred) : 0);
  auto phred_1_minus_2_times_delta_j = -10 * std::log10(1 - 2 * delta_j_real);
  }
  for (auto i = size_t{2}; i <= haplotype_size; i++)
    for (auto j = size_t{1}; j <= read_size; j++) {
      A.set_cell(i, j, std::min({s(i, j),
        DeltaEp.get_cell(i - 1, j) + cH(j),
        DeltaFp.get_cell(i, j - 1) + cV(j)}));
      DeltaH.set_cell(i, j, A.get_cell(i, j) - DeltaV.get_cell(i - 1, j));
      DeltaV.set_cell(i, j, A.get_cell(i, j) - DeltaH.get_cell(i, j - 1));
      DeltaEp.set_cell(i, j, std::min({A.get_cell(i, j) + oH(j),
        DeltaEp.get_cell(i - 1, j) + eH(j)}) - DeltaH.get_cell(i, j - 1));
      DeltaFp.set_cell(i, j, std::min({A.get_cell(i, j) + oV(j),
        DeltaFp.get_cell(i, j - 1) + eV(j)}) - DeltaV.get_cell(i - 1, j));
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
        S.set_cell(i, j, S.get_cell(i, j - 1) + DeltaV.get_cell(i, j));
    auto infty = std::numeric_limits<double>::max();
    T ans = infty;
    for (int i = 1; i <= haplotype_size; i++) {
      ans = std::min(ans, S.get_cell(i, read_size));
    }
    return ans;
  }
  table::ProbabilityTable<T> get_A() { return A; }
  table::ProbabilityTable<T> get_DeltaH() { return DeltaH; }
  table::ProbabilityTable<T> get_DeltaV() { return DeltaV; }
  table::ProbabilityTable<T> get_DeltaEp() { return DeltaEp; }
  table::ProbabilityTable<T> get_DeltaFp() { return DeltaFp; }
};

template class SuzukiPairHMM<double>;
} // namespace pairhmm