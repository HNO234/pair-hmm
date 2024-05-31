#pragma once

#include "pairhmm/pairhmm.hpp"
#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include "pairhmm/nw_pairhmm.hpp"
#include "table/probability_table.hpp"
#include "utils/constant.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>

namespace pairhmm {
template <typename T> class NWPairHMM : public PairHMM<T> {
private:
  table::ProbabilityTable<T> S, E, F;
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
  NWPairHMM() {}
  NWPairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
                    table::STRTable<T> gop_, table::STRTable<T> gcp_,
                    std::vector<int> score_)
    : PairHMM<T>(haplotype_, read_, gop_, gcp_) {
  score = score_;
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  S = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  E = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  F = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  auto [period_read_, repeat_read_] = this->get_read_best_repeat();
  period_read = period_read_;
  repeat_read = repeat_read_;
}
  biovoltron::Cigar get_cigar() {
  throw std::logic_error("Function not yet implemented");
  return biovoltron::Cigar{};
}

  void run_alignment() {
  // throw std::logic_error("Function not yet implemented");
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  auto infty = std::numeric_limits<double>::max();
  for (auto i = size_t{}; i <= haplotype_size; i++) {
    S.set_cell(i, 0, initial_value());
    E.set_cell(i, 0, infty);
    F.set_cell(i, 0, infty);
  }
  for (auto j = size_t{1}; j <= read_size; j++) {
    S.set_cell(0, j, infty);
    E.set_cell(0, j, infty);
    F.set_cell(0, j, infty);
  }

  for (auto i = size_t{1}; i <= haplotype_size; i++)
    for (auto j = size_t{1}; j <= read_size; j++) {
      // Update M
      auto S_match = S.get_cell(i - 1, j - 1) + s(i, j);
      auto S_deletion = E.get_cell(i - 1, j) + cH(j);
      auto S_insertion = F.get_cell(i, j - 1) + cV(j);
      S.set_cell(i, j, std::min({S_match, S_deletion, S_insertion}));
      // Update D
      auto E_gap_open = S.get_cell(i, j) + oH(j);
      auto E_gap_continue = E.get_cell(i - 1, j) + eH(j);
      E.set_cell(i, j, std::min({E_gap_open, E_gap_continue}));
      // Update I
      auto F_gap_open = S.get_cell(i, j) + oV(j);
      auto F_gap_continue = F.get_cell(i, j - 1) + eV(j);
      F.set_cell(i, j, std::min(F_gap_open, F_gap_continue));
  }
}

  T get_align_score() {
    auto haplotype_size = this->haplotype.size();
    auto read_size = this->read.size();
    run_alignment();

    auto infty = std::numeric_limits<double>::max();
    T ans = infty;
    for (int i = 0; i <= haplotype_size; i++) {
      ans = std::min(ans, S.get_cell(i, read_size));
    }
    if (ans > 1000000) {
      for (int i = 1; i<= haplotype_size; i++)
        for (int j = 1; j<= read_size; j++) {
          std::cout << i << ' '<< j << ' ' << S.get_cell(i, j) << ' '
          << S.get_cell(i - 1, j - 1) << ' ' << s(i, j) << ' '
          << E.get_cell(i - 1, j) << ' ' << cH(j) << ' ' <<
          F.get_cell(i, j - 1) << ' ' << cV(j) << ' ' << period_read[j + 1] << ' ' << repeat_read[j + 1] << std::endl;
        }
    }
    return ans;
  }
  table::ProbabilityTable<T> get_S() { return S; }
  table::ProbabilityTable<T> get_E() { return E; }
  table::ProbabilityTable<T> get_F() { return F; }
};

template class NWPairHMM<double>;
} // namespace pairhmm