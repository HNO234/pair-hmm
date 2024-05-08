#pragma once

#include "pairhmm/pairhmm.hpp"
#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include "pairhmm/naive_pairhmm.hpp"
#include "table/probability_table.hpp"
#include "utils/constant.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>

namespace pairhmm {
template <typename T> class NaivePairHMM : public PairHMM<T> {
private:
  table::ProbabilityTable<T> M, D, I;

public:
  // using PairHMM<T>::PairHMM;
  NaivePairHMM() {}
  NaivePairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
                    table::STRTable<T> gop_, table::STRTable<T> gcp_)
    : PairHMM<T>(haplotype_, read_, gop_, gcp_) {
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  M = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  D = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  I = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
}
  biovoltron::Cigar get_cigar() {
  // throw std::logic_error("Function not yet implemented");
  auto is_close = [](T x, T y) { return std::abs(x - y) < 1e-9; };
  // preparation
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  auto [period_read, repeat_read] = this->get_read_best_repeat();
  auto neg_infty = std::numeric_limits<double>::lowest();
  auto i = haplotype_size, j = read_size;
  auto state = 'M';
  auto return_cigar = biovoltron::Cigar{};

  while (i || j) {
    auto best_period_j = period_read[j];
      auto best_period_j_minus_1 = (j ? period_read[j - 1] : 0);
      auto best_repeat_j = repeat_read[j];
      auto best_repeat_j_minus_1 = (j ? repeat_read[j - 1] : 0);
      auto delta_j = std::pow(10, -0.1 * this->gop.get_cell(best_period_j -  1, best_repeat_j - 1));
      auto delta_j_minus_1 = (j ? std::pow(10, -0.1 * this->gop.get_cell(best_period_j_minus_1 - 1, best_repeat_j_minus_1 - 1)) : 0);
      auto epsilon_j = std::pow(10, -0.1 * this->gcp.get_cell(best_period_j - 1, best_repeat_j - 1));
      auto epislon_j_minus_1 = (j ? std::pow(10, -0.1 * this->gcp.get_cell(best_period_j_minus_1 - 1, best_repeat_j_minus_1 - 1)) : 0);
    if (state == 'M') { // M
      auto match_prob = (i && j ? base_match_prob[this->haplotype[i - 1]][this->read[j - 1]] : 0);
      auto M_match = (i && j ? M.get_cell(i - 1, j - 1) + std::log10(1 - 2 * delta_j_minus_1) - 0.1 * match_prob : neg_infty);
      auto M_deletion = (i ? D.get_cell(i - 1, j) + std::log10(1 - epsilon_j) - 0.1 * mismatch_prob : neg_infty);
      auto M_insertion = (j ? I.get_cell(i, j - 1) + std::log10(1 - epislon_j_minus_1) - 0.1 * mismatch_prob : neg_infty);
      if (is_close(M_match, M.get_cell(i, j))) {
        i--, j--;
        state = 'M';
        return_cigar.emplace_back(1, 'M');
      } else if (is_close(M_deletion, M.get_cell(i, j))) {
        i--;
        state = 'D';
        return_cigar.emplace_back(1, 'D');
      } else if (is_close(M_insertion, M.get_cell(i, j))) {
        j--;
        state = 'I';
        return_cigar.emplace_back(1, 'I');
      }
    } else if (state == 'I') { // I
      auto I_gap_open = M.get_cell(i, j) + std::log10(delta_j);
      auto I_gap_continue = (j ? I.get_cell(i, j - 1) + std::log10(epislon_j_minus_1) - 0.1 * mismatch_prob : neg_infty);
      if (is_close(I_gap_open, I.get_cell(i, j))) {
        state = 'M';
      } else if (is_close(I_gap_continue, I.get_cell(i, j))) {
        j--;
        state = 'I';
        return_cigar.emplace_back(1, 'I');
      }
    } else if (state == 'D') { // D
      auto D_gap_open = M.get_cell(i, j) + std::log10(delta_j);
      auto D_gap_continue = (i ? D.get_cell(i - 1, j) + std::log10(epsilon_j) - 0.1 * mismatch_prob : neg_infty);
      if (is_close(D_gap_open, D.get_cell(i, j))) {
        state = 'M';
      } else if (is_close(D_gap_continue, D.get_cell(i, j))) {
        i--;
        state = 'D';
        return_cigar.emplace_back(1, 'D');
      }
    }
  }
  return_cigar.reverse();
  return_cigar.compact();
  return return_cigar;
}
  void run_alignment() {
  // throw std::logic_error("Function not yet implemented");
  // preparation
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  auto [period_read, repeat_read] = this->get_read_best_repeat();
  auto neg_infty = std::numeric_limits<double>::lowest();
  for (auto i = size_t{}; i <= haplotype_size; i++)
    for (auto j = size_t{}; j <= read_size; j++) {
      auto best_period_j = period_read[j];
      auto best_period_j_minus_1 = (j ? period_read[j - 1] : 0);
      auto best_repeat_j = repeat_read[j];
      auto best_repeat_j_minus_1 = (j ? repeat_read[j - 1] : 0);
      auto delta_j = std::pow(10, -0.1 * this->gop.get_cell(best_period_j -  1, best_repeat_j - 1));
      auto delta_j_minus_1 = (j ? std::pow(10, -0.1 * this->gop.get_cell(best_period_j_minus_1 - 1, best_repeat_j_minus_1 - 1)) : 0);
      auto epsilon_j = std::pow(10, -0.1 * this->gcp.get_cell(best_period_j - 1, best_repeat_j - 1));
      auto epislon_j_minus_1 = (j ? std::pow(10, -0.1 * this->gcp.get_cell(best_period_j_minus_1 - 1, best_repeat_j_minus_1 - 1)) : 0);
      // Update M
      if (i == 0 && j == 0) {
        M.set_cell(i, j, 0);
      } else {
        auto match_prob = (i && j ? base_match_prob[this->haplotype[i - 1]][this->read[j - 1]] : 0);
        auto M_match = (i && j ? M.get_cell(i - 1, j - 1) + std::log10(1 - 2 * delta_j_minus_1) - 0.1 * match_prob : neg_infty);
        auto M_deletion = (i ? D.get_cell(i - 1, j) + std::log10(1 - epsilon_j) - 0.1 * mismatch_prob : neg_infty);
        auto M_insertion = (j ? I.get_cell(i, j - 1) + std::log10(1 - epislon_j_minus_1) - 0.1 * mismatch_prob : neg_infty);
        M.set_cell(i, j, std::max({M_match, M_deletion, M_insertion}));
      }
      // Update D
      auto D_gap_open = M.get_cell(i, j) + std::log10(delta_j);
      auto D_gap_continue = (i ? D.get_cell(i - 1, j) + std::log10(epsilon_j) - 0.1 * mismatch_prob : neg_infty);
      D.set_cell(i, j, std::max({D_gap_open, D_gap_continue}));
      // Update I
      auto I_gap_open = M.get_cell(i, j) + std::log10(delta_j);
      auto I_gap_continue = (j ? I.get_cell(i, j - 1) + std::log10(epislon_j_minus_1) - 0.1 * mismatch_prob : neg_infty);
      I.set_cell(i, j, std::max(I_gap_open, I_gap_continue));
  }
}
  table::ProbabilityTable<T> get_M() { return M; }
  table::ProbabilityTable<T> get_D() { return D; }
  table::ProbabilityTable<T> get_I() { return I; }
};

template class NaivePairHMM<double>;
} // namespace pairhmm