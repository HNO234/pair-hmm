#include "pairhmm/nw_pairhmm.hpp"
#include "table/probability_table.hpp"
#include "utils/constant.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pairhmm {
template <typename T>
NWPairHMM<T>::NWPairHMM() {}
template <typename T>
NWPairHMM<T>::NWPairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
                    table::STRTable<T> gop_, table::STRTable<T> gcp_)
    : PairHMM<T>(haplotype_, read_, gop_, gcp_) {
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  S = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  E = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  F = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  auto [period_read_, repeat_read_] = this->get_read_best_repeat();
  period_read = period_read_;
  repeat_read = repeat_read_;
}
template <typename T>
T NWPairHMM<T>::s(size_t i, size_t j) {
  auto match_prob = (i && j ? base_match_prob[this->haplotype[i - 1]][this->read[j - 1]] : 0);
  auto best_period_j_minus_1 = (j ? period_read[j - 1] : 0);
  auto best_repeat_j_minus_1 = (j ? repeat_read[j - 1] : 0);
  auto delta_j_minus_1_real = (j ? std::pow(10, -0.1 * this->gop.get_cell(best_period_j_minus_1 - 1, best_repeat_j_minus_1 - 1)) : 0);
  auto phred_1_minus_2_times_delta_j_minus_1 = -10 * std::log10(1 - 2 * delta_j_minus_1_real);
  return phred_1_minus_2_times_delta_j_minus_1 + match_prob;
}
template <typename T>
T NWPairHMM<T>::oH(size_t j) {
  auto best_period_j = period_read[j];
  auto best_repeat_j = repeat_read[j];
  auto delta_j = this->gop.get_cell(best_period_j -  1, best_repeat_j - 1);
  return delta_j;
}
template <typename T>
T NWPairHMM<T>::oV(size_t j) {
  auto best_period_j = period_read[j];
  auto best_repeat_j = repeat_read[j];
  auto delta_j = this->gop.get_cell(best_period_j -  1, best_repeat_j - 1);
  return delta_j;
}
template <typename T>
T NWPairHMM<T>::eH(size_t j) {
  auto best_period_j = period_read[j];
  auto best_repeat_j = repeat_read[j];
  auto epsilon_j = this->gcp.get_cell(best_period_j - 1, best_repeat_j - 1);
  return epsilon_j + mismatch_prob;
}
template <typename T>
T NWPairHMM<T>::eV(size_t j) {
  auto best_period_j_minus_1 = (j ? period_read[j - 1] : 0);
  auto best_repeat_j_minus_1 = (j ? repeat_read[j - 1] : 0);
  auto epislon_j_minus_1 = (j ? this->gcp.get_cell(best_period_j_minus_1 - 1, best_repeat_j_minus_1 - 1) : 0);
  return epislon_j_minus_1 + mismatch_prob;
}
template <typename T>
T NWPairHMM<T>::cH(size_t j) {
  auto best_period_j = period_read[j];
  auto best_repeat_j = repeat_read[j];
  auto epsilon_j = this->gcp.get_cell(best_period_j - 1, best_repeat_j - 1);
  auto phred_1_minus_epislon_j = -10 * std::log10(1 - std::pow(10, -0.1 * epsilon_j));
  return phred_1_minus_epislon_j + mismatch_prob;
}
template <typename T>
T NWPairHMM<T>::cV(size_t j) {
  auto best_period_j_minus_1 = (j ? period_read[j - 1] : 0);
  auto best_repeat_j_minus_1 = (j ? repeat_read[j - 1] : 0);
  auto epislon_j_minus_1 = (j ? this->gcp.get_cell(best_period_j_minus_1 - 1, best_repeat_j_minus_1 - 1) : 0);
  auto phred_1_minus_epislon_j_minus_1 = -10 * std::log10(1 - std::pow(10, -0.1 * epislon_j_minus_1));
  return phred_1_minus_epislon_j_minus_1 + mismatch_prob;
}
template <typename T> biovoltron::Cigar NWPairHMM<T>::get_cigar() {
  auto is_close = [](T x, T y) { return std::abs(x - y) < 1e-9; };
  // preparation
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  auto infty = std::numeric_limits<double>::lowest();
  auto i = haplotype_size, j = read_size;
  auto state = 'M';
  auto return_cigar = biovoltron::Cigar{};

  while (i || j) {
    if (state == 'M') { // M
      auto S_match = (i && j ? S.get_cell(i - 1, j - 1) + s(i, j) : infty);
      auto S_deletion = (i ? E.get_cell(i - 1, j) + cH(j) : infty);
      auto S_insertion = (j ? F.get_cell(i, j - 1) + cV(j) : infty);
      if (is_close(S_match, S.get_cell(i, j))) {
        i--, j--;
        state = 'M';
        return_cigar.emplace_back(1, 'M');
      } else if (is_close(S_deletion, S.get_cell(i, j))) {
        i--;
        state = 'D';
        return_cigar.emplace_back(1, 'D');
      } else if (is_close(S_insertion, S.get_cell(i, j))) {
        j--;
        state = 'I';
        return_cigar.emplace_back(1, 'I');
      }
    } else if (state == 'I') { // I
      auto E_gap_open = S.get_cell(i, j) + oH(j);
      auto E_gap_continue = (i ? E.get_cell(i - 1, j) + eH(j) : infty);
    if (is_close(E_gap_open, E.get_cell(i, j))) {
        state = 'M';
      } else if (is_close(E_gap_continue, E.get_cell(i, j))) {
        j--;
        state = 'I';
        return_cigar.emplace_back(1, 'I');
      }
    } else if (state == 'D') { // D
      auto F_gap_open = S.get_cell(i, j) + oV(j);
      auto F_gap_continue = (j ? F.get_cell(i, j - 1) + eV(j) : infty);
      if (is_close(F_gap_open, F.get_cell(i, j))) {
        state = 'M';
      } else if (is_close(F_gap_continue, F.get_cell(i, j))) {
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
template <typename T> void NWPairHMM<T>::run_alignment() {
  // throw std::logic_error("Function not yet implemented");
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  auto infty = std::numeric_limits<double>::max();
  for (auto i = size_t{}; i <= haplotype_size; i++)
    for (auto j = size_t{}; j <= read_size; j++) {
      if (i == 0 && j == 0) {
        S.set_cell(i, j, 0);
        E.set_cell(i, j, infty);
        F.set_cell(i, j, infty);
        continue;
      }
      // Update M
      auto S_match = (i && j ? S.get_cell(i - 1, j - 1) + s(i, j) : infty);
      auto S_deletion = (i ? E.get_cell(i - 1, j) + cH(j) : infty);
      auto S_insertion = (j ? F.get_cell(i, j - 1) + cV(j) : infty);
      S.set_cell(i, j, std::min({S_match, S_deletion, S_insertion}));
      // Update D
      auto E_gap_open = S.get_cell(i, j) + oH(j);
      auto E_gap_continue = (i ? E.get_cell(i - 1, j) + eH(j) : infty);
      E.set_cell(i, j, std::min({E_gap_open, E_gap_continue}));
      // Update I
      auto F_gap_open = S.get_cell(i, j) + oV(j);
      auto F_gap_continue = (j ? F.get_cell(i, j - 1) + eV(j) : infty);
      F.set_cell(i, j, std::min(F_gap_open, F_gap_continue));
  }
}
template <typename T> table::ProbabilityTable<T> NWPairHMM<T>::get_S() { return S; }
template <typename T> table::ProbabilityTable<T> NWPairHMM<T>::get_E() { return E; }
template <typename T> table::ProbabilityTable<T> NWPairHMM<T>::get_F() { return F; }
} // namespace pairhmm