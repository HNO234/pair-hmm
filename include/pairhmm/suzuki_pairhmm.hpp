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
    else if (j == this->read.size()) {
      return 45;
    } else {
      auto best_period_j = period_read[j];
      auto best_repeat_j = repeat_read[j];
      auto delta_j = this->gop.get_cell(best_period_j -  1, best_repeat_j - 1);
      return std::max((T)40, delta_j);
    }
  }

  T get_gcp_phred (size_t j) {
    if (j == 0)
      return 0;
    else if (j == this->read.size()) {
      return 10;
    } else {
      auto best_period_j = period_read[j];
      auto best_repeat_j = repeat_read[j];
      auto epislon_j = this->gcp.get_cell(best_period_j -  1, best_repeat_j - 1);
      return epislon_j;
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
  auto gop_phred = get_gop_phred(j - 1);
  auto delta_j_minus_1_real = (j ? std::pow(10, -0.1 * gop_phred) : 0);
  auto phred_1_minus_2_times_delta_j_minus_1 = -10 * std::log10(1 - 2 * delta_j_minus_1_real);
  return phred_1_minus_2_times_delta_j_minus_1 + match_prob;
}
  T oH(size_t j) {
    return get_gop_phred(j);
}
  T oV(size_t j) {
    return get_gop_phred(j);
}
  T eH(size_t j) {
  return get_gcp_phred(j) + misalign();
}
  T eV(size_t j) {
  return get_gcp_phred(j - 1) + misalign();
}
  T cH(size_t j) {
  auto epsilon_j = get_gcp_phred(j);
  auto phred_1_minus_epislon_j = -10 * std::log10(1 - std::pow(10, -0.1 * epsilon_j));
  return phred_1_minus_epislon_j + misalign();
}
  T cV(size_t j) {
  auto epislon_j_minus_1 = get_gcp_phred(j - 1);
  auto phred_1_minus_epislon_j_minus_1 = -10 * std::log10(1 - std::pow(10, -0.1 * epislon_j_minus_1));
  return phred_1_minus_epislon_j_minus_1 + misalign();
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
  // throw std::logic_error("Function not yet implemented");
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  auto S = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  auto E = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  auto F = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  for (auto j = size_t{1}; j <= read_size; j++) {
    S.set_cell(0, j, S.get_cell(0, j - 1) + DeltaV.get_cell(0, j));
  }
  for (auto i = size_t{1}; i <= haplotype_size; i++)
    for (auto j = size_t{}; j <= read_size; j++)
      S.set_cell(i, j, S.get_cell(i - 1, j) + DeltaH.get_cell(i, j));
  for (auto i = size_t{}; i <= haplotype_size; i++)
    for (auto j = size_t{1}; j <= read_size; j++)
      E.set_cell(i, j, DeltaEp.get_cell(i, j) + S.get_cell(i, j - 1));
  for (auto i = size_t{1}; i <= haplotype_size; i++)
    for (auto j = size_t{}; j <= read_size; j++)
      F.set_cell(i, j, DeltaFp.get_cell(i, j) + S.get_cell(i - 1, j));
  
  // auto neg_infty = std::numeric_limits<double>::lowest();
  // auto i = haplotype_size, j = read_size;
  // auto state = 'M';
  // auto return_cigar = biovoltron::Cigar{};

  // while (i || j) {
  //   if (i == 0) {
  //     j--;
  //     state = 'I';
  //     return_cigar.emplace_back(1, 'I');
  //     continue;
  //   }
  //   if (j == 0) {
  //     i--;
  //     state = 'D';
  //     return_cigar.emplace_back(1, 'D');
  //     continue;
  //   }
  //   auto best_period_j = period_read[j];
  //     auto best_period_j_minus_1 = (j ? period_read[j - 1] : 0);
  //     auto best_repeat_j = repeat_read[j];
  //     auto best_repeat_j_minus_1 = (j ? repeat_read[j - 1] : 0);
  //     auto delta_j = std::pow(10, -0.1 * this->gop.get_cell(best_period_j -  1, best_repeat_j - 1));
  //     auto delta_j_minus_1 = (j ? std::pow(10, -0.1 * this->gop.get_cell(best_period_j_minus_1 - 1, best_repeat_j_minus_1 - 1)) : 0);
  //     auto epsilon_j = std::pow(10, -0.1 * this->gcp.get_cell(best_period_j - 1, best_repeat_j - 1));
  //     auto epislon_j_minus_1 = (j ? std::pow(10, -0.1 * this->gcp.get_cell(best_period_j_minus_1 - 1, best_repeat_j_minus_1 - 1)) : 0);
  //   if (state == 'M') { // M
  //     auto match_prob = (i && j ? base_match_prob[this->haplotype[i - 1]][this->read[j - 1]] : 0);
  //     auto M_match = (i && j ? M.get_cell(i - 1, j - 1) + std::log(1 - 2 * delta_j_minus_1) - 0.1 * match_prob : neg_infty);
  //     auto M_deletion = (i ? D.get_cell(i - 1, j) + std::log(1 - epsilon_j) - 0.1 * misalign() : neg_infty);
  //     auto M_insertion = (j ? I.get_cell(i, j - 1) + std::log(1 - epislon_j_minus_1) - 0.1 * misalign() : neg_infty);
  //     if (is_close(M_match, M.get_cell(i, j))) {
  //       i--, j--;
  //       state = 'M';
  //       return_cigar.emplace_back(1, 'M');
  //     } else if (is_close(M_deletion, M.get_cell(i, j))) {
  //       i--;
  //       state = 'D';
  //       return_cigar.emplace_back(1, 'D');
  //     } else if (is_close(M_insertion, M.get_cell(i, j))) {
  //       j--;
  //       state = 'I';
  //       return_cigar.emplace_back(1, 'I');
  //     }
  //   } else if (state == 'I') { // I
  //     auto I_gap_open = M.get_cell(i, j) + std::log(delta_j);
  //     auto I_gap_continue = (j ? I.get_cell(i, j - 1) + std::log(epislon_j_minus_1) - 0.1 * misalign() : neg_infty);
  //     if (is_close(I_gap_open, I.get_cell(i, j))) {
  //       state = 'M';
  //     } else if (is_close(I_gap_continue, I.get_cell(i, j))) {
  //       j--;
  //       state = 'I';
  //       return_cigar.emplace_back(1, 'I');
  //     }
  //   } else if (state == 'D') { // D
  //     auto D_gap_open = M.get_cell(i, j) + std::log(delta_j);
  //     auto D_gap_continue = (i ? D.get_cell(i - 1, j) + std::log(epsilon_j) - 0.1 * misalign() : neg_infty);
  //     if (is_close(D_gap_open, D.get_cell(i, j))) {
  //       state = 'M';
  //     } else if (is_close(D_gap_continue, D.get_cell(i, j))) {
  //       i--;
  //       state = 'D';
  //       return_cigar.emplace_back(1, 'D');
  //     }
  //   }
  // }
  // return_cigar.reverse();
  // return_cigar.compact();
  // return return_cigar;
  return biovoltron::Cigar{};
}

  void run_alignment() {
  // throw std::logic_error("Function not yet implemented");
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  auto neg_infty = std::numeric_limits<double>::lowest();
  // calculate values for j = 0
  auto last_gap_penalty_j_0 = T{};
  auto accumulated_gap_penalty_j_0 = T{};
  for (auto i = size_t{1}; i <= haplotype_size; i++) {
    // Calculate DeltaH
    if (i == 1)
      accumulated_gap_penalty_j_0 = (oH(0) + cH(0));
    else
      accumulated_gap_penalty_j_0 += (-cH(0) + eH(0) + cH(0));
    DeltaH.set_cell(i, 0, accumulated_gap_penalty_j_0 - last_gap_penalty_j_0);
    last_gap_penalty_j_0 = accumulated_gap_penalty_j_0;
    // Calculate DeltaFp
    DeltaFp.set_cell(i, 0, DeltaH.get_cell(i, 0) + oV(0));
  }
  // calculate values for i = 0
  auto last_gap_penalty_i_0 = T{};
  auto accumulated_gap_penalty_i_0 = T{};
  for (auto j = size_t{1}; j <= read_size; j++) {
    // Calculate DeltaV
    if (j == 1)
      accumulated_gap_penalty_i_0 = (oV(j - 1) + cV(j));
    else
      accumulated_gap_penalty_i_0 += (-cV(j - 1) + eV(j - 1) + cV(j));
    DeltaV.set_cell(0, j, accumulated_gap_penalty_i_0 - last_gap_penalty_i_0);
    last_gap_penalty_i_0 = accumulated_gap_penalty_i_0;
    // Calculate DeltaEp
    DeltaEp.set_cell(0, j, DeltaV.get_cell(0, j) + oH(j));
    int i = 0;
  }
  for (auto i = size_t{1}; i <= haplotype_size; i++)
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
    for (auto j = size_t{1}; j <= read_size; j++) {
      S.set_cell(0, j, S.get_cell(0, j - 1) + DeltaV.get_cell(0, j));
    }
    for (auto i = size_t{1}; i <= haplotype_size; i++)
      for (auto j = size_t{}; j <= read_size; j++)
        S.set_cell(i, j, S.get_cell(i - 1, j) + DeltaH.get_cell(i, j));
    return S.get_cell(haplotype_size, read_size);
    // return A.get_cell(haplotype_size, read_size);
  }
  table::ProbabilityTable<T> get_A() { return A; }
  table::ProbabilityTable<T> get_DeltaH() { return DeltaH; }
  table::ProbabilityTable<T> get_DeltaV() { return DeltaV; }
  table::ProbabilityTable<T> get_DeltaEp() { return DeltaEp; }
  table::ProbabilityTable<T> get_DeltaFp() { return DeltaFp; }
};

template class SuzukiPairHMM<double>;
} // namespace pairhmm