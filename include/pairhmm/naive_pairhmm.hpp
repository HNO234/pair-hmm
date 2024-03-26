#pragma once

#include "pairhmm/pairhmm.hpp"
#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>

namespace pairhmm {
template <typename T> class NaivePairHMM : public PairHMM<T> {
private:
  table::ProbabilityTable<T> M, D, I;

public:
  // using PairHMM<T>::PairHMM;
  NaivePairHMM();
  NaivePairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
          table::STRTable<T> gop_, table::STRTable<T> gcp_);
  biovoltron::Cigar get_cigar();
  void run_alignment();
  table::ProbabilityTable<T> get_M();
  table::ProbabilityTable<T> get_D();
  table::ProbabilityTable<T> get_I();
};

template class NaivePairHMM<double>;
} // namespace pairhmm