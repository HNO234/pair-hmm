#pragma once

#include "pairhmm/pairhmm.hpp"
#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>

namespace pairhmm {
template <typename T> class SuzukiPairHMM : public PairHMM<T> {
private:
  table::ProbabilityTable<T> A, DeltaH, DeltaV, DeltaEp, DeltaFp;

public:
  using PairHMM<T>::PairHMM;
  biovoltron::Cigar get_cigar();
  table::ProbabilityTable<T> get_A();
  table::ProbabilityTable<T> get_DeltaH();
  table::ProbabilityTable<T> get_DeltaV();
  table::ProbabilityTable<T> get_DeltaEp();
  table::ProbabilityTable<T> get_DeltaFp();
};

template class SuzukiPairHMM<double>;
} // namespace pairhmm