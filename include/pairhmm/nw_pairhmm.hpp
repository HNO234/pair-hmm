#pragma once

#include "pairhmm/pairhmm.hpp"
#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>

namespace pairhmm {
template <typename T> class NWPairHMM : public PairHMM<T> {
private:
  table::ProbabilityTable<T> S, E, F;

public:
  using PairHMM<T>::PairHMM;
  biovoltron::Cigar get_cigar();
  void run_alignment();
  table::ProbabilityTable<T> get_S();
  table::ProbabilityTable<T> get_E();
  table::ProbabilityTable<T> get_F();
};

template class NWPairHMM<double>;
} // namespace pairhmm