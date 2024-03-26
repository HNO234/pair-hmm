#pragma once

#include "pairhmm/pairhmm.hpp"
#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>

namespace pairhmm {
template <typename T> class SuzukiPairHMM : public PairHMM<T> {
private:
  table::ProbabilityTable<T> A, DeltaH, DeltaV, DeltaEp, DeltaFp;
  std::vector<size_t> period_read, repeat_read;
  T s(size_t i, size_t j);
  T oH(size_t j);
  T oV(size_t j);
  T eH(size_t j);
  T eV(size_t j);
  T cH(size_t j);
  T cV(size_t j);

public:
  // using PairHMM<T>::PairHMM;
  SuzukiPairHMM();
  SuzukiPairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
                    table::STRTable<T> gop_, table::STRTable<T> gcp_);
  biovoltron::Cigar get_cigar();
  void run_alignment();
  table::ProbabilityTable<T> get_A();
  table::ProbabilityTable<T> get_DeltaH();
  table::ProbabilityTable<T> get_DeltaV();
  table::ProbabilityTable<T> get_DeltaEp();
  table::ProbabilityTable<T> get_DeltaFp();
};

template class SuzukiPairHMM<double>;
} // namespace pairhmm