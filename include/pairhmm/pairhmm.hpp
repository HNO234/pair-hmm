#pragma once

#include "table/STR_table.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <biovoltron/utility/istring.hpp>
#include <tuple>

namespace pairhmm {
template <typename T> class PairHMM {
protected:
  biovoltron::istring haplotype, read;
  table::STRTable<T> gop, gcp;

public:
  PairHMM();
  PairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
          table::STRTable<T> gop_, table::STRTable<T> gcp_);
// vector of period, vector of repeat size
  std::tuple<std::vector<size_t>, std::vector<size_t>> get_read_best_repeat();
  virtual biovoltron::Cigar get_cigar() = 0;
  virtual void run_alignment() = 0;
};

template class PairHMM<double>;
} // namespace pairhmm