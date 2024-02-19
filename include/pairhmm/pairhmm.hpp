#pragma once

#include "table/STR_table.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <biovoltron/utility/istring.hpp>

namespace pairhmm {
template <typename T> class PairHMM {
private:
  biovoltron::istring haplotype, read;
  table::STRTable<T> gop, gcp;

public:
  PairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
          table::STRTable<T> gop_, table::STR_table<T> gcp_);
  virtual biovoltron::Cigar get_cigar() = 0;
}
} // namespace pairhmm