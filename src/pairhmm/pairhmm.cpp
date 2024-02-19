#include "pairhmm/pairhmm.hpp"
#include "table/STR_table.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <biovoltron/utility/istring.hpp>

namespace pairhmm {
PairHMM::PairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
          table::STRTable<T> gop_, table::STR_table<T> gcp_) {}
} // namespace pairhmm