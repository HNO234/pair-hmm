#pragma once

#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>

namespace pairhmm {
template <typename T> class SuzukiPairHMM: {
private:
  table::ProbablityTable A, DeltaH, DeltaV, DeltaEp, DeltaFp;

public:
  biovoltron::Cigar get_cigar();
  auto get_A();
  auto get_DeltaH();
  auto get_DeltaV();
  auto get_DeltaEp();
  auto get_DeltaFp();
}
} // namespace pairhmm