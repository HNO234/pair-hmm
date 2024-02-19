#pragma once

#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>

namespace pairhmm {
template <typename T> class NaivePairHMM: {
private:
  table::ProbablityTable S, E, F;

public:
  biovoltron::Cigar get_cigar();
  auto get_S();
  auto get_E();
  auto get_F();
}
} // namespace pairhmm