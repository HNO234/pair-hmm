#include "pairhmm/nw_pairhmm.hpp"
#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <stdexcept>

namespace pairhmm {
template <typename T> biovoltron::Cigar NWPairHMM<T>::get_cigar() {
  throw std::logic_error("Function not yet implemented");
  return biovoltron::Cigar{};
}
template <typename T> void NWPairHMM<T>::run_alignment() {
  throw std::logic_error("Function not yet implemented");
}
template <typename T> table::ProbabilityTable<T> NWPairHMM<T>::get_S() { return S; }
template <typename T> table::ProbabilityTable<T> NWPairHMM<T>::get_E() { return E; }
template <typename T> table::ProbabilityTable<T> NWPairHMM<T>::get_F() { return F; }
} // namespace pairhmm