#include "pairhmm/naive_pairhmm.hpp"
#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>

namespace pairhmm {
template <typename T> biovoltron::Cigar NaivePairHMM<T>::get_cigar() {
  return biovoltron::Cigar{};
}
template <typename T> table::ProbabilityTable<T> NaivePairHMM<T>::get_S() { return S; }
template <typename T> table::ProbabilityTable<T> NaivePairHMM<T>::get_E() { return E; }
template <typename T> table::ProbabilityTable<T> NaivePairHMM<T>::get_F() { return F; }
} // namespace pairhmm