#include "pairhmm/suzuki_pairhmm.hpp"
#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>

namespace pairhmm {
template <typename T> biovoltron::Cigar SuzukiPairHMM<T>::get_cigar() {
  return biovoltron::Cigar{};
}
template <typename T> table::ProbabilityTable<T> SuzukiPairHMM<T>::get_A() { return A; }
template <typename T> table::ProbabilityTable<T> SuzukiPairHMM<T>::get_DeltaH() { return DeltaH; }
template <typename T> table::ProbabilityTable<T> SuzukiPairHMM<T>::get_DeltaV() { return DeltaV; }
template <typename T> table::ProbabilityTable<T> SuzukiPairHMM<T>::get_DeltaEp() { return DeltaEp; }
template <typename T> table::ProbabilityTable<T> SuzukiPairHMM<T>::get_DeltaFp() { return DeltaFp; }
} // namespace pairhmm