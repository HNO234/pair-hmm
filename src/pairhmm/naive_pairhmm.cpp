#include "pairhmm/naive_pairhmm.hpp"
#include "table/probability_table.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <stdexcept>

namespace pairhmm {
template <typename T>
NaivePairHMM<T>::NaivePairHMM() {}
template <typename T>
NaivePairHMM<T>::NaivePairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
                    table::STRTable<T> gop_, table::STRTable<T> gcp_)
    : PairHMM<T>(haplotype_, read_, gop_, gcp_) {
  auto haplotype_size = this->haplotype.size();
  auto read_size = this->read.size();
  M = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  D = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
  I = table::ProbabilityTable<T>(haplotype_size + 1, read_size + 1);
}
template <typename T> biovoltron::Cigar NaivePairHMM<T>::get_cigar() {
  throw std::logic_error("Function not yet implemented");
  return biovoltron::Cigar{};
}
template <typename T> void NaivePairHMM<T>::run_alignment() {
  throw std::logic_error("Function not yet implemented");
}
template <typename T> table::ProbabilityTable<T> NaivePairHMM<T>::get_M() { return M; }
template <typename T> table::ProbabilityTable<T> NaivePairHMM<T>::get_D() { return D; }
template <typename T> table::ProbabilityTable<T> NaivePairHMM<T>::get_I() { return I; }
} // namespace pairhmm