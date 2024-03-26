#include "pairhmm/pairhmm.hpp"
#include "table/STR_table.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <biovoltron/utility/istring.hpp>

namespace pairhmm {
template <typename T>
PairHMM<T>::PairHMM() : haplotype(), read(), gop(), gcp() {}
template <typename T>
PairHMM<T>::PairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
                    table::STRTable<T> gop_, table::STRTable<T> gcp_)
    : haplotype(haplotype_), read(read_), gop(gop_), gcp(gcp_) {}
template <typename T>
std::tuple<std::vector<size_t>, std::vector<size_t>>
PairHMM<T>::get_read_best_repeat() {
  std::vector<size_t> best_period(haplotype.size()), best_size(read.size());
  for (auto &period: best_period)
    period = rand() % 8;
  for (auto &size: best_size)
    size = rand() % 20;
  return std::make_tuple(best_period, best_size);
}
} // namespace pairhmm