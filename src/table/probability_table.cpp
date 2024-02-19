#include "table/probability_table.hpp"

namespace table {
template <typename T> class ProbabilityTable : public Table {
  ProbabilityTable ProbabilityTable::shifted_down() {}
  ProbabilityTable ProbabilityTable::shifted_right() {}
  ProbabilityTable::operator-(const ProbabilityTable &right_table) const {}
  bool is_close(ProbabilityTable &table1, ProbabilityTable &table2,
                T tolerance) {}
}
} // namespace table