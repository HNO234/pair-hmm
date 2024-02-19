#pragma once

#include "table/table.hpp"

namespace table {
template <typename T> class ProbabilityTable : public Table {
public:
  ProbabilityTable shifted_down();
  ProbabilityTable shifted_right();
  ProbabilityTable operator-(const ProbabilityTable &right_table) const;
  friend bool is_close(ProbabilityTable &table1, ProbabilityTable &table2,
                       T tolerance);
}
} // namespace table