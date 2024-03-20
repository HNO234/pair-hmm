#pragma once

#include "table/table.hpp"

namespace table {
template <typename T> class ProbabilityTable : public Table<T> {
public:
  using Table<T>::Table;
  ProbabilityTable<T> shifted_down();
  ProbabilityTable<T> shifted_right();
  ProbabilityTable<T> operator-(const ProbabilityTable<T> &right_table) const;
  ProbabilityTable<T> operator+(const ProbabilityTable<T> &right_table) const;
};

template <typename T>
bool is_close(ProbabilityTable<T> &table1, ProbabilityTable<T> &table2,
              T tolerance);

template class ProbabilityTable<double>;
template <> bool is_close(ProbabilityTable<double> &table1, ProbabilityTable<double> &table2,
              double tolerance);
} // namespace table