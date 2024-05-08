#pragma once

#include "table/table.hpp"

namespace table {
template <typename T> class ProbabilityTable : public Table<T> {
public:
  using Table<T>::Table;
  ProbabilityTable<T> shifted_down() {
  ProbabilityTable<T> ret(this->rows, this->columns);
  for (auto i = size_t{}; i < this->rows; i++)
    for (auto j = size_t{}; j < this->columns; j++)
      if (i != 0) {
        ret.set_cell(i, j, this->table[i - 1][j]);
      } else {
        ret.set_cell(i, j, 0);
      }
  return ret;
}
  ProbabilityTable<T> shifted_right() {
  ProbabilityTable<T> ret(this->rows, this->columns);
  for (auto i = size_t{}; i < this->rows; i++)
    for (auto j = size_t{}; j < this->columns; j++)
      if (j != 0) {
        ret.set_cell(i, j, this->table[i][j - 1]);
      } else {
        ret.set_cell(i, j, 0);
      }
  return ret;
}
  ProbabilityTable<T> operator-(const ProbabilityTable &right_table) const {
  ProbabilityTable<T> ret(this->rows, this->columns);
  for (auto i = size_t{}; i < this->rows; i++)
    for (auto j = size_t{}; j < this->columns; j++)
        ret.set_cell(i, j, this->table[i][j] - right_table.get_cell(i, j));
  return ret;
}
  ProbabilityTable<T> operator+(const ProbabilityTable &right_table) const {
  ProbabilityTable<T> ret(this->rows, this->columns);
  for (auto i = size_t{}; i < this->rows; i++)
    for (auto j = size_t{}; j < this->columns; j++)
        ret.set_cell(i, j, this->table[i][j] + right_table.get_cell(i, j));
  return ret;
}
  ProbabilityTable<T> operator*(const T &rhs) const {
  ProbabilityTable<T> ret(this->rows, this->columns);
  for (auto i = size_t{}; i < this->rows; i++)
    for (auto j = size_t{}; j < this->columns; j++)
        ret.set_cell(i, j, this->table[i][j] * rhs);
  return ret;
}
};

template <typename T>
bool is_close(ProbabilityTable<T> &table1, ProbabilityTable<T> &table2,
              T tolerance) {
  if (table1.get_rows() != table2.get_rows())
    return false;
  if (table1.get_columns() != table2.get_columns())
    return false;
  for (auto i = size_t{1}; i < table1.get_rows(); i++)
    for (auto j = size_t{1}; j < table1.get_columns(); j++)
      if (abs(table1.get_cell(i, j) - table2.get_cell(i, j)) > tolerance)
        return false;
  return true;
}

// template class ProbabilityTable<double>;
// template <> bool is_close(ProbabilityTable<double> &table1, ProbabilityTable<double> &table2,
//               double tolerance);
} // namespace table