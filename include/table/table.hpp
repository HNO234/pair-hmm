#pragma once

#include <vector>

namespace table {
template <typename T> class Table {
protected:
  size_t rows, columns;
  std::vector<std::vector<T>> table;

public:
  Table();
  Table(size_t rows_, size_t columns_);
  void set_cell(size_t i, size_t j, T element);
  T get_cell(size_t i, size_t j) const;
  void set_table(std::vector<std::vector<T>> &table_);
  std::vector<std::vector<T>> get_table();
  size_t get_rows();
  size_t get_columns();
};

template class Table<double>;
} // namespace table