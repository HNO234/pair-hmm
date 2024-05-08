#pragma once

#include <vector>

namespace table {
template <typename T> class Table {
protected:
  size_t rows, columns;
  std::vector<std::vector<T>> table;

public:
  Table(): rows(0), columns(0), table() {}
  Table(size_t rows_, size_t columns_): rows(rows_), columns(columns_), table(rows_, std::vector<T>(columns_)) {}
  void set_cell(size_t i, size_t j, T element) {
    table[i][j] = element;
  }
  T get_cell(size_t i, size_t j) const {
    return table[i][j];
  }
  void set_table(std::vector<std::vector<T>> &table_) {
    table = table_;
    rows = table.size();
    columns = table[0].size();
  }

  std::vector<std::vector<T>> get_table() {
    return table;
  }
  size_t get_rows() { return this->rows; }
  size_t get_columns() { return this->columns; }
};

template class Table<double>;
} // namespace table