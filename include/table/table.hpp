#pragma once

#include <vector>

namespace table {
template <typename T> class Table {
protected:
  uint32_t rows, columns;
  std::vector<std::vector<T>> table;

public:
  Table(uint32_t rows_, uint32_t columns_);
  void set_cell(uint32_t i, uint32_t j);
  auto get_cell(uint32_t i, uint32_t j);
  auto set_table(std::vector<std::vector<T>> &table_);
  auto get_table();
}
} // namespace table