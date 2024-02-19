#include "table/table.hpp"
#include <vector>

namespace table {
template <typename T> class Table {
  Table::Table(uint32_t rows_, uint32_t columns_) {}
  void Table::set_cell(uint32_t i, uint32_t j) {}
  auto Table::get_cell(uint32_t i, uint32_t j) {}
  auto Table::set_table(std::vector<std::vector<T>> &table_) {}
  auto Table::get_table() {}
}
} // namespace table