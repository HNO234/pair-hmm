#include "table/table.hpp"
#include <vector>

namespace table {
template <typename T> Table<T>::Table() : rows(0), columns(0), table() {}
template <typename T>
Table<T>::Table(size_t rows_, size_t columns_)
    : rows(rows_), columns(columns_), table(rows_, std::vector<T>(columns_)) {}
template <typename T> void Table<T>::set_cell(size_t i, size_t j, T element) {}
template <typename T> T Table<T>::get_cell(size_t i, size_t j) {
  return table[i][j];
}
template <typename T>
void Table<T>::set_table(std::vector<std::vector<T>> &table_) {
    table = table_;
}
template <typename T> std::vector<std::vector<T>> Table<T>::get_table() {
  return table;
}
} // namespace table