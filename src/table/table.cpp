#include "table/table.hpp"
#include <vector>

namespace table {
template <typename T> Table<T>::Table() : rows(0), columns(0), table() {}
template <typename T>
Table<T>::Table(size_t rows_, size_t columns_)
    : rows(rows_), columns(columns_), table(rows_, std::vector<T>(columns_)) {}
template <typename T> void Table<T>::set_cell(size_t i, size_t j, T element) {
  table[i][j] = element;
}
template <typename T> T Table<T>::get_cell(size_t i, size_t j) const {
  return table[i][j];
}
template <typename T>
void Table<T>::set_table(std::vector<std::vector<T>> &table_) {
    table = table_;
    rows = table.size();
    columns = table[0].size();
}
template <typename T> std::vector<std::vector<T>> Table<T>::get_table() {
  return table;
}

template <typename T>
size_t Table<T>::get_rows() { return this->rows; }

template <typename T>
size_t Table<T>::get_columns() { return this->columns; }
} // namespace table