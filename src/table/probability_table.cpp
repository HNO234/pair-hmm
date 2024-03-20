#include "table/probability_table.hpp"

namespace table {
template <typename T> ProbabilityTable<T> ProbabilityTable<T>::shifted_down() {
  return (*this);
}
template <typename T> ProbabilityTable<T> ProbabilityTable<T>::shifted_right() {
  return (*this);
}
template <typename T>
ProbabilityTable<T>
ProbabilityTable<T>::operator-(const ProbabilityTable &right_table) const {
  return (*this);
}
template <typename T>
ProbabilityTable<T>
ProbabilityTable<T>::operator+(const ProbabilityTable &right_table) const {
  return (*this);
}

template <typename T>
bool is_close(ProbabilityTable<T> &table1, ProbabilityTable<T> &table2,
              T tolerance) {
  return false;
}

template <>
bool is_close(ProbabilityTable<double> &table1, ProbabilityTable<double> &table2,
              double tolerance) {
  return false;
}
} // namespace table