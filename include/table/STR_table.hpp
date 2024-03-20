#pragma once

#include "table/table.hpp"

namespace table {
template <typename T> class STRTable : public Table<T> {
public:
  using Table<T>::Table;
  
};

template class STRTable<double>;
} // namespace table