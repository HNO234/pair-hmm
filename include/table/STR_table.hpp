#pragma once

#include "table/table.hpp"

namespace table {
template <typename T> class STRTable : public Table<T> {
public:
  using Table<T>::Table;
  
};
} // namespace table