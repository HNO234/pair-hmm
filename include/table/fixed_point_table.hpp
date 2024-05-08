#pragma once

#include "table/table.hpp"
#include "utils/fixed_point.hpp"

namespace table {
template <typename T> class FixedPointTable : public Table<T> {
public:
  using Table<T>::Table;
};
} // namespace table