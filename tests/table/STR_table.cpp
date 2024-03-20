#include "table/STR_table.hpp"
#include <catch_amalgamated.hpp>
#include <cmath>
#include <vector>

TEST_CASE("table") {
  table::STRTable<double> table(2, 2);
  auto is_close = [](double a, double b) {
    return std::abs(a - b) < 1e-8;
  };
  SECTION("set/get cell") {
    table.set_cell(0, 0, 1.2);
    auto val = table.get_cell(0, 0);
    REQUIRE(is_close(val, 1.2));

    table.set_cell(1, 1, -0.9);
    val = table.get_cell(1, 1);
    REQUIRE(is_close(val, -0.9));

    val = table.get_cell(0, 1);
    REQUIRE(is_close(val, 0));
  }
  SECTION("set/get table") {
    std::vector<std::vector<double>> new_table = {{1.2, 2.3, 3.4},
                                                  {4.5, 5.6, 6.7}};
    table.set_table(new_table);
    auto new_table_ret = table.get_table();
    REQUIRE(new_table_ret == new_table);
  }
}