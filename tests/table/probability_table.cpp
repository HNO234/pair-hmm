#include "table/probability_table.hpp"
#include <catch_amalgamated.hpp>
#include <cmath>
#include <vector>

TEST_CASE("table") {
  auto is_close = [](double a, double b) {
    return std::abs(a - b) < 1e-8;
  };
  table::ProbabilityTable<double> table(2, 2);
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

TEST_CASE("shift") {
  auto is_close = [](double a, double b) {
    return std::abs(a - b) < 1e-8;
  };
  std::vector<std::vector<double>> new_table = {{1.2, 2.3, 3.4},
                                                {4.5, 5.6, 6.7}};
  table::ProbabilityTable<double> table(2, 2);
  SECTION("shifted_down") {
    table.set_table(new_table);
    auto table_shifted_down = table.shifted_down();
    REQUIRE(is_close(table_shifted_down.get_cell(1, 1), 2.3));
    REQUIRE(is_close(table_shifted_down.get_cell(0, 2), 0));
  }
  SECTION("shifted_right") {
    table.set_table(new_table);
    auto table_shifted_down = table.shifted_right();
    REQUIRE(is_close(table_shifted_down.get_cell(1, 1), 4.5));
    REQUIRE(is_close(table_shifted_down.get_cell(0, 0), 0));
  }
}

TEST_CASE("+/- operator") {
  auto is_close = [](double a, double b) {
    return std::abs(a - b) < 1e-8;
  };
  std::vector<std::vector<double>> new_table_1 = {{1.2, 2.3, 3.4},
                                                  {4.5, 5.6, 6.7}};
  std::vector<std::vector<double>> new_table_2 = {{1, 2, 3}, {4, 5, 6}};
  SECTION("+ operator") {
    table::ProbabilityTable<double> table1(2, 2);
    table::ProbabilityTable<double> table2(2, 2);
    table1.set_table(new_table_1);
    table2.set_table(new_table_2);
    auto table3 = table1 + table2;
    REQUIRE(is_close(table3.get_cell(1, 1), 10.6));
    REQUIRE(is_close(table3.get_cell(0, 0), 2.2));
  }
  SECTION("- operator") {
    table::ProbabilityTable<double> table1(2, 2);
    table::ProbabilityTable<double> table2(2, 2);
    table1.set_table(new_table_1);
    table2.set_table(new_table_2);
    auto table3 = table1 - table2;
    REQUIRE(is_close(table3.get_cell(1, 1), 0.6));
    REQUIRE(is_close(table3.get_cell(0, 0), 0.2));
  }
}

TEST_CASE("is_close function") {
  auto is_close = [](double a, double b) {
    return std::abs(a - b) < 1e-8;
  };
  std::vector<std::vector<double>> new_table_1 = {{1.2, 2.3, 3.4},
                                                  {4.5, 5.6, 6.7}};
  std::vector<std::vector<double>> new_table_2 = {{1, 2, 3}, {4, 5, 6}};
  SECTION("is_close") {
    table::ProbabilityTable<double> table1(2, 2);
    table::ProbabilityTable<double> table2(2, 2);
    table1.set_table(new_table_1);
    table2.set_table(new_table_2);
    REQUIRE(table::is_close(table1, table2, double{1}));
    table1.set_cell(1, 2, 100);
    REQUIRE(!table::is_close(table1, table2, double{1}));
    std::vector<std::vector<double>> new_table_3 = {{1.2, 2.3, 3.4, 0},
                                                    {4.5, 5.6, 6.7, 0}};
    table1.set_table(new_table_3);
    REQUIRE(!table::is_close(table1, table2, double{1}));
  }
}