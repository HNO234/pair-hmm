#pragma once
#include "table/STR_table.hpp"
#include <vector>
namespace pairhmm {

enum PairHMMAlgorithm { NAIVE, SUZUKI_KASAHARA };

extern std::vector<std::vector<double>> standard_gop_vector;

table::STRTable<double> transform_into_STRTable(std::vector<std::vector<double>> table);

extern table::STRTable<double> standard_gop;

extern std::vector<std::vector<double>> standard_gcp_vector;

extern table::STRTable<double> standard_gcp;

extern std::vector<std::vector<double>> base_match_prob;

} // namespace pairhmm