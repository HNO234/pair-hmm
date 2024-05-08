#pragma once
#include "utils/constant.hpp"
#include "utils/options.hpp"
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/file_io/bam.hpp>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <variant>
#include <tuple>

namespace pairhmm {
    auto get_haplotype(std::string haplotype_filename) {
        auto haplotype_fa =
        std::ifstream{haplotype_filename};
        auto haplotype_refs = 
            std::ranges::istream_view<biovoltron::FastaRecord<true>>(haplotype_fa);
        biovoltron::istring haplotype;
        for (auto &ref : haplotype_refs)
            haplotype += ref.seq;
        std::ranges::transform(haplotype, haplotype.begin(),
            [](auto &c) { return c % 4; });
        return haplotype;
    }

    auto get_reads(std::string read_bam) {
        int cnt = 0;
        biovoltron::IBamStream fin(read_bam);
        biovoltron::SamHeader header;
        std::vector<biovoltron::SamRecord<false>> records;
        biovoltron::SamRecord<false> record;
        fin >> header;
        while (fin >> record) {
            records.emplace_back(record);
            cnt++;
            if (cnt > 10)
                break;
        }
        return records;
    }

    auto get_gop_gcp(std::string tables_file_name) {
        std::vector<std::vector<long double>> 
            gop_vector(8, std::vector<long double>(20)), 
            gcp_vector(8, std::vector<long double>(20));
        std::ifstream tables_file;
        tables_file.open(tables_file_name);
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 20; j++)
            tables_file >> gop_vector[i][j];
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 20; j++)
            tables_file >> gcp_vector[i][j];
        auto gop = table::STRTable<long double>();
        gop.set_table(gop_vector);
        auto gcp = table::STRTable<long double>();
        gcp.set_table(gcp_vector);
        tables_file.close();
        return std::tuple(gop, gcp);
    }

    // 2 13 23 37
    auto score_to_phred(std::string quality) {
        std::vector<int> res;
        for (char c : quality) {
            res.emplace_back(c - 33);
        }
        return res;
    }
} // pairhmm