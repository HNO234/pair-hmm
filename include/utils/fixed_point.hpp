#pragma once

#include "utils/constant.hpp"
#include <cmath>
#include <iostream>

namespace pairhmm {

template <size_t PRECISION> struct FixedPoint {
    int64_t integer, fractional;

    FixedPoint(): integer(0), fractional(0) {}
    FixedPoint(int x, int y): integer(x), fractional(y) {}
    FixedPoint(int64_t x, int64_t y): integer(x), fractional(y) {}
    FixedPoint(int x): integer(x), fractional(0) {}
    FixedPoint(int64_t x): integer(x), fractional(0) {}
    FixedPoint(double x) {
        integer = std::floor(x);
        x -= std::floor(x);

        fractional = 0;
        for (int i = 1; i <= PRECISION; i++) {
            x *= 2;
            if (std::floor(x) > 0) {
                fractional |= (1ULL << (FLOAT_START_INDEX - i));
            }
            x -= std::floor(x);
        }
    }
    FixedPoint(long double x) {
        integer = std::floor(x);
        x -= std::floor(x);

        fractional = 0;
        for (int i = 1; i <= PRECISION; i++) {
            x *= 2;
            if (std::floor(x) > 0) {
                fractional |= (1ULL << (FLOAT_START_INDEX - i));
            }
            x -= std::floor(x);
        }
    }
template <size_t PRECISION2>
    FixedPoint(FixedPoint<PRECISION2> x) {
        FixedPoint<PRECISION> ret;

        if (PRECISION >= PRECISION2) {
            ret.integer = x.integer;
            ret.fractional = x.fractional;
        } else {
            ret.integer = x.integer;
            auto mask = ((1ULL << FLOAT_START_INDEX) - 1) ^ 
                        ((1ULL << (FLOAT_START_INDEX - PRECISION)) - 1);
            ret.fractional = x.fractional & (mask);
        }
    }

template <size_t PRECISION2>
    auto operator=(FixedPoint<PRECISION2> x) {
        FixedPoint<PRECISION> ret;

        if (PRECISION >= PRECISION2) {
            ret.integer = x.integer;
            ret.fractional = x.fractional;
        } else {
            ret.integer = x.integer;
            auto mask = ((1ULL << FLOAT_START_INDEX) - 1) ^ 
                        ((1ULL << (FLOAT_START_INDEX - PRECISION)) - 1);
            ret.fractional = x.fractional & (mask);
        }
        return ret;
    }

    long double to_float() {
        double ret = integer + (long double)fractional / (1ULL << FLOAT_START_INDEX);
        return ret;
    }
};

template<size_t T>
std::ostream &operator<<(std::ostream &out, FixedPoint<T> x) {
  out << x.to_float();
  return out;
}

template<size_t T1, size_t T2 >
auto operator+(
    FixedPoint<T1> const left,
    FixedPoint<T2> const right
) -> FixedPoint<std::max(T1, T2)> 
{
    auto ret = FixedPoint<std::max(T1, T2)>(
        left.integer + right.integer,
        left.fractional + right.fractional
    );

    if (ret.fractional >= (1ULL << FLOAT_START_INDEX)) {
        ret.integer++;
        ret.fractional -= (1ULL << FLOAT_START_INDEX);
    }

    return ret;
}

template<size_t T1, size_t T2 >
auto operator-(
    FixedPoint<T1> const left,
    FixedPoint<T2> const right
) -> FixedPoint<std::max(T1, T2)> 
{
    auto ret = FixedPoint<std::max(T1, T2)>(
        left.integer - right.integer,
        left.fractional - right.fractional
    );

    if (ret.fractional < 0) {
        ret.integer--;
        ret.fractional += (1ULL << FLOAT_START_INDEX);
    }

    return ret;
}

template<size_t T>
auto operator-(
    FixedPoint<T> const num
) -> FixedPoint<T> 
{
    auto ret = FixedPoint<T>(0.0) - num;
    return ret;
}

template<size_t T1, size_t T2 >
auto min(
    FixedPoint<T1> const left,
    FixedPoint<T2> const right
) -> FixedPoint<std::max(T1, T2)> 
{
    if (left.integer < right.integer) {
        return left;
    } else if (left.integer > right.integer) {
        return right;
    } else {
        if (left.fractional < right.fractional) {
            return left;
        } else {
            return right;
        }
    }
}

template class FixedPoint<16>;

} // namespace pairhmm