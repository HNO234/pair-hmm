#pragma once

#include "table/STR_table.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <biovoltron/utility/istring.hpp>
#include <tuple>
#include "pairhmm/pairhmm.hpp"
#include "table/STR_table.hpp"
#include <biovoltron/file_io/cigar.hpp>
#include <biovoltron/utility/istring.hpp>
#include <iostream>

namespace pairhmm {
template <typename T> class PairHMM {
protected:
  biovoltron::istring haplotype, read;
  table::STRTable<T> gop, gcp;
  std::vector<size_t> calculate_1() {
  std::vector<size_t> output(read.size());
  int seqLength = read.size();
  auto& bases = read;
  const int rightMargin = seqLength - 1;
  auto last = bases[rightMargin];

  // backward phase:
  int carryBack = output[rightMargin] = 1;
  for (int position = rightMargin - 1; position >= 0; position--) {
      const auto next = bases[position];
      output[position] = next == last ? ++carryBack : (carryBack = 1);
      last = next;
  }
  // forward phase:
  // last = sequence[0]; // already true.
  int prevRunLength = output[0];
  for (int position = 1; position <= rightMargin; position++) {
      const auto next = bases[position];
      if (next == last) {
          output[position] = prevRunLength;
      } else {
          const int thisRunLength = output[position];
          if (prevRunLength < thisRunLength) { // so we propagate this position run-length to the prev position if longer.
              output[position - 1] = thisRunLength;
          }
          last = next;
          prevRunLength = thisRunLength;
      }
  }
  return output;
}
  std::vector<size_t> calculate_2_or_above(size_t period) {
  std::vector<size_t> output(read.size());
  if (read.size() < period)
    return output;
  
  std::vector<size_t> runLengthBuffer(read.size() + 1);
  std::vector<size_t> heap(9);
  auto& bases = read;

  auto fixHeapUp = [&](std::vector<size_t> &heap, size_t idx) {
    const int value = heap[idx];
    do {
        const int upIdx = (idx - 1) >> 1;
        const int upValue = heap[upIdx];
        if (upValue >= value) {
            break;
        } else {
            heap[idx] = upValue;
            idx = upIdx;
        }
    } while (idx > 0);
    heap[idx] = value;
  };

  auto fixHeapDown = [&](std::vector<size_t> &heap, int idx, const int heapSize) {
      const int value = heap[idx];
      while (true) {
          const int idxRight = (idx + 1) << 1;
          const int rightValue = idxRight < heapSize ? heap[idxRight] : -1;
          const int idxLeft = idxRight - 1;
          const int leftValue = idxLeft < heapSize ? heap[idxLeft] : -1;
          if (rightValue > value) {
              if (rightValue > leftValue) {
                  heap[idx] = rightValue;
                  idx = idxRight;
              } else {
                  heap[idx] = leftValue;
                  idx = idxLeft;
              }
          } else if (leftValue > value) {
              heap[idx] = leftValue;
              idx = idxLeft;
          } else {
              break;
          }
      }
      heap[idx] = value;
  };

  auto fixHeap = [&](std::vector<size_t> &heap, size_t idx, size_t heapSize) {
     if (idx == 0 || heap[(idx - 1) >> 1] > heap[idx]) {
        fixHeapDown(heap, idx, heapSize);
    } else {
        fixHeapUp(heap, idx);
    }
  };

  int position, matchedCycles, cycleIndex, positionPlusPeriod;
  int seqLength = read.size();

  // suffixes of zero run-lengths at the end of the sequence (period - 1).
  for (position = seqLength - 1, cycleIndex = period; cycleIndex > 1; position--, cycleIndex--) {
      runLengthBuffer[position] = 0;
  }

  // backward phase.
  // at the end of this phase runLength[x] will have the total number of equal repeats downstream with length starting at position X
  // inclusive.
  int carryBack = runLengthBuffer[position--] = 1; //prevValue holds the num of repeats reported in the previous (+1) position.
  for (positionPlusPeriod = position + period, matchedCycles = 0; position >= 0; position--, positionPlusPeriod--) {
      if (bases[position] == bases[positionPlusPeriod]) { // we keep matching repeat unit bases.
          if (++matchedCycles == period) { // we got a new full repeat matched so the run length increases:
              runLengthBuffer[position] = ++carryBack;
              matchedCycles = 0; // we reset the match-run-length to 0 as we start over (a new repeat).
          } else { // if we haven't completed a repeat we simply copy the run length from the +1 position base.
              runLengthBuffer[position] = carryBack;
          }
      } else { // we bump into a mismatch that ends the run.
          carryBack = runLengthBuffer[position] = 1;
          matchedCycles = 0; // we reset the match-run-length to 0.
      }
  }

  // Now we propagate forward the number of equivalent repeats up stream:
  // So at the end runLength[X] == the repeats (up and down-stream) of the
  // unit of length period that starts at X.
  // The outside for-loop iterates over different repeat unit start offsets (here named cycles):
  for (cycleIndex = 0; cycleIndex < period; cycleIndex++) {
      // The left most repeated unit runLength[i] contains the actual run length for all the matching
      // units to the right. We copy that value forward to the other run-length units.
      // We do this by iterating over consecutive repeat runs.
      for (position = cycleIndex; position < seqLength; position += period) {
          const int totalRunLength = runLengthBuffer[position];
          for (int repeatInRun = 1; repeatInRun < totalRunLength; repeatInRun++) {
              runLengthBuffer[position += period] = totalRunLength;
          }
      }
  }

  // NOTE: this line below is added solely to replicate an oddity of DRAGEN's algorithm that discounts one repeat
  // only at the beginning of the sequence for period 2 or above (e.g. ^CACATG would yield periods 122211 repeats
  // 12221 when the natural/intuitive solution is period 222211 reps 222211). This is true for the original Matlab
  // and latest DRAGEN so we must add this here:
  if (runLengthBuffer[0] > 1) runLengthBuffer[0]--;
  // Finally we need to propagate the best run-lengths to neighbor positions
  // so that a position has the maximum run-length of all the possible
  // units that contain the position + the unit that follow up-stream.
  const int heapSize = period + 1;
  for (int i = 0; i < heapSize; i++)
    heap[i] = 0;
  int currentMax = heap[0] = runLengthBuffer[0];
  // the first period length prefix is trivial:
  // a position's max run-length is the largest seen so far and the following position's
  // We use this opportunity to populate the
  // run-length max-heap
  const int stop0 = std::min((int)period, seqLength - 1);
  for (position = 0; position < stop0; ) {
      output[position] = currentMax = std::max(currentMax, (int)(heap[position + 1] = runLengthBuffer[position + 1]));
      position++;
      fixHeap(heap, position, heapSize);
  }

  runLengthBuffer[seqLength] = 1; // we add this 1 as the runLength after the last sequence position which in fact does not exist.
                                  // this allows us to save a conditional to avoid a index-out-of-range.
  for (int outPosition = 0; position < seqLength;) {
      const int valueOut = runLengthBuffer[outPosition++]; // value leaving the heap.
      const int valueIn = runLengthBuffer[++position]; // value entering the heap.
      if (valueIn != valueOut) { // these two are the same (often the case) we don't need to do a heap updating at all.
          const int fixHeapIdx = std::find(heap.begin(), heap.end(), valueOut) - heap.begin(); // O(n) search although usually n is very small.
          heap[fixHeapIdx] = valueIn;
          fixHeap(heap, fixHeapIdx, heapSize);
          currentMax = heap[0];
      }
      output[position - 1] = currentMax; // we use the heap's max as the final run-length for the prev position.
  }

  return output;
}

// int calculateExpectedRepeatLength(biovoltron::istring& sequence, const int position, const int period) {
//         if (period > (int)sequence.size()) {
//             return 0;
//         }
//         // The simple story:
//         // we must calculate the maxium number of adjacent repeats for any of the overlapping k-mers (k = period)
//         // at that position + the k-mer that immediately follows it
//         // Details:
//         //     The outer loop goes thru the start of each of those k-mers. start = [position - period ... position + 1]
//         //        Then two inner loops simply calculate the number of copies downstream (forward) and upstream (backwards)
//         //     There is an adjustment a the end due to a quirkiness of DRAGEN's STR analyzer in where a copy is discounted
//         //        in the first k-1 position of the sequence.
//         // int start = Math.max(0, position - period + 1);
//         int start = std::max(0, position - period + 1);
//         int max = 0;
//         for (int end = start + period; start <= position + 1 && end <= (int)sequence.size(); start++, end++) {
//             biovoltron::istring unit(end-start,0);
//             std::copy(sequence.begin() + start, sequence.begin() + end, unit.begin());
//             // final byte[] unit = Arrays.copyOfRange(sequence, start, end);
//             int forward = 1; // = 1 since we must count the unit starting at position.
//             for (int offset = end; offset <= (int)sequence.size() - period; offset += period) {
//                 biovoltron::istring other(period,0);
//                 std::copy(sequence.begin() + offset, sequence.begin() + offset + period, other.begin());
//                 std::cerr << position << ' ' << offset << ' ' << offset + period << ' ' << unit << ' ' << other << std::endl;
//                 // final byte[] other = Arrays.copyOfRange(sequence, offset, offset + period);
//                 if (unit != other) {
//                     break;
//                 }
//                 forward++;
//             }
//             int backward = 0;
//             for (int offset = start - period; offset >= 0; offset -= period) {
//                 biovoltron::istring other(period,0);
//                 std::copy(sequence.begin() + offset, sequence.begin() + offset + period, other.begin());
//                 // final byte[] other = Arrays.copyOfRange(sequence, offset, offset + period);
//                 if (unit != other) {
//                     break;
//                 }
//                 backward++;
//             }
//             int candidate = forward + backward;
//             // Matching strange behavior in DRAGEN where repeat lengths at first position for period larger than one
//             // are reduced by one (minimum one. (but no less than 1.
//             if (position < period - 1 && start == 0) {
//                 // candidate = Math.max(1, candidate - 1);
//                 candidate = std::max(1, candidate - 1);
//             }

//             // max = Math.max(max, candidate);
//             max = std::max(max, candidate);
//         }
//         return max;
//     }

public:
  PairHMM() : haplotype(), read(), gop(), gcp() {}
  PairHMM(biovoltron::istring haplotype_, biovoltron::istring read_,
                    table::STRTable<T> gop_, table::STRTable<T> gcp_)
    : haplotype(haplotype_), read(read_), gop(gop_), gcp(gcp_) {}
// vector of period, vector of repeat size
  std::tuple<std::vector<size_t>, std::vector<size_t>> get_read_best_repeat() {
  std::vector<size_t> best_period(read.size() + 1), best_size(read.size() + 1);

  std::vector<std::vector<size_t>> repeats_by_period_and_position(8);

  repeats_by_period_and_position[0] = calculate_1();
  for (size_t period = 2; period <= 8; period++)
    repeats_by_period_and_position[period - 1] = calculate_2_or_above(period);

  for (size_t pos = size_t{}; pos < read.size(); pos++)
    best_period[pos + 1] = 1, best_size[pos + 1] = repeats_by_period_and_position[0][pos];

  for (size_t period = 2; period <= 8; period++) {
    for (size_t pos = size_t{}; pos < read.size(); pos++) {
      if (repeats_by_period_and_position[period - 1][pos] > best_size[pos + 1]) {
        best_period[pos + 1] = period;
        best_size[pos + 1] = repeats_by_period_and_position[period - 1][pos];
      }
    }
  }
  best_period[0] = 1, best_size[0] = 1;
  best_period.back() = 1;
  best_size.back() = 1;
  return std::make_tuple(best_period, best_size);
}
  virtual biovoltron::Cigar get_cigar() = 0;
  virtual void run_alignment() = 0;
};

template class PairHMM<double>;
} // namespace pairhmm