#pragma once
#include "data.h"
#include <cstdint>
#include <format>
#include <iostream>
#include <random>

template <int max_length> class BnbOptimizer
{
public:
    BnbOptimizer(encodingData<max_length> &ed, int numKeys, int pruneBound, bool print)
        : encoding(ed), numKeys(numKeys), print(print) {
        minDupCnt = pruneBound;
        bestMapping.resize(ed.mapping.numRadical + 1, 0);
        pruneStats.resize(ed.mapping.numRadical + 2, 0);
    }
    bool solve(const std::vector<int> &useRadical) {
        startTime = std::chrono::steady_clock::now();
        bool success = search(0, useRadical);
        if (print)
            std::cout << "\nOptimization Complete. Min Duplicates: " << minDupCnt << std::endl;
        return success;
    }
    std::vector<uint16_t> bestMapping;
    int minDupCnt;

private:
    void printProgress() {
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();
        std::cout << std::format("\n[Time: {}s | Nodes: {} | Leaves: {} | Best: {}]\n Pruned at: ", elapsed, numSearched,
                                 numLeaf, minDupCnt);
        for (int i = 1; i < pruneStats.size(); ++i) {
            if (pruneStats[i] > 0) {
                std::cout << std::format("{}:{} \t", i, pruneStats[i]);
            }
        }
        std::cout << std::flush;
    }

    bool search(const size_t pos, const std::vector<int> &useRadical) {
        numSearched++;
        if (print && numSearched % 1000000 == 0) {
            printProgress();
        }

        // 终止
        if (pos >= useRadical.size()) {
            if (encoding.dupCnt < minDupCnt) {
                numLeaf++;
                minDupCnt = encoding.dupCnt;
                bestMapping = encoding.mapping.radicalToKey;
                return true;
            }
            return false;
        }

        // 剪枝
        if (encoding.dupCnt >= minDupCnt) {
            pruneStats[pos + 1]++;
            return false;
        }

        const int &radIdx = useRadical[pos];
        bool success = false;
        // 分支
        const auto &affectedHanzi = encoding.mapping.radicalToHanzi[radIdx];
        std::vector<std::array<uint16_t, max_length>> savedEncodings;
        savedEncodings.reserve(affectedHanzi.size());
        for (int hanziIdx : affectedHanzi)
            savedEncodings.push_back(encoding.keyEncoding[hanziIdx]);

        for (uint16_t k = 1; k <= numKeys; ++k) {
            encoding.modifyMapping(radIdx, k);
            success = success || search(pos + 1, useRadical);
            encoding.modifyMapping(radIdx, 0);
        }
        return success;
    }
    encodingData<max_length> &encoding;
    int numKeys;
    std::mt19937 gen;
    bool print;
    long long numSearched = 0;
    long long numLeaf = 0;
    std::vector<long long> pruneStats;
    std::chrono::steady_clock::time_point startTime;
};