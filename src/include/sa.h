#pragma once
// #include "data.h"
#include "data_numerical_hashonly.h"
#include <cmath>
#include <iostream>
#include <random>

template <int maxLength, int numKeys> class SAOptimizer
{
public:
    struct SAParameters
    {
        double T_start;
        int bestBound = 0;
        int pos = -1;
    };
    SAOptimizer(mapping<maxLength, numKeys> &ed, int preAllocRadical, bool print)
        : encoding(ed), preAllocRadical(preAllocRadical), gen(std::random_device{}()), print(print) {}

    template <typename nextTempType> int solve(SAParameters param, nextTempType nextTemp) {
        if (print) {
            if (param.pos != -1)
                std::cout << std::format("\nSA | Optimization Start. Duplicates: {}, Pos: {}", encoding.dupCnt, param.pos)
                          << std::flush;
            else
                std::cout << std::format("\nSA | Optimization Start. Duplicates: {}", encoding.dupCnt) << std::flush;
        }
        double T = param.T_start;
        bestMapping = encoding.radicalToKey;
        int currentEnergy = encoding.dupCnt;
        int bestEnergy = currentEnergy;
        int numItr = 0;

        std::uniform_int_distribution<int> radicalGenNoPreAlloc(preAllocRadical + 1, encoding.numRadical);
        std::uniform_int_distribution<int> radicalGenWithPreAlloc(1, encoding.numRadical);
        std::uniform_int_distribution<int> radicalGenOnlyPreAlloc(1, preAllocRadical);
        std::uniform_int_distribution<int> keyGen(1, numKeys );
        std::uniform_real_distribution<double> probGen(0.0, 1.0);

        auto startTime = std::chrono::steady_clock::now();

        int test = 0;
        auto pointwiseModification = [&]() {
            // 1. 随机选择一个字根和新按键
            int radIdx = radicalGenNoPreAlloc(gen);
            uint16_t oldKey = encoding.radicalToKey[radIdx], newKey = keyGen(gen);
            while (oldKey == newKey)
                newKey = keyGen(gen);

            // 2. 尝试修改
            encoding.modifyRadical(radIdx, newKey);
            int deltaE = encoding.dupCnt - currentEnergy;
            // 3. Metropolis准则
            if (deltaE < 0 || probGen(gen) < std::exp(-deltaE / T)) {
                // 接受新解
                currentEnergy = encoding.dupCnt;
                if (currentEnergy < bestEnergy) {
                    bestEnergy = currentEnergy;
                    bestMapping = encoding.radicalToKey;
                }
            }
            else {
                // 拒绝新解
                encoding.modifyRadical(radIdx, oldKey);
            }
        };

        auto exchangeModification = [&]() {
            // 1. 随机选择两个字根
            uint16_t radIdx1, radIdx2;
            auto getRandomPair = [&]() {
                radIdx1 = radicalGenWithPreAlloc(gen);
                if (radIdx1 > preAllocRadical) {
                    radIdx2 = radicalGenNoPreAlloc(gen);
                    while (radIdx2 == radIdx1)
                        radIdx2 = radicalGenNoPreAlloc(gen);
                }
                else {
                    radIdx2 = radicalGenOnlyPreAlloc(gen);
                    while (radIdx2 == radIdx1)
                        radIdx2 = radicalGenOnlyPreAlloc(gen);
                }
            };
            do {
                getRandomPair();
            } while (encoding.radicalToKey[radIdx1] == encoding.radicalToKey[radIdx2]);

            // 2. 交换两个字根的键位
            auto exchangeKey = [&](uint16_t a, uint16_t b) {
                auto tmp = encoding.radicalToKey[b];
                encoding.modifyRadical(b, encoding.radicalToKey[a]);
                encoding.modifyRadical(a, tmp);
            };
            exchangeKey(radIdx1, radIdx2);

            // 3. Metropolis准则
            int deltaE = encoding.dupCnt - currentEnergy;
            if (deltaE < 0 || probGen(gen) < std::exp(-deltaE / T)) {
                currentEnergy = encoding.dupCnt;
                if (currentEnergy < bestEnergy) {
                    bestEnergy = currentEnergy;
                    bestMapping = encoding.radicalToKey;
                }
            }
            else {
                exchangeKey(radIdx1, radIdx2);
            }
        };

        while (T > 0) {
            ++numItr;
            double threshold;
            if (T <= 0.35)
                threshold = 0.8;
            else if (T <= 0.85)
                threshold = 0.8;
            else if (T <= 1.5)
                threshold = 0.8;
            else
                threshold = 0.75;
            if (probGen(gen) <= threshold)
                pointwiseModification();
            else
                exchangeModification();
            if (encoding.dupCnt <= param.bestBound)
                break;
            T = nextTemp(T, numItr, encoding.dupCnt);
            if (encoding.dupCnt <= param.bestBound)
                break;
            if (print && numItr % 10000 == 0) {
                auto now = std::chrono::steady_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
                std::cout << std::format("\r[T: {:.4f} | Best: {} | Current: {} | Time: {:.3f}s | Steps: {}]", T, bestEnergy,
                                         currentEnergy, elapsed/1000.0, numItr)
                          << std::flush;
            }
        }
        if (print) {
            if (param.pos != -1)
                std::cout << std::format("\nSA | Optimization Complete. Min Duplicates: {}, Pos: {}", bestEnergy, param.pos)
                          << std::flush;
            else
                std::cout << std::format("\nSA | Optimization Complete. Min Duplicates: {}", bestEnergy) << std::flush;
        }
        return bestEnergy;
    }
    void applyBest() {
        encoding.radicalToKey = bestMapping;
        encoding.createEncodingFromMapping();
        encoding.buildHash();
    }

private:
    mapping<maxLength, numKeys> &encoding;
    int preAllocRadical;
    std::vector<uint16_t> bestMapping;
    std::mt19937 gen;
    bool print;
};
