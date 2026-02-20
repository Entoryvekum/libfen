#pragma once
#include "data.h"
#include <vector>
#include <random>
#include <iostream>
#include <format>
#include <chrono>
#include <limits>

template <int max_length> class TabuOptimizer
{
public:
    TabuOptimizer(encodingData<max_length> &ed, int numKeys, int preAllocRadical, bool print)
        : encoding(ed), numKeys(numKeys), preAllocRadical(preAllocRadical), 
          gen(0 /*std::random_device{}()*/), print(print) {
        // tabuList[radicalIdx][keyIdx] 存储该移动解禁的迭代轮次
        tabuList.resize(encoding.numRadical + 1, std::vector<int>(numKeys + 1, 0));
    }
    
    void solve(int maxIterations, int tabuTenure, int neighborSize, int bestBound) {
        bestMapping = encoding.mapping.radicalToKey;
        int currentEnergy = encoding.dupCnt;
        int bestEnergy = currentEnergy;

        // 随机生成器初始化
        std::uniform_int_distribution<int> radicalGen(preAllocRadical + 1, encoding.numRadical);
        std::uniform_int_distribution<int> keyGen(1, numKeys);

        auto startTime = std::chrono::steady_clock::now();

        for (int iter = 1; iter <= maxIterations; ++iter) {
            if (currentEnergy <= bestBound)
                break;

            struct Move {
                int radIdx;
                uint16_t oldKey;
                uint16_t newKey;
                int delta;
            };
            
            // 随机采样 neighborSize个移动，从中选出最好的
            Move bestLocalMove = {-1, 0, 0, std::numeric_limits<int>::max()};
            bool foundMove = false;

            for (int k = 0; k < neighborSize; ++k) {
                // 1. 生成随机移动
                int radIdx = radicalGen(gen);
                uint16_t oldKey = encoding.mapping.radicalToKey[radIdx], newKey = keyGen(gen);
                while (oldKey == newKey) {
                    newKey = keyGen(gen);
                }

                // 2. 试探性修改并计算 Delta
                encoding.modifyMapping(radIdx, newKey);
                int newTotal = encoding.dupCnt;
                int delta = newTotal - currentEnergy;
                
                // 3. 判断是否由禁忌或渴望准则接受
                if (iter > tabuList[radIdx][newKey] || newTotal < bestEnergy) {
                    if (delta < bestLocalMove.delta) {
                        bestLocalMove = {radIdx, oldKey, newKey, delta};
                        foundMove = true;
                    }
                }

                // 4. 回溯（撤销修改），准备测试下一个邻居
                encoding.modifyMapping(radIdx, oldKey);
            }

            if (foundMove) {
                // 应用最佳修改
                encoding.modifyMapping(bestLocalMove.radIdx, bestLocalMove.newKey);
                currentEnergy = encoding.dupCnt;

                // 更新全局最优
                if (currentEnergy < bestEnergy) {
                    bestEnergy = currentEnergy;
                    bestMapping = encoding.mapping.radicalToKey;
                }

                // 更新禁忌表
                tabuList[bestLocalMove.radIdx][bestLocalMove.oldKey] = iter + tabuTenure;
            }
            if (print && iter % 100 == 0) {
                auto now = std::chrono::steady_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();
                std::cout << std::format("\r[Iter: {} | Best: {} | Current: {} | Time: {}s]", 
                    iter, bestEnergy, currentEnergy, elapsed) << std::flush;
            }
        }
        if (print)
            std::cout << "\nTS | Optimization Complete. Min Duplicates: " << bestEnergy << std::endl;
        encoding.mapping.radicalToKey = bestMapping;
        encoding.createEncodingFromMapping();
        encoding.buildHash();
    }

private:
    encodingData<max_length> &encoding;
    int numKeys;
    int preAllocRadical;
    std::vector<uint16_t> bestMapping;
    std::vector<std::vector<int>> tabuList; 
    std::mt19937 gen;
    bool print;
};