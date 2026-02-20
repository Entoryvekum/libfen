#pragma once
#include "BS_thread_pool.hpp"
#include "data.h"
#include "sa.h"
#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

template <int max_length> class MemeticOptimizer
{
public:
    struct Individual
    {
        std::vector<uint16_t> mapping; // radicalToKey
        int fitness;                   // dupCnt
    };

    struct optimizationPlan
    {
        optimizationPlan(const mappingData<max_length> &mp, int numKeys, int preAllocRadical)
            : mapping(mp), encoding(mapping), sa(encoding, numKeys, preAllocRadical, false) {}
        mappingData<max_length> mapping;
        encodingData<max_length> encoding;
        SAOptimizer<max_length> sa;
    };

    MemeticOptimizer(encodingData<max_length> &ed, int numKeys, int preAllocRadical, int popSize)
        : encoding(ed), numKeys(numKeys), preAllocRadical(preAllocRadical), popSize(popSize), population(popSize),
          gen(std::random_device{}()) {
        plans.reserve(popSize);
        for (int i = 0; i < popSize; ++i) {
            plans.push_back(std::make_unique<optimizationPlan>(encoding.mapping, numKeys, preAllocRadical));
        }
    }

    // 初始化种群
    void initialize(double mutationRate = -1) {
        std::uniform_int_distribution<int> keyDist(1, numKeys);
        for (auto &ind : population) {
            ind.mapping = encoding.mapping.radicalToKey; // 复制初始状态
            std::uniform_real_distribution<double> probGen(0.0, 1.0);
            for (int i = preAllocRadical + 1; i <= encoding.numRadical; ++i) {
                if (probGen(gen) < mutationRate)
                    ind.mapping[i] = keyDist(gen);
            }
            evaluate(ind, *plans[0]);
            ind.fitness = plans[0]->encoding.dupCnt;
        }
        std::sort(population.begin(), population.end(),
                  [](const Individual &a, const Individual &b) { return a.fitness < b.fitness; });
        bestInd = population[0];
    }

    // 主运行循环
    void run(int generations, int eliteCount, double mutationRate, std::function<double(int)> sa_T_start,
             std::function<double(double, int, int,int)> sa_nextTemp,std::function<double(int)> crossoverRatio) {
        auto startTime = std::chrono::steady_clock::now();
        BS::thread_pool pool;
        std::uniform_real_distribution<double> probGen(0.0, 1.0);
        for (int genIdx = 0; genIdx < generations; ++genIdx) {
            std::vector<Individual> newPop;
            newPop.resize(popSize);

            // 1. 精英保留
            for (int i = 0; i < eliteCount; ++i) {
                newPop[i] = population[i];
            }

            for (int pos = eliteCount; pos < popSize; pos++) {
                pool.detach_task([this,&newPop, crossoverRatio,genIdx,mutationRate,&probGen, sa_T_start, sa_nextTemp, pos, eliteCount]() {
                    auto &plan = *plans[pos];
                    // 2. 交叉与变异
                    int q1=tournamentSelect(3),q2=tournamentSelect(3);
                    while(q1==q2)
                        q2=tournamentSelect(3);
                    const auto &p1 = population[q1];
                    const auto &p2 = population[q2];

                    Individual &child = newPop[pos];
                    child.mapping = crossoverGPX(p1.mapping, p2.mapping,crossoverRatio(genIdx));
                    if (genIdx >= 3 && probGen(gen) < mutationRate)
                        mutate(child.mapping);

                    // 3. 局部搜索 (SA) - 拉马克进化
                    evaluate(child, plan);
                    int cnt = plan.encoding.dupCnt;
                    auto sa_nextTemp_withGenIdx=[sa_nextTemp,genIdx](double T, int step, int best){
                        return sa_nextTemp(T,step,best,genIdx);
                    };
                    plan.sa.solve(sa_T_start(cnt), 0, sa_nextTemp_withGenIdx, pos - eliteCount + 1);
                    child.mapping = plan.mapping.radicalToKey;
                    child.fitness = plan.encoding.dupCnt;
                    std::cout << std::format("Pos: {} | DupCnt: {} -> {}  ...\n", pos, cnt, child.fitness);
                });
            }
            pool.wait();
            population = std::move(newPop);
            std::sort(population.begin(), population.end(),
                      [](const Individual &a, const Individual &b) { return a.fitness < b.fitness; });
            if (population[0].fitness < bestInd.fitness) {
                bestInd = population[0];
            }
            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();
            std::cout << std::format("[Gen {} | Best: {} | Avg: {:.1f} | Time: {}s]\n", genIdx + 1, bestInd.fitness,
                                     calculateAvgFitness(), elapsed);
        }

        // 最后将最优解应用到 encoding
        encoding.mapping.radicalToKey = bestInd.mapping;
        encoding.createEncodingFromMapping();
        encoding.buildHash();
    }

private:
    encodingData<max_length> &encoding;
    int numKeys;
    int preAllocRadical;
    int popSize;
    std::vector<Individual> population;
    std::vector<std::unique_ptr<optimizationPlan>> plans;
    Individual bestInd;
    std::mt19937 gen;

    // 评估函数，重建哈希表
    void evaluate(Individual &ind, optimizationPlan &plan) {
        plan.mapping.radicalToKey = ind.mapping;
        plan.encoding.createEncodingFromMapping();
        plan.encoding.buildHash();
    }

    int tournamentSelect(int k) {
        std::uniform_int_distribution<int> dist(0, popSize - 1);
        int bestIdx = dist(gen);
        for (int i = 1; i < k; ++i) {
            int idx = dist(gen);
            if (population[idx].fitness < population[bestIdx].fitness) {
                bestIdx = idx;
            }
        }
        return bestIdx;
    }

    void mutate(std::vector<uint16_t> &mapping) {
        std::uniform_int_distribution<int> radDist(preAllocRadical + 1, encoding.numRadical);
        std::uniform_int_distribution<int> keyDist(1, numKeys);
        int mutations = 15;
        for (int i = 0; i < mutations; ++i) {
            mapping[radDist(gen)] = keyDist(gen);
        }
    }

    double calculateAvgFitness() {
        double sum = 0;
        for (const auto &ind : population)
            sum += ind.fitness;
        return sum / popSize;
    }

    // 交叉算子，输入两个字根到键的映射方案，输出交叉结果
    std::vector<uint16_t> crossoverGPX(const std::vector<uint16_t> &p1, const std::vector<uint16_t> &p2,double crossoverRatio) {
        std::vector<uint16_t> childMap(encoding.numRadical + 1, 0);
        std::vector<bool> radicalAssigned(encoding.numRadical + 1, false);

        // 辅助结构：预计算每个父本中 [按键 -> 字根列表]
        // 及其对应的剩余未分配字根数量 (动态更新)
        std::vector<std::vector<uint16_t>> keyToRadical1(numKeys + 1), keyToRadical2(numKeys + 1);
        std::vector<int> cnt1(numKeys + 1, 0), cnt2(numKeys + 1, 0);
        std::vector<double> weight1(numKeys + 1, 0), weight2(numKeys + 1, 0);

        for (int r = 1; r <= encoding.numRadical; ++r) {
            keyToRadical1[p1[r]].push_back(r);
            cnt1[p1[r]]++;
            weight1[p1[r]] += encoding.mapping.radicalWeight[r];
            keyToRadical2[p2[r]].push_back(r);
            cnt2[p2[r]]++;
            weight2[p2[r]] += encoding.mapping.radicalWeight[r];
        }

        // 存储提取出来的字根组
        std::vector<std::vector<uint16_t>> fixedGroups;    // 包含固定字根的组
        std::vector<std::vector<uint16_t>> floatingGroups; // 不含固定字根的组（原固定字根已被删除）

        std::uniform_real_distribution<double> probGen(0.0, 1.0);
        int extractedKeysCount = 0;

        // --- 第一阶段：贪心抽取 ---
        while (extractedKeysCount < numKeys) {
            // A. 随机选择父本
            bool useP1 = (probGen(gen)<=crossoverRatio);
            auto &curCnt = useP1 ? cnt1 : cnt2;
            auto &otherCnt = useP1 ? cnt2 : cnt1;
            auto &curWeight = useP1 ? weight1 : weight2;
            auto &otherWeight = useP1 ? weight2 : weight1;
            const auto &curKeyToRadical = useP1 ? keyToRadical1 : keyToRadical2;
            const auto &otherRadicalToKey = useP1 ? p2 : p1;

            // B. 贪心选择剩余字根最多的按键
            int bestKey = 0;
            int cntRecord = 0;
            double maxWeight = 0;
            for (int k = 1; k <= numKeys; ++k) {
                if (curWeight[k] > maxWeight) {
                    cntRecord = curCnt[k];
                    maxWeight = curWeight[k];
                    bestKey = k;
                }
            }
            // 如果没有剩余字根可取，提前结束
            if (bestKey == 0)
                break;

            // C. 提取该按键下的所有未删除字根
            std::vector<uint16_t> group;
            group.reserve(cntRecord);
            bool hasFixedRadical = false;

            for (auto r : curKeyToRadical[bestKey]) {
                if (!radicalAssigned[r]) {
                    radicalAssigned[r] = true;
                    group.push_back(r);

                    // 检查是否包含固定字根
                    if (r <= preAllocRadical) {
                        hasFixedRadical = true;
                    }
                    otherCnt[otherRadicalToKey[r]]--;
                    otherWeight[otherRadicalToKey[r]] -= encoding.mapping.radicalWeight[r];
                }
            }
            // 当前选中的按键计数归零
            curCnt[bestKey] = 0;
            curWeight[bestKey] = 0;

            // D. 分类暂存
            if (hasFixedRadical) {
                fixedGroups.push_back(std::move(group));
            }
            else {
                floatingGroups.push_back(std::move(group));
            }
            extractedKeysCount++;
        }

        // --- 3. 顺序分配 ---
        int currentKey = 1;

        // A. 先分配带固定字根的组
        for (const auto &grp : fixedGroups) {
            for (auto r : grp)
                childMap[r] = currentKey;
            currentKey++;
        }

        // 从keyForFloatingGroups开始的按键，目前都没有固定字根
        int keyForFloatingGroups = currentKey;

        // B. 继续往后填不带固定字根的组
        for (const auto &grp : floatingGroups) {
            for (auto r : grp)
                childMap[r] = currentKey;
            currentKey++;
        }

        // --- 4. 处理剩余固定字根  ---
        std::vector<uint16_t> leftoverFixed;
        for (int r = 1; r <= preAllocRadical; ++r) {
            if (!radicalAssigned[r])
                leftoverFixed.push_back(r);
        }

        if (!leftoverFixed.empty()) {
            std::shuffle(leftoverFixed.begin(), leftoverFixed.end(), gen);
            // 从刚才记录的位置开始填，因为之前的键位已经被 fixedGroups 占了固定名额
            // 由于preAllocRadical<=numKey，直接往后分配
            int fixKeyIter = keyForFloatingGroups;
            for (auto r : leftoverFixed) {
                childMap[r] = fixKeyIter;
                radicalAssigned[r] = true;
                fixKeyIter++;
            }
        }

        // --- 5. 处理其余所有未分配字根---
        std::uniform_int_distribution<int> keyDist(1, numKeys);
        for (int r = preAllocRadical + 1; r <= encoding.numRadical; ++r) {
            if (!radicalAssigned[r]) {
                childMap[r] = keyDist(gen);
            }
        }
        return childMap;
    }
};