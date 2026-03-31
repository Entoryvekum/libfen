#pragma once
#include "data_numerical_basic.h"
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <chrono>
#include <functional>
#include <format>
#include <string>
#include <algorithm>

const double SCORE_BEST = 10.0, SCORE_BETTER = 5.0, SCORE_ACCEPT = 2.0;

template <int maxLength, int numKeys> 
class ALNSOptimizer {
public:
    using DestroyOp = std::function<std::vector<int>(mapping<maxLength,numKeys>&)>;
    using RepairOp = std::function<int(mapping<maxLength,numKeys>&,const std::vector<int>&)>;

    struct ALNSParameters {
        int maxIterations;
        double T_start;
        double coolingRate;
        int segmentLength;
        double reactionFactor;
        double strengthPenalty;
        double T_Delay_Coefficient;
        std::function<double(double,int)> costFunction;
    };

    ALNSOptimizer(mapping<maxLength, numKeys> &ed, bool print)
        : encoding(ed), print(print), gen(std::random_device{}()) {}

    void addDestroyOperator(DestroyOp op, std::string name, double strength=-1, int delayRounds = 0, int cooldownRounds = 0) {
        destroyPool.add(op, name, strength,delayRounds, cooldownRounds);
    }

    void addRepairOperator(RepairOp op, std::string name, double strength=-1, int delayRounds = 0, int cooldownRounds = 0) {
        repairPool.add(op, name, strength,delayRounds, cooldownRounds);
    }

    void solve(ALNSParameters param) {
        if (repairPool.empty() || destroyPool.empty()) {
            std::cerr << "ALNS | Error: Operators not fully initialized." << std::endl;
            return;
        }
        if (print) {
            std::cout << "\nALNS | Optimization Start. Initial Duplicates: " << encoding.dupCnt << std::endl;
        }

        bestMapping = encoding.radicalToKey;
        int currentEnergy = encoding.dupCnt;
        int bestEnergy = currentEnergy;
        double T = param.T_start;
        auto startTime = std::chrono::steady_clock::now();

        for (int iter = 1; iter <= param.maxIterations; ++iter) {
            auto backupMapping = encoding.radicalToKey;
            int energyBefore = currentEnergy;

            // 1. 自适应选择算子 (排除冷却中和等待结算的算子)
            int destIdx = destroyPool.select(gen,param);
            int repIdx = repairPool.select(gen,param,destroyPool.strengthConfig[destIdx]);

            // 2. 执行并计时：破坏算子
            auto t1 = std::chrono::steady_clock::now();
            std::vector<int> destroyedRadicals = destroyPool.ops[destIdx](encoding);
            auto t2 = std::chrono::steady_clock::now();
            double destroyOpTime=std::chrono::duration<double, std::milli>(t2 - t1).count()/1000;

            // 3. 执行并计时：修复算子
            auto t3 = std::chrono::steady_clock::now();
            int newEnergy = repairPool.ops[repIdx](encoding,destroyedRadicals);
            auto t4 = std::chrono::steady_clock::now();
            double repairTime=std::chrono::duration<double, std::milli>(t4 - t3).count()/1000;

            // 4. Metropolis接受准则
            bool accepted = false;
            if (newEnergy < bestEnergy) {
                bestEnergy = newEnergy;
                bestMapping = encoding.radicalToKey;
                accepted = true;
            } else if (newEnergy < currentEnergy) {
                accepted = true;
            } else {
                std::uniform_real_distribution<double> probDist(0.0, 1.0);
                double T1=destroyPool.delayConfig[destIdx]*param.T_Delay_Coefficient+T;
                if (probDist(gen) < std::exp((currentEnergy - newEnergy) / T1)) {
                    accepted = true;
                }
            }

            // 5. 状态更新或回滚
            if (accepted) {
                currentEnergy = newEnergy;
            } else {
                encoding.radicalToKey = backupMapping;
                encoding.createEncodingFromMapping();
                encoding.buildHash();
            }

            // 6. 维护池中已有算子的延迟状态，并结算到期的算子
            destroyPool.tick(currentEnergy, T, bestEnergy, gen);
            repairPool.tick(currentEnergy, T, bestEnergy, gen);

            // 7. 将本次使用的算子加入跟踪/冷却系统
            destroyPool.startTracking(destIdx, energyBefore, newEnergy, T, bestEnergy,accepted,destroyOpTime, gen);
            repairPool.startTracking(repIdx, energyBefore, newEnergy, T, bestEnergy, accepted,destroyOpTime, gen);

            // 8. 周期性更新权重
            if (iter % param.segmentLength == 0) {
                destroyPool.updateWeights(param.reactionFactor, param);
                repairPool.updateWeights(param.reactionFactor, param);
            }

            T *= param.coolingRate;
            if (print) {
                auto now = std::chrono::steady_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count() / 1000.0;
                std::cout << std::format("\n[Iter: {} | Best: {} | {} -> {} {} | Time: {:.2f}s | Operator: {}, {} ]", 
                             iter, bestEnergy, energyBefore,newEnergy,(accepted?"A":"N"),  elapsed, destroyPool.names[destIdx],repairPool.names[repIdx] ) << std::flush;
                if (iter % param.segmentLength == 0) {
                    std::cout << std::format("\nALNS | Segment {} Weights Update:", iter / param.segmentLength) << std::endl;
                    // 输出破坏算子权重
                    std::cout << "Destroy -> | ";
                    for (size_t i = 0; i < destroyPool.names.size(); ++i) {
                        std::cout << std::format("{}: {:.2f} | ", destroyPool.names[i], destroyPool.weights[i]);
                    }
                    // 输出修复算子权重
                    std::cout << "\nRepair  -> | ";
                    for (size_t i = 0; i < repairPool.names.size(); ++i) {
                        std::cout << std::format("{}: {:.2f} | ", repairPool.names[i], repairPool.weights[i]);
                    }
                    std::cout << std::endl;
                }
            }
        }

        encoding.radicalToKey = bestMapping;
        encoding.createEncodingFromMapping();
        encoding.buildHash();
        if (print) std::cout << "\nALNS | Optimization Complete. Best: " << bestEnergy << std::endl;
    }

private:
    mapping<maxLength, numKeys> &encoding;
    bool print;
    std::mt19937 gen;
    std::vector<uint16_t> bestMapping;

    template <typename OpType>
    struct OperatorPool {
        std::vector<OpType> ops;
        std::vector<std::string> names;
        
        // 配置数据
        std::vector<int> delayConfig;
        std::vector<int> cooldownConfig;
        std::vector<double> strengthConfig;
        
        // 核心评价数据
        std::vector<double> weights;
        std::vector<double> scores;
        std::vector<int> counts;
        std::vector<double> accumulatedTime;
        
        // 状态跟踪数据
        std::vector<int> currentDelay;
        std::vector<int> currentCooldown;
        std::vector<int> startEnergyRecord;
        std::vector<int> minEnergyRecord;
        std::vector<double> timeRecord;
        std::vector<double> tempRecord;

        bool empty() const { return ops.empty(); }
        
        void add(OpType op, std::string name, double strength, int delay, int cooldown) {
            ops.push_back(op);
            names.push_back(name);
            strengthConfig.push_back(strength);
            delayConfig.push_back(delay);
            cooldownConfig.push_back(cooldown);
            
            weights.push_back(1.0);
            scores.push_back(0.0);
            counts.push_back(0);
            accumulatedTime.push_back(0.0);
            
            currentDelay.push_back(0);
            currentCooldown.push_back(0);
            startEnergyRecord.push_back(0);
            minEnergyRecord.push_back(0);
            timeRecord.push_back(0);
            tempRecord.push_back(0);
        }

        int select(std::mt19937& rng,ALNSParameters &param,double trgtStrength=-1) {
            double total = 0;
            std::vector<int> available;
            
            // 仅筛选不在冷却中、且不在等待结算的算子
            auto penalty=[trgtStrength,param](double strength){
                if(trgtStrength>0)
                    return std::exp(-param.strengthPenalty*std::abs(trgtStrength-strength));
                else
                    return 1.0;
            };

            for (size_t i = 0; i < weights.size(); ++i) {
                if (currentDelay[i] == 0 && currentCooldown[i] == 0) {
                    total += weights[i]*penalty(strengthConfig[i]);
                    available.push_back(i);
                }
            }
            if (available.empty()) {
                throw std::runtime_error("No Available Operator");
            }

            std::uniform_real_distribution<double> dist(0.0, total);
            double r = dist(rng);
            double cum = 0;
            for (int idx : available) {
                cum += weights[idx]*penalty(strengthConfig[idx]);
                if (r <= cum)
                    return idx;
            }
            return available.back();
        }

        // 推进延迟轮次，并在到期时进行虚拟Metropolis结算
        void tick(int currentE, double currentT, int globalBestE, std::mt19937& rng) {
            for (size_t i = 0; i < ops.size(); ++i) {
                if (currentCooldown[i] > 0) {
                    currentCooldown[i]--;
                }
                if (currentDelay[i] > 0) {
                    // 更新追踪期内观察到的最小能量
                    minEnergyRecord[i] = std::min(minEnergyRecord[i], currentE);
                    currentDelay[i]--;
                    if (currentDelay[i] == 0) {
                        //在结算时计算代价，防止在运行后在结算新开segment
                        counts[i]++;
                        accumulatedTime[i]+=timeRecord[i];
                        evaluateScore(i, tempRecord[i], globalBestE, rng);
                    }
                }
            }
        }

        // 开始追踪新一轮被选中的算子
        void startTracking(int idx, int startE, int immediateE, double currentT, int globalBestE, bool accepted, double runningTime, std::mt19937& rng) {
            currentCooldown[idx] = cooldownConfig[idx];
            if(accepted&&delayConfig[idx] != 0) {// 只有被接受且算子有延迟结算时启动延迟追踪
                currentDelay[idx] = delayConfig[idx];
                startEnergyRecord[idx] = startE;
                minEnergyRecord[idx] = immediateE;
                timeRecord[idx]=runningTime;
                tempRecord[idx]=currentT;
            }
            else {
                counts[idx]++;
                accumulatedTime[idx]+=runningTime;
                if (immediateE < globalBestE) {
                    scores[idx] += SCORE_BEST;
                } else if (immediateE < startE) {
                    scores[idx] += SCORE_BETTER;
                } else if (accepted) {
                    scores[idx] += SCORE_ACCEPT;
                }
            }
        }

        // 虚拟 Metropolis 准则判定
        void evaluateScore(int idx, double currentT, int globalBestE, std::mt19937& rng) {
            double score = 0;
            int minE = minEnergyRecord[idx];
            int startE = startEnergyRecord[idx];
            
            if (minE < globalBestE) {
                score = SCORE_BEST;
            } else if (minE < startE) {
                score = SCORE_BETTER;
            } else {
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                if (dist(rng) < std::exp((startE - minE) / currentT)) {
                    score = SCORE_ACCEPT;
                }
            }
            scores[idx] += score;
        }

        void updateWeights(double reactionFactor, ALNSParameters &param) {
            for (size_t i = 0; i < ops.size(); ++i) {
                if (counts[i] > 0) {
                    double denominator = param.costFunction(accumulatedTime[i],counts[i]);
                    if (denominator > 0) { // 防除0
                        double averageScore = scores[i] / denominator;
                        weights[i] = (1 - reactionFactor) * weights[i] + reactionFactor * averageScore;
                    }
                }
                // 清空统计数据，迎接下一个 Segment
                scores[i] = 0;
                counts[i] = 0;
                accumulatedTime[i] = 0.0;
            }
        }
    };

    OperatorPool<DestroyOp> destroyPool;
    OperatorPool<RepairOp> repairPool;
};