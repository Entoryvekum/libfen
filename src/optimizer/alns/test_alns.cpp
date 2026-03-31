#include "BS_thread_pool.hpp"
#include "data_numerical_basic.h"
#include "alns.h"
#include "local_sa.h"
#include "local_tabu.h"
#include <algorithm>

int main() {
    mapping<3, 26> encoding;
    encoding.readRadicalEncoding("test_data/niwang_sanma/radical_data.txt");
    encoding.readMapping("test_data/niwang_sanma/output/26-02-19-03-37-15.6752058-mapping.first_below_90.txt");

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> dist(1, 26);

    int preAllocRadical = 26;
    // 初始化尚未分配的字根
    for (int i = 1; i <= encoding.numRadical; i++) {
        if (i <= preAllocRadical)
            encoding.radicalToKey[i] = i;
        else
            encoding.radicalToKey[i] = dist(gen);
    }
    encoding.createEncodingFromMapping();
    encoding.buildHash();

    ALNSOptimizer<3, 26>::ALNSParameters alnsParam{500, 25.0, 0.99, 20,
                                                   0.3, 0.2,  0.15, [](double t, int c) {
            return 0.2*std::sqrt(t)+c; }};
    ALNSOptimizer<3, 26> alns(encoding, true);

    std::vector<int> available;
    for (int i = 26; i <= encoding.numRadical; ++i) {
        available.push_back(i);
    }

    auto Random_Destroy_20 = [preAllocRadical, &available](mapping<3, 26> &enc) {
        static std::mt19937 rng(std::random_device{}());
        std::shuffle(available.begin() + preAllocRadical, available.end(), rng);
        std::vector<int> destroyIdx(available.begin() + preAllocRadical, available.begin() + preAllocRadical + 20);
        std::uniform_int_distribution<int> dist(1, 26);
        for (auto i : destroyIdx)
            enc.modifyRadical(i, dist(rng));
        return destroyIdx;
    };

    auto Random_Destroy_40 = [preAllocRadical, &available](mapping<3, 26> &enc) {
        static std::mt19937 rng(std::random_device{}());
        std::shuffle(available.begin() + preAllocRadical, available.end(), rng);
        std::vector<int> destroyIdx(available.begin() + preAllocRadical, available.begin() + preAllocRadical + 40);
        std::uniform_int_distribution<int> dist(1, 26);
        for (auto i : destroyIdx)
            enc.modifyRadical(i, dist(rng));
        return destroyIdx;
    };

    auto Random_Destroy_80 = [preAllocRadical, &available](mapping<3, 26> &enc) {
        static std::mt19937 rng(std::random_device{}());
        std::shuffle(available.begin() + preAllocRadical, available.end(), rng);
        std::vector<int> destroyIdx(available.begin() + preAllocRadical, available.begin() + preAllocRadical + 80);
        std::uniform_int_distribution<int> dist(1, 26);
        for (auto i : destroyIdx)
            enc.modifyRadical(i, dist(rng));
        return destroyIdx;
    };

    auto Random_Destroy_120 = [preAllocRadical, &available](mapping<3, 26> &enc) {
        static std::mt19937 rng(std::random_device{}());
        std::shuffle(available.begin() + preAllocRadical, available.end(), rng);
        std::vector<int> destroyIdx(available.begin() + preAllocRadical, available.begin() + preAllocRadical + 120);
        std::uniform_int_distribution<int> dist(1, 26);
        for (auto i : destroyIdx)
            enc.modifyRadical(i, dist(rng));
        return destroyIdx;
    };

    auto Greedy_Destroy_10_40 = [preAllocRadical, &available](mapping<3, 26> &enc) -> std::vector<int> {
        std::vector<uint64_t> contribution(enc.numRadical + 1, 0);
        for (size_t i = 0; i < enc.keyEncoding.size(); ++i) {
            uint32_t idx = enc.hashIndex(enc.keyEncoding[i]);
            if (enc.hashTable[idx] > 1)
                for (int j = 0; j < 3; ++j)
                    contribution[enc.radicalEncoding[i][j]] += enc.hashTable[idx];
        }
        std::vector<std::pair<int, uint64_t>> sortedRadicals;
        for (int i = 1; i <= enc.numRadical; ++i)
            if (contribution[i] > 0 && i > preAllocRadical)
                sortedRadicals.push_back({i, contribution[i]});
        std::sort(sortedRadicals.begin(), sortedRadicals.end(),
                  [](const auto &a, const auto &b) { return a.second > b.second; });
        std::vector<int> destroyIdx;
        for (int i = 0; i < std::min(size_t(10), sortedRadicals.size()); i++)
            destroyIdx.push_back(sortedRadicals[i].first);
        static std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(1, 26);
        std::shuffle(available.begin() + preAllocRadical, available.end(), rng);
        for (int i = 0; i < 30; i++)
            destroyIdx.push_back(available[i + preAllocRadical]);
        for (auto i : destroyIdx)
            enc.modifyRadical(i, dist(rng));
        return destroyIdx;
    };

    auto Greedy_Destroy_20_80 = [preAllocRadical, &available](mapping<3, 26> &enc) -> std::vector<int> {
        std::vector<uint64_t> contribution(enc.numRadical + 1, 0);
        for (size_t i = 0; i < enc.keyEncoding.size(); ++i) {
            uint32_t idx = enc.hashIndex(enc.keyEncoding[i]);
            if (enc.hashTable[idx] > 1)
                for (int j = 0; j < 3; ++j)
                    contribution[enc.radicalEncoding[i][j]] += enc.hashTable[idx];
        }
        std::vector<std::pair<int, uint64_t>> sortedRadicals;
        for (int i = 1; i <= enc.numRadical; ++i)
            if (contribution[i] > 0 && i > preAllocRadical)
                sortedRadicals.push_back({i, contribution[i]});
        std::sort(sortedRadicals.begin(), sortedRadicals.end(),
                  [](const auto &a, const auto &b) { return a.second > b.second; });
        std::vector<int> destroyIdx;
        for (int i = 0; i < std::min(size_t(20), sortedRadicals.size()); i++)
            destroyIdx.push_back(sortedRadicals[i].first);
        static std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(1, 26);
        std::shuffle(available.begin() + preAllocRadical, available.end(), rng);
        for (int i = 0; i < 60; i++)
            destroyIdx.push_back(available[i + preAllocRadical]);
        for (auto i : destroyIdx)
            enc.modifyRadical(i, dist(rng));
        return destroyIdx;
    };

    auto Greedy_Swap_10_40 = [preAllocRadical, &available](mapping<3, 26> &enc) -> std::vector<int> {
        std::vector<uint64_t> contribution(enc.numRadical + 1, 0);
        for (size_t i = 0; i < enc.keyEncoding.size(); ++i) {
            uint32_t idx = enc.hashIndex(enc.keyEncoding[i]);
            if (enc.hashTable[idx] > 1)
                for (int j = 0; j < 3; ++j)
                    contribution[enc.radicalEncoding[i][j]] += enc.hashTable[idx];
        }
        std::vector<std::pair<int, uint64_t>> sortedRadicals;
        for (int i = 1; i <= enc.numRadical; ++i)
            if (contribution[i] > 0 && i > preAllocRadical)
                sortedRadicals.push_back({i, contribution[i]});
        std::sort(sortedRadicals.begin(), sortedRadicals.end(),
                  [](const auto &a, const auto &b) { return a.second > b.second; });
        std::vector<int> destroyIdx;
        for (int i = 0; i < std::min(size_t(10), sortedRadicals.size()); i++)
            destroyIdx.push_back(sortedRadicals[i].first);
        static std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(1, 26);
        std::shuffle(available.begin() + preAllocRadical, available.end(), rng);
        auto exchangeKey = [&](uint16_t a, uint16_t b) {
            auto tmp = enc.radicalToKey[b];
            enc.modifyRadical(b, enc.radicalToKey[a]);
            enc.modifyRadical(a, tmp);
        };
        for (int i = 0; i < destroyIdx.size(); i++)
            exchangeKey(destroyIdx[i], available[i + preAllocRadical]);
        for (int i = 0; i < 30; i++)
            destroyIdx.push_back(available[i + preAllocRadical]);
        return destroyIdx;
    };

    auto Greedy_Swap_20_80 = [preAllocRadical, &available](mapping<3, 26> &enc) -> std::vector<int> {
        std::vector<uint64_t> contribution(enc.numRadical + 1, 0);
        for (size_t i = 0; i < enc.keyEncoding.size(); ++i) {
            uint32_t idx = enc.hashIndex(enc.keyEncoding[i]);
            if (enc.hashTable[idx] > 1)
                for (int j = 0; j < 3; ++j)
                    contribution[enc.radicalEncoding[i][j]] += enc.hashTable[idx];
        }
        std::vector<std::pair<int, uint64_t>> sortedRadicals;
        for (int i = 1; i <= enc.numRadical; ++i)
            if (contribution[i] > 0 && i > preAllocRadical)
                sortedRadicals.push_back({i, contribution[i]});
        std::sort(sortedRadicals.begin(), sortedRadicals.end(),
                  [](const auto &a, const auto &b) { return a.second > b.second; });
        std::vector<int> destroyIdx;
        for (int i = 0; i < std::min(size_t(20), sortedRadicals.size()); i++)
            destroyIdx.push_back(sortedRadicals[i].first);
        static std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(1, 26);
        std::shuffle(available.begin() + preAllocRadical, available.end(), rng);
        auto exchangeKey = [&](uint16_t a, uint16_t b) {
            auto tmp = enc.radicalToKey[b];
            enc.modifyRadical(b, enc.radicalToKey[a]);
            enc.modifyRadical(a, tmp);
        };
        for (int i = 0; i < destroyIdx.size(); i++)
            exchangeKey(destroyIdx[i], available[i + preAllocRadical]);
        for (int i = 0; i < 60; i++)
            destroyIdx.push_back(available[i + preAllocRadical]);
        return destroyIdx;
    };

    auto Local_SA_Repair_0_5 = [](mapping<3, 26> &enc, const std::vector<int> &destroyedRadicals) -> int {
        std::vector<mapping<3, 26>> encodings(18, enc);
        BS::thread_pool pool;
        for (int i = 0; i < 16; i++) {
            pool.detach_task([&encodings, &enc, &destroyedRadicals, i]() {
                LocalSAOptimizer<3, 26> saOpt(encodings[i], false);
                typename LocalSAOptimizer<3, 26>::SAParameters params{30};
                struct nextTempType
                {
                    double operator()(double T, int numItr, double currentEnergy) {
                        if (T > 7.5)
                            return T * 0.999999;
                        else if (T > 1.5)
                            return T * 0.99999;
                        else if (T > 0.25)
                            return T * 0.99999;
                        else if (T > 0.1)
                            return T * 0.999975;
                        else if (T > 0.005)
                            return T * 0.9999;
                        else
                            return -1;
                    }
                } nextTemp;
                saOpt.solve<nextTempType>(params, nextTemp, destroyedRadicals);
            });
        }
        pool.wait();
        int bestIdx = 0;
        for (int i = 0; i < 18; i++) {
            if (encodings[i].dupCnt < enc.dupCnt)
                enc = encodings[i];
        }
        return enc.dupCnt;
    };

    auto Local_SA_Repair_1 = [](mapping<3, 26> &enc, const std::vector<int> &destroyedRadicals) -> int {
        std::vector<mapping<3, 26>> encodings(18, enc);
        BS::thread_pool pool;
        for (int i = 0; i < 16; i++) {
            pool.detach_task([&encodings, &enc, &destroyedRadicals, i]() {
                LocalSAOptimizer<3, 26> saOpt(encodings[i], false);
                typename LocalSAOptimizer<3, 26>::SAParameters params{30};
                struct nextTempType
                {
                    double operator()(double T, int numItr, double currentEnergy) {
                        if (T > 7.5)
                            return T * 0.999999;
                        else if (T > 1.5)
                            return T * 0.999995;
                        else if (T > 0.25)
                            return T * 0.99999925;
                        else if (T > 0.1)
                            return T * 0.99999;
                        else if (T > 0.005)
                            return T * 0.9999;
                        else
                            return -1;
                    }
                } nextTemp;
                saOpt.solve<nextTempType>(params, nextTemp, destroyedRadicals);
            });
        }
        pool.wait();
        int bestIdx = 0;
        for (int i = 0; i < 18; i++) {
            if (encodings[i].dupCnt < enc.dupCnt)
                enc = encodings[i];
        }
        return enc.dupCnt;
    };

    auto Local_SA_Repair_2 = [](mapping<3, 26> &enc, const std::vector<int> &destroyedRadicals) -> int {
        std::vector<mapping<3, 26>> encodings(18, enc);
        BS::thread_pool pool;
        for (int i = 0; i < 16; i++) {
            pool.detach_task([&encodings, &enc, &destroyedRadicals, i]() {
                LocalSAOptimizer<3, 26> saOpt(encodings[i], false);
                typename LocalSAOptimizer<3, 26>::SAParameters params{20};
                struct nextTempType
                {
                    double operator()(double T, int numItr, double currentEnergy) {
                        if(T>0.21)
                            return 18252303454 * std::exp(-20.6318 * std::pow(static_cast<double>(numItr), 1.0/80));
                        else
                            return -1;
                    }
                } nextTemp;
                saOpt.solve<nextTempType>(params, nextTemp, destroyedRadicals);
            });
        }
        pool.wait();
        int bestIdx = 0;
        for (int i = 0; i < 18; i++) {
            if (encodings[i].dupCnt < enc.dupCnt)
                enc = encodings[i];
        }
        return enc.dupCnt;
    };

    auto Local_TABU_Repair = [](mapping<3, 26> &enc, const std::vector<int> &destroyedRadicals) -> int {
        LocalTabuOptimizer<3, 26> tabuOpt(enc, false);
        typename LocalTabuOptimizer<3, 26>::TabuParameters params{
            5000,
            200,
            500,
        };
        tabuOpt.solve(params, destroyedRadicals);
        tabuOpt.applyBest();
        return enc.dupCnt;
    };
    alns.addDestroyOperator(Random_Destroy_20, "Random_Destroy_20", 1);
    alns.addDestroyOperator(Random_Destroy_40, "Random_Destroy_40", 2);
    alns.addDestroyOperator(Random_Destroy_80, "Random_Destroy_80", 4, 2);
    alns.addDestroyOperator(Random_Destroy_120, "Random_Destroy_120", 6.5, 3);
    alns.addDestroyOperator(Greedy_Destroy_10_40, "Greedy_Destroy_10_40", 3, 1);
    alns.addDestroyOperator(Greedy_Destroy_20_80, "Greedy_Destroy_20_80", 6, 3);
    alns.addDestroyOperator(Greedy_Swap_10_40, "Greedy_Swap_10_40", 3, 1);
    alns.addDestroyOperator(Greedy_Swap_20_80, "Greedy_Swap_20_80", 6, 3);
    alns.addRepairOperator(Local_SA_Repair_0_5, "Local_SA_Repair_0_5", 0.5);
    alns.addRepairOperator(Local_SA_Repair_1, "Local_SA_Repair_1", 1.5);
    alns.addRepairOperator(Local_SA_Repair_2, "Local_SA_Repair_2", 5);
    alns.addRepairOperator(Local_TABU_Repair, "Local_TABU_Repair", 1);

    alns.solve(alnsParam);

    // 可在此处添加保存最终 encoding 映射的代码
    // encoding.writeMapping("...");

    return 0;
}