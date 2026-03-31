#include "BS_thread_pool.hpp"
#include "data_numerical_hashonly.h"
#include "sa.h"

int main() {
    BS::thread_pool pool;
    for (int i = 1; i <= 400; i++) {
        pool.detach_task([i]() {
            mapping<3, 26> encoding;
            encoding.readRadicalEncoding("test_data/niwang_sanma/radical_data.txt");

            std::mt19937 gen(std::random_device{}());
            std::uniform_int_distribution<int> dist(1, 26);

            int preAllocRadical = 26;
            for (int i = 1; i <= encoding.numRadical; i++) {
                if (i <= preAllocRadical)
                    encoding.radicalToKey[i] = i;
                else
                    encoding.radicalToKey[i] = dist(gen);
            }

            auto startTime = std::chrono::steady_clock::now();
            encoding.createEncodingFromMapping();
            encoding.buildHash();
            int cnt = encoding.dupCnt;

            SAOptimizer<3, 26> optimizer(encoding, preAllocRadical, false);
            SAOptimizer<3, 26>::SAParameters param{
                22.5,
            };
            struct nextTempType
            {
                double operator()(double T, int numItr, double currentEnergy) {
                    if (T > 7.5)
                        return T * 0.999999;
                    else if (T > 1.5)
                        return T * 0.9999999;
                    else if (T > 0.25)
                        return T * 0.9999999625;
                    else if (T > 0.1) {
                        return T * 0.9999999;
                    }
                    else if (T > 0.005) {
                        return T * 0.999999;
                    }
                    else {
                        return -1;
                    }
                }
            } nextTemp;
            optimizer.solve<nextTempType>(param, nextTemp);
            optimizer.applyBest();
            encoding.writeMapping("testdata-triple/output/sa3");
            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
            std::cout << std::format("Pos: {} | DupCnt: {} -> {} | Time: {:.3f}s\n", i, cnt, encoding.dupCnt,elapsed/1000.0);
        });
    }
    pool.wait();
    return 0;
}
