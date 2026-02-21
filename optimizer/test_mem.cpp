#include "include/data.h"
#include "include/memetic.h"
#include <cmath>
#include <iostream>
#include <random>

int main() {
    mappingData<3> mp;
    mp.readRadicalEncoding("testdata-triple/3d_data.txt");
    mp.readMapping("testdata-triple/26-02-17-06-26-01.4573490-mapping.txt");

    for (auto &v : mp.radicalWeight) {
        v = std::log(v * 0.1 + 1) * 0.2 + 1;
        // v=1;
    }

    std::mt19937 gen(0);
    std::uniform_int_distribution<int> dist(1, 26);

    int preAllocRadical = 26;
    int numKeys = 26;

    // for (int i = 1; i <= mp.numRadical; i++) {
    //     if (i <= preAllocRadical)
    //         mp.radicalToKey[i] = i;
    //     else
    //         mp.radicalToKey[i] = dist(gen);
    // }

    // 2. 构建编码环境
    encodingData<3> encoding(mp);
    encoding.createEncodingFromMapping();
    encoding.buildHash();

    std::cout << "Initial Duplicate Count: " << encoding.dupCnt << std::endl;

    // 实例化优化器
    int popSize = 25;
    MemeticOptimizer<3> optimizer(encoding, numKeys, preAllocRadical, popSize);

    // 初始化种群
    optimizer.initialize();

    // 3. 配置模因算法参数
    struct nextTempType
    {
        double operator()(double T, int numItr, double bestEnergy) {
            if (T > 1)
                return T * 0.9999999;
            else if (T > 0.2) {
                if(genIdx<=3) {
                    if(genIdx==1&&T<0.5||genIdx==2&&T<0.4||genIdx==3&&T<0.3){
                        resetCnt();
                        return -1;
                    }
                    return T * 0.99999995;
                }
                else
                    return T * 0.9999999;
            }
            else {
                if (bestEnergy > 100) {
                    // 若最优解仍大于100，则重新升温
                    if (cntReheatForReject < 1) {
                        ++cntReheatForReject;
                        return T + 1;
                    }
                    else {
                        resetCnt();
                        return -1;
                    }
                }
                else if (bestEnergy < 95) {
                    // 若最优解小于95，则重新升温
                    if (cntReheatForRefine < 1) {
                        ++cntReheatForRefine;
                        return T + 0.5;
                    }
                    else if (cntReheatForRefine < 2 && genIdx > 3) {
                        ++cntReheatForRefine;
                        return T + 0.3;
                    }
                    else {
                        resetCnt();
                        return -1;
                    }
                }
                else {
                    resetCnt();
                    return -1;
                }
            }
        }
        void setGenIdx(int x) { genIdx = x; }
        void resetCnt() {
            cntReheatForReject = 0;
            cntReheatForRefine = 0;
        }
        int cntReheatForReject = 0;
        int cntReheatForRefine = 0;
        int genIdx;
    } SANextTemp;
    MemeticOptimizer<3>::memeParameter<nextTempType> param{10000,
                                                           4,
                                                           0.4,
                                                           [&gen](int genIdx) -> double {
                                                               if (genIdx <= 2)
                                                                   return 0.5;
                                                               else if (genIdx <= 3) {
                                                                   std::uniform_real_distribution<double> probGen(0.6, 0.8);
                                                                   return probGen(gen);
                                                               }
                                                               else {
                                                                   std::uniform_real_distribution<double> probGen(0.75, 0.85);
                                                                   return probGen(gen);
                                                               }
                                                           },
                                                           [](int dup) -> double {
                                                               if (dup > 500)
                                                                   return 5;
                                                               else if (dup > 300)
                                                                   return 4;
                                                               else if (dup > 200)
                                                                   return 3;
                                                               else if (dup > 150)
                                                                   return 1.5;
                                                               else if (dup > 110)
                                                                   return 1;
                                                               else
                                                                   return 0.85;
                                                           },
                                                           SANextTemp};

    // 4. 运行优化
    std::cout << "Starting Memetic Optimization..." << std::endl;
    optimizer.run(param);

    // 5. 输出最终结果
    // run() 结束时会将最优解应用回 encoding 对象
    std::cout << "Optimization Finished." << std::endl;
    std::cout << "Final Duplicate Count: " << encoding.dupCnt << std::endl;

    // 可选：保存结果
    // mp.writeMapping(".");
    // encoding.writeKeyEncoding(".");

    return 0;
}