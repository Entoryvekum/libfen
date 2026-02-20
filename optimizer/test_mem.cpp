#include "include/data.h"
#include "include/memetic.h"
#include <iostream>
#include <random>
#include <cmath>

int main() {
    mappingData<3> mp;
    mp.readRadicalEncoding("testdata-triple/3d_data.txt");
    mp.readMapping("testdata-triple/26-02-17-06-26-01.4573490-mapping.txt");

    for(auto &v:mp.radicalWeight) {
        v=std::log(v*0.1+1)*0.5+1;
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

    // 3. 配置模因算法参数
    int popSize = 25;       // 种群大小
    int generations = 10000;  // 进化代数
    int eliteCount = 4;     // 精英保留数量
    double mutationRate=0.4;
    auto crossoverRatio=[&gen](int genIdx){
        if(genIdx<=2)
            return 0.5;
        else if(genIdx<=3) {
            std::uniform_real_distribution<double> probGen(0.6, 0.8);
            return probGen(gen);
        }
        else {
            std::uniform_real_distribution<double> probGen(0.75, 0.85);
            return probGen(gen);
        }
    };

    // 实例化优化器
    MemeticOptimizer<3> optimizer(encoding, numKeys, preAllocRadical, popSize);

    // 初始化种群
    optimizer.initialize();

    // SA PARAMETER
    auto localSearchSchedule = [](double T, int step, int best,int genIdx) -> double {
        static int cnt=0;
        if(T>1)
            return T * 0.99999925;
        else if(T>0.2)
            return T * 0.9999999;
        else {
            if(best<95&&cnt<2) {
                ++cnt;
                return T+0.15;
            }
            else if(cnt==2&&T>0.175) {
                return T*0.9999999;
            }
            else {
                cnt=0;
                return -1;
            }
        }
    };
    auto sa_T_start = [](int dup) -> double {
        if(dup>500)
            return 5;
        else if(dup>300)
            return 4;
        else if(dup>200)
            return 3;
        else if(dup>150)
            return 1.5;
        else if(dup>110)
            return 1;
        else
            return 0.85;
    };

    // 4. 运行优化
    std::cout << "Starting Memetic Optimization..." << std::endl;
    optimizer.run(
        generations, 
        eliteCount, 
        mutationRate,
        sa_T_start,
        localSearchSchedule,
        crossoverRatio
    );

    // 5. 输出最终结果
    // run() 结束时会将最优解应用回 encoding 对象
    std::cout << "Optimization Finished." << std::endl;
    std::cout << "Final Duplicate Count: " << encoding.dupCnt << std::endl;
    
    // 可选：保存结果
    // mp.writeMapping(".");
    // encoding.writeKeyEncoding(".");

    return 0;
}