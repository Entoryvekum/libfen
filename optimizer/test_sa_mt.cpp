#include "include/data.h"
#include "include/sa.h"
#include "BS_thread_pool.hpp"

int main() {
    BS::thread_pool pool;
    for(int i=1;i<=1000;i++) {
        pool.detach_task([i]() {
            mappingData<3> mp;
            mp.readRadicalEncoding("testdata-triple/3d_data.txt");
            mp.readMapping("testdata-triple/26-02-17-06-26-01.4573490-mapping.txt");

            int preAllocRadical = 26;

            encodingData<3> encoding(mp);
            encoding.createEncodingFromMapping();
            encoding.buildHash();
            int cnt=encoding.dupCnt;

            SAOptimizer<3> optimizer(encoding, 26, preAllocRadical, false);
            SAOptimizer<3>::SAParameters param{
                1, 
                0, 
            };
            struct nextTempType{
                double operator()(double T,int numItr,double bestEnergy){
                    if (T > 0.2) {
                        return T * 0.9999999;
                    }
                    else {
                        return -1;
                    }
                }
            } nextTemp;
            optimizer.solve<nextTempType>(param,nextTemp);
            std::cout << std::format("Pos: {} | DupCnt: {} -> {} \n", i, cnt, encoding.dupCnt);
        });
    }
    return 0;
}
