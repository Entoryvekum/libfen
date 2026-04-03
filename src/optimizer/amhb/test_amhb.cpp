#include "BS_thread_pool.hpp"
#include "data_numerical_hashonly.h"
#include "amhb.h"

using OpResult = std::variant<PointwiseModificationResultType, ExchangeModificationResultType>;
using OpType = std::variant<PointwiseModificationOperator<3, 26>, ExchangeModificationOperator<3, 26>>;

int main() {
    mapping<3, 26> encoding;
    encoding.readRadicalEncoding("test_data/niwang_sanma/radical_data.txt");

    pcg32 gen(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, 25);

    int preAllocRadical = 26;
    for (int i = 1; i <= encoding.numRadical; i++) {
        if (i <= preAllocRadical)
            encoding.radicalToKey[i] = i;
        else
            encoding.radicalToKey[i] = dist(gen);
    }
    encoding.createEncodingFromMapping();
    encoding.buildHash();

    int numWorker = 4;
    int totalNeighbors = 240;
    int stealThreshold = 0;

    AMHBOptimizer<OpType, OpResult, 3, 26> optimizer(encoding, numWorker, true, totalNeighbors, stealThreshold);

    // 注册算子
    optimizer.LiteOperatorPool.addOperator(PointwiseModificationOperator<3, 26>(preAllocRadical), "Pointwise", 1.0);
    optimizer.LiteOperatorPool.addOperator(ExchangeModificationOperator<3, 26>(preAllocRadical), "Exchange", 2.0);

    // 设置运行参数
    AMHBOptimizer<OpType, OpResult, 3, 26>::AMHBParameters param{100000000, 40};

    struct nextTempType
    {
        double operator()(double T, int numItr, double currentEnergy) {
            if (T > 7.5)
                return T * 0.99999;
            else if (T > 1.75)
                return T * 0.999999;
            else if (T > 1.25)
                return T * 0.9999995;
            else if (T > 0.35)
                return T * 0.9999999625;
            else {
                return -1; // 标志位
            }
        }
    } nextTemp;

    // 启动优化
    optimizer.solve<nextTempType>(param, nextTemp);

    return 0;
}