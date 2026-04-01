#include "BS_thread_pool.hpp"
#include "data_numerical_hashonly.h"
#include "pcg_random.hpp"
#include "platform.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

template <typename OperatorType, typename OperatorResultType, int maxLength, int numKeys>
class
    AMHBOptimizer // OperatorType为std::variant，应为返回OperatorResultType的仿函数；OperatorResultType为std::variant，应具有deltaE、index、apply成员
{
public:
    struct AMHBParameters
    {
        int maxIterations;
        double TempStart;
    };

    struct AMHBWorker;

    class LiteOperatorPoolType
    {
    public:
        LiteOperatorPoolType(int totalNeighbors) : totalNeighbors(totalNeighbors) {
            taskBuffer.resize(totalNeighbors);
            EXP3Discount.resize(totalNeighbors);
        }

        void addOperator(OperatorType op, std::string name, double cost) {
            operators.push_back(op);
            nameConfig.push_back(name);
            costConfig.push_back(cost);
            tempCoef.push_back(1);
            taskCount.push_back(0);
            weight.push_back(1);
            VaR.push_back(1);
            EXP3Weight.push_back(1);
        }

        void updateStats(int index, int deltaE) {
            int optIdx = taskBuffer[index];
            // 更根VaR_(1/2)(delta E^+)
            if (deltaE > 0)
                VaR[optIdx] =
                    std::max(std::numeric_limits<double>::epsilon(), VaR[optIdx] + ((deltaE > VaR[optIdx]) ? etaVaR : -etaVaR));
            // 更新EXP3权重
            const double &r = 1.0 / (1.0 + (deltaE < 0 ? std::exp(deltaE) : std::exp(deltaE / VaR[0]))) / weight[optIdx];
            EXP3Weight[optIdx] += etaEXP3 * r * EXP3Discount[index];
        }

        int cal(std::vector<std::unique_ptr<AMHBWorker>> &allWorkers) {
            // LogSumExp计算基础权重
            double maxL = -std::numeric_limits<double>::infinity();
            for (const auto &i : EXP3Weight)
                maxL = std::max(maxL, i);
            double sumExp = 0.0, totalWeight = 0.0;
            for (size_t i = 0; i < operators.size(); ++i) {
                weight[i] = std::exp(EXP3Weight[i] - maxL);
                sumExp += weight[i];
            }
            for (size_t i = 0; i < operators.size(); ++i) {
                weight[i] = (1.0 - gamma) * (weight[i] / sumExp) + (gamma / operators.size());
                totalWeight += weight[i];
            }
            // 多项分布抽样
            pcg32 gen(std::random_device{}());
            double remainingWeight = totalWeight;
            int remainingSamples = totalNeighbors;
            for (size_t i = 0; i < operators.size(); ++i) {
                if (remainingSamples <= 0) {
                    taskCount[i] = 0;
                    continue;
                }
                if (i == operators.size() - 1)
                    taskCount[i] = remainingSamples;
                else {
                    double p = weight[i] / remainingWeight;
                    p = std::clamp(p, 0.0, 1.0);
                    std::binomial_distribution<int> dist(remainingSamples, p);
                    taskCount[i] = dist(gen);
                }
                remainingSamples -= taskCount[i];
                remainingWeight -= weight[i];
            }
            // 任务生成
            for (int i = 0, pos = 0; i < operators.size(); ++i) {
                double curDiscount = 1;
                for (int j = 0; j < taskCount[i]; ++j, ++pos, curDiscount *= lambda) {
                    taskBuffer[pos] = i;
                    EXP3Discount[pos] = curDiscount;
                }
                EXP3Weight[i] *= curDiscount;
            }
            // 任务分配
            // std::shuffle(taskBuffer.begin(),taskBuffer.end(),gen);
            int totalTasks = taskBuffer.size();
            int numWorkers = allWorkers.size();
            int baseCount = totalTasks / numWorkers;
            int extra = totalTasks % numWorkers;
            int currentOffset = 0;
            for (int i = 0; i < numWorkers; ++i) {
                int assigned = baseCount + (i < extra ? 1 : 0);
                allWorkers[i]->left.store(currentOffset, std::memory_order_relaxed);
                allWorkers[i]->right.store(currentOffset + assigned, std::memory_order_relaxed);
                currentOffset += assigned;
            }
            return totalTasks;
        }
        std::vector<OperatorType> operators;
        std::vector<int> taskBuffer;

        int totalNeighbors;
        std::vector<std::string> nameConfig;
        std::vector<double> costConfig;
        std::vector<double> tempCoef;
        std::vector<int> taskCount;
        std::vector<double> weight;
        std::vector<double> VaR;          // VaR估值
        std::vector<double> EXP3Discount; // EXP3遗忘权重
        std::vector<double> EXP3Weight;   // EXP3对数权重
    private:
        double alpha = 0;
        double gamma = 0.05;    // EXP3 探索率
        double lambda = 0.9999; // EXP3 遗忘因子
        double etaEXP3 = 0.02;  // EXP3 学习率
        double etaVaR = 0.01;   // VaR学习率
    };

    struct AMHBWorker
    {
        AMHBWorker(int stealThreshold, std::atomic<bool> &allowRun, std::atomic<int> &usedWorker, std::atomic<int> &workerDone,
                   std::atomic<int> &totalTask, std::atomic<bool> &allTaskComplete, LiteOperatorPoolType &liteOperatorPool,
                   std::vector<std::unique_ptr<AMHBWorker>> &workerPointer, const size_t bufferSize, const double &T)
            : stealThreshold(stealThreshold), globalAllowRun(allowRun), globalUsedWorker(usedWorker),
              globalWorkerDone(workerDone), globalTotalTask(totalTask), allTaskComplete(allTaskComplete),
              liteOperatorPool(liteOperatorPool), allWorkers(workerPointer), liteOperatorRst(bufferSize), T(T),
              localGen(std::random_device{}()) {}
        void run() {
            while (true) {
                // 等待允许开始指令
                while (!globalAllowRun.load(std::memory_order_acquire)) {
                    if (allTaskComplete.load(std::memory_order_acquire))
                        return;
                    cpu_pause();
                }
                // 开始运行
                int currentLeft = left.load(std::memory_order_seq_cst);
                int currentRight = right.load(std::memory_order_seq_cst);
                // TODO: 重算子使用无锁队列
                // 没有任务
                if (currentRight - currentLeft == 0) {
                    // 直接跳到等allowRun变false 的阶段
                    while (globalAllowRun.load(std::memory_order_acquire))
                        cpu_pause();
                    continue;
                }
                // 存在任务
                globalUsedWorker.fetch_add(1, std::memory_order_relaxed);
                availableRst = 0;
                maxGumbelScore = -std::numeric_limits<double>::infinity();
                // 循环运行和偷窃
                for (bool checkAllowRun = true; checkAllowRun;) {
                    runLocalTasks();
                    while (!runSteal()) {
                        if (!globalAllowRun.load(std::memory_order_acquire)) {
                            checkAllowRun = false;
                            break;
                        }
                    }
                }
                // 汇报本线程局部结算完成
                globalWorkerDone.fetch_add(1, std::memory_order_release);
            }
        }
        // 线程数据
        mapping<maxLength, numKeys> encoding;
        pcg32 localGen;
        // Gumbel-Max数据
        double maxGumbelScore;
        OperatorResultType *bestCandidate = nullptr;
        // 运行数据
        int availableRst = 0;
        alignas(64) std::atomic<int> left{0};
        alignas(64) std::atomic<int> right{0};
        alignas(64) std::atomic<bool> allowSteal{false};
        alignas(64) std::atomic_flag stealLock = ATOMIC_FLAG_INIT;
        int localLeft;
        int pos;
        bool localAllowSteal;
        std::vector<std::unique_ptr<AMHBWorker>> &allWorkers;
        // 控制变量
        int stepSize;             // 禁止偷取阈值
        const int stealThreshold; // 禁止偷窃阈值
        const double &T;          // 全局温度
        std::atomic<bool> &globalAllowRun;
        std::atomic<int> &globalUsedWorker;
        std::atomic<int> &globalWorkerDone;
        std::atomic<int> &globalTotalTask;
        std::atomic<bool> &allTaskComplete;
        // 算子数据
        LiteOperatorPoolType &liteOperatorPool;
        std::vector<OperatorResultType> liteOperatorRst;

    private:
        int getStepSize(int l, int r) {
            const int &n = r - l;
            return std::max(1, std::min(32, static_cast<int>(std::floor(std::sqrt(n)))));
        }
        void runLocalTasks() {
            int localCompleted = 0;
            // pos为当前任务，localLeft、left为已锁定的任务，stepSize根据left和right初始化
            pos = left.load(std::memory_order_seq_cst);
            int localRight = right.load(std::memory_order_seq_cst);
            stepSize = getStepSize(pos, localRight);
            // 判断是否拥有足够的任务以允许偷窃
            if (localRight - pos >= stepSize + stealThreshold) {
                // 有足够的任务以供偷取
                localLeft = std::min(pos + stepSize, localRight);
                left.store(localLeft, std::memory_order_seq_cst);
                localAllowSteal = true;
                allowSteal.store(true, std::memory_order_release);
            }
            else {
                // 没有足够的任务以供偷取，锁定剩下所有任务，关闭偷窃
                localLeft = localRight;
                left.store(localLeft, std::memory_order_seq_cst);
                localAllowSteal = false;
                allowSteal.store(false, std::memory_order_release);
            }
            // 开始运行
            while (true) {
                // 执行pos处的任务
                liteOperatorRst[availableRst] =
                    std::visit([this](auto &opt) -> OperatorResultType { return opt(encoding, pos, localGen); },
                               liteOperatorPool.operators[liteOperatorPool.taskBuffer[pos]]);
                // 计算Gumbel-Max Trick
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                double gumbelNoise = -std::log(-std::log(std::clamp(dist(localGen), std::numeric_limits<double>::min(),
                                                                    1.0 - std::numeric_limits<double>::epsilon())));
                const int &deltaE = std::visit([](auto &v) { return v.deltaE; }, liteOperatorRst[availableRst]);
                const int &index = std::visit([](auto &v) { return v.index; }, liteOperatorRst[availableRst]);
                double score =
                    -(static_cast<double>(deltaE) / (liteOperatorPool.tempCoef[liteOperatorPool.taskBuffer[index]] * T)) +
                    gumbelNoise;
                if (score > maxGumbelScore) {
                    maxGumbelScore = score;
                    bestCandidate = &(liteOperatorRst[availableRst]);
                }
                // 修改状态
                availableRst++;
                localCompleted++;
                // 自增pos以获取下一个任务
                pos++;
                // 若pos>=left，则开始更改left
                if (pos >= localLeft) {
                    auto lockFinalBatch = [&]() {
                        // 加锁并等待加锁完成
                        while (stealLock.test_and_set(std::memory_order_acquire))
                            cpu_pause();
                        // 赋值left=right以锁定剩下所有任务，并更改allowSteal为false
                        localRight = right.load(std::memory_order_seq_cst);
                        localLeft = localRight;
                        left.store(localLeft, std::memory_order_seq_cst);
                        localAllowSteal = false;
                        allowSteal.store(false, std::memory_order_relaxed);
                        stealLock.clear(std::memory_order_release);
                        return pos != localLeft;
                    };
                    localRight = right.load(std::memory_order_seq_cst);
                    stepSize = getStepSize(pos, localRight);
                    // 余量足够以当前步长自增，且不会达到锁定阀值
                    if (localRight - localLeft >= stepSize + stealThreshold) {
                        localLeft += stepSize;
                        left.store(localLeft, std::memory_order_seq_cst);
                        // 重新判断right-left>=k
                        if (right.load(std::memory_order_seq_cst) - localLeft < stealThreshold)
                            if (!lockFinalBatch())
                                break;
                    }
                    // 余量不够，锁定最后一批
                    else if (localAllowSteal) {
                        if (!lockFinalBatch())
                            break;
                    }
                    // 所有任务都已经完成
                    else
                        break;
                }
            }
            // 判断本次任务是不是最后一批
            if (localCompleted > 0) {
                int prevTotal = globalTotalTask.fetch_sub(localCompleted, std::memory_order_acq_rel);
                if (prevTotal == localCompleted)
                    globalAllowRun.store(false, std::memory_order_release);
            }
        }
        bool runSteal() {
            const int &numWorkers = allWorkers.size();
            if (numWorkers <= 1)
                return false;
            // 随机挑选一个受害者
            std::uniform_int_distribution<int> dist(0, numWorkers - 1);
            int targetIdx = dist(localGen);
            AMHBWorker &victim = *(allWorkers[targetIdx]);
            // 不能偷自己
            if (&victim == this)
                return false;
            // 检查受害者是否允许被偷
            if (!victim.allowSteal.load(std::memory_order_relaxed))
                return false;
            // 尝试加锁，如果锁已经被占用，直接放弃
            if (victim.stealLock.test_and_set(std::memory_order_acquire))
                return false;
            // 再次确认allowSteal，防止在第一次检查到加锁间，受害者已经更改allowSteal
            if (!victim.allowSteal.load(std::memory_order_relaxed)) {
                victim.stealLock.clear(std::memory_order_release);
                return false;
            }
            // 加锁成功，获取数据，此处的值只用来计算偷取范围，可以用放松内存序
            int vLeft = victim.left.load(std::memory_order_relaxed);
            int vRight = victim.right.load(std::memory_order_relaxed);
            // 不够一次偷窃
            if (vRight - vLeft < 2 * stealThreshold || vRight - vLeft < 1) {
                victim.stealLock.clear(std::memory_order_release);
                return false;
            }
            int middle = std::min((vLeft + vRight + 1) / 2, vRight - 1);
            // 修改right，偷取后半部分
            victim.right.store(middle, std::memory_order_seq_cst);
            // 回滚检查：重新判断 left <= middle - k
            int newLeft = victim.left.load(std::memory_order_seq_cst);
            if (newLeft > middle - stealThreshold) {
                // 即使此时受害者如果察觉到干扰，会阻塞在自旋锁上，这里可以使用放松内存序回滚
                victim.right.store(vRight, std::memory_order_relaxed);
                victim.stealLock.clear(std::memory_order_release);
                return false;
            }
            // 成功截获了[middle, vRight)
            this->left.store(middle, std::memory_order_seq_cst);
            this->right.store(vRight, std::memory_order_seq_cst);
            // 释放受害者的锁
            victim.stealLock.clear(std::memory_order_release);
            return true;
        }
    };

    AMHBOptimizer(mapping<maxLength, numKeys> &ed, int numWorker, bool print, int totalNeighbors, int stealThreshold)
        : print(print), LiteOperatorPool(totalNeighbors), allowRun(false), allTaskComplete(false) {
        for (int i = 0; i < numWorker; ++i) {
            // 将所有需要的全局变量和引用传入 Worker
            workers.push_back(std::make_unique<AMHBWorker>(stealThreshold, allowRun, usedWorker, workerDone, totalTask,
                                                           allTaskComplete, LiteOperatorPool, workers, totalNeighbors,
                                                           this->T));
            workers.back()->encoding = ed;
        }
    }
    template <typename nextTempType> void solve(AMHBParameters param, nextTempType nextTemp) {
        auto &encoding = workers[0]->encoding;
        if (print) {
            std::cout << "\n AMHB | Optimization Start. Initial Duplicates: " << encoding.dupCnt << std::endl;
        }
        auto startTime = std::chrono::steady_clock::now();
        T = param.TempStart;
        // 线程杂拉起所有worker
        BS::thread_pool pool(workers.size());
        for (auto &workerPtr : workers) {
            pool.detach_task([&workerPtr]() { workerPtr->run(); });
        }
        // 准备数据
        int bestEnergy = encoding.dupCnt, energyBefore = encoding.dupCnt, energyCur = encoding.dupCnt;
        bestMapping = encoding.radicalToKey;
        // 主循环
        for (int iter = 1; iter <= param.maxIterations; ++iter) {
            // 准备任务
            int totalGeneratedTasks = LiteOperatorPool.cal(workers);
            OperatorResultType *bestCandidate = nullptr;
            double maxGumbelScore = -std::numeric_limits<double>::infinity();
            totalTask.store(totalGeneratedTasks, std::memory_order_release);
            usedWorker.store(0, std::memory_order_relaxed);
            workerDone.store(0, std::memory_order_relaxed);
            // 开始计算任务，所有worker察觉allowRun==true后脱离自旋
            allowRun.store(true, std::memory_order_release);
            // 等待所有邻域计算和窃取结束，最后一个完成任务的worker会将allowRun改为false
            while (allowRun.load(std::memory_order_acquire))
                cpu_pause();
            // 等待有效worker完成局部Gumbel-Max采样，只有实际执行过任务的worker才会进入usedWorker
            while (usedWorker.load(std::memory_order_acquire) != workerDone.load(std::memory_order_acquire))
                cpu_pause();
            // 结果汇总、算子分数更新
            OperatorResultType *globalBestCandidate = nullptr;
            double globalMaxScore = -std::numeric_limits<double>::infinity();
            for (int thread = 0; thread < workers.size(); thread++) {
                if (workers[thread]->bestCandidate != nullptr) {
                    if (workers[thread]->maxGumbelScore > globalMaxScore) {
                        globalMaxScore = workers[thread]->maxGumbelScore;
                        globalBestCandidate = workers[thread]->bestCandidate;
                    }
                    workers[thread]->bestCandidate = nullptr;
                    for (int i = 0; i < workers[thread]->availableRst; i++) {
                        const int &deltaE = std::visit([](auto &v) { return v.deltaE; }, workers[thread]->liteOperatorRst[i]);
                        const int &index = std::visit([](auto &v) { return v.index; }, workers[thread]->liteOperatorRst[i]);
                        LiteOperatorPool.updateStats(index, deltaE);
                    }
                }
            }
            // 应用选中更改
            if (globalBestCandidate != nullptr)
                std::visit([this](auto &rst) { return rst.template apply<std::vector<std::unique_ptr<AMHBWorker>>>(workers); },
                           *globalBestCandidate);
            energyBefore = energyCur;
            energyCur = encoding.dupCnt;
            if (energyCur < bestEnergy) {
                bestEnergy = energyCur;
                bestMapping = encoding.radicalToKey;
            }
            // 输出
            if (print && iter % 10000 == 0) {
                auto now = std::chrono::steady_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count() / 1000.0;
                std::cout << std::format("\n[Iter: {} | Temp: {} | Best: {} | Current: {} | Time: {:.2f}s ]", iter, T,
                                         bestEnergy, energyCur, elapsed)
                          << std::endl;
                for (int i = 0; i < LiteOperatorPool.operators.size(); i++)
                    std::cout << std::format("[ {}: tempCoef: {:.4f} VaR: {:.4f} Prob: {:.4f} ] ", i,
                                             LiteOperatorPool.tempCoef[i], LiteOperatorPool.VaR[i], LiteOperatorPool.weight[i]);
                std::cout << std::endl;
            }
            T = nextTemp(T, iter, energyCur);
            if (T < 0)
                break;
        }
        allTaskComplete.store(true, std::memory_order_release);
        encoding.radicalToKey = bestMapping;
        encoding.createEncodingFromMapping();
        encoding.buildHash();
        encoding.writeMapping("test_data/niwang_sanma/output/AMHB/");
        if (print)
            std::cout << "\n AMHB | Optimization Complete. Best: " << bestEnergy << std::endl;
        pool.wait();
    }
    LiteOperatorPoolType LiteOperatorPool;

private:
    double T;
    std::vector<std::unique_ptr<AMHBWorker>> workers;
    std::atomic<int> totalTask;
    std::atomic<int> usedWorker;
    std::atomic<int> workerDone;
    std::atomic<bool> allowRun;
    std::atomic<bool> allTaskComplete;
    bool print;
    std::vector<uint16_t> bestMapping;
};

struct PointwiseModificationResultType
{
    int deltaE;
    int index;
    int radIdx = 0;
    uint16_t newKey = 0;
    template <typename WorkerList> void apply(WorkerList &workers) {
        for (auto &w : workers)
            w->encoding.modifyRadical(radIdx, newKey);
    }
};

template <int maxLength, int numKeys> struct PointwiseModificationOperator
{
    int preAllocRadical;

    PointwiseModificationOperator(int preAlloc) : preAllocRadical(preAlloc) {}

    PointwiseModificationResultType operator()(mapping<maxLength, numKeys> &encoding, int taskIndex, pcg32 &gen) const {
        std::uniform_int_distribution<int> radicalGenNoPreAlloc(preAllocRadical + 1, encoding.numRadical);
        std::uniform_int_distribution<int> keyGen(1, numKeys);

        // 1. 随机选择一个字根和新按键
        int radIdx = radicalGenNoPreAlloc(gen);
        uint16_t oldKey = encoding.radicalToKey[radIdx];
        uint16_t newKey = keyGen(gen);
        while (oldKey == newKey) {
            newKey = keyGen(gen);
        }

        // 2. 尝试修改并计算能量差 (deltaE)
        int currentEnergy = encoding.dupCnt;
        encoding.modifyRadical(radIdx, newKey);
        int deltaE = encoding.dupCnt - currentEnergy;

        // 3. 撤销修改 (在AMHB机制中，算子仅负责探测邻域，并不直接应用状态)
        encoding.modifyRadical(radIdx, oldKey);

        // 返回一个符合 AMHB 要求的样本结果
        return PointwiseModificationResultType{deltaE, taskIndex, radIdx, newKey};
    }
};

struct ExchangeModificationResultType
{
    int deltaE;
    int index;
    int radIdx1 = 0;
    int radIdx2 = 0;
    template <typename WorkerList> void apply(WorkerList &workers) {
        for (auto &w : workers) {
            auto tmp = w->encoding.radicalToKey[radIdx2];
            w->encoding.modifyRadical(radIdx2, w->encoding.radicalToKey[radIdx1]);
            w->encoding.modifyRadical(radIdx1, tmp);
        }
    }
};

template <int maxLength, int numKeys> struct ExchangeModificationOperator
{
    int preAllocRadical;
    ExchangeModificationOperator(int preAlloc) : preAllocRadical(preAlloc) {}
    ExchangeModificationResultType operator()(mapping<maxLength, numKeys> &encoding, int taskIndex, pcg32 &gen) const {
        std::uniform_int_distribution<int> radicalGenNoPreAlloc(preAllocRadical + 1, encoding.numRadical);
        std::uniform_int_distribution<int> radicalGenWithPreAlloc(1, encoding.numRadical);
        std::uniform_int_distribution<int> radicalGenOnlyPreAlloc(1, preAllocRadical);
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

        auto exchangeKey = [&](uint16_t a, uint16_t b) {
            auto tmp = encoding.radicalToKey[b];
            encoding.modifyRadical(b, encoding.radicalToKey[a]);
            encoding.modifyRadical(a, tmp);
        };

        int currentEnergy = encoding.dupCnt;
        exchangeKey(radIdx1, radIdx2);
        int deltaE = encoding.dupCnt - currentEnergy;
        exchangeKey(radIdx1, radIdx2);
        return ExchangeModificationResultType{deltaE, taskIndex, radIdx1, radIdx2};
    }
};

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

    // 5. 注册算子
    optimizer.LiteOperatorPool.addOperator(PointwiseModificationOperator<3, 26>(preAllocRadical), "Pointwise", 1.0);
    optimizer.LiteOperatorPool.addOperator(ExchangeModificationOperator<3, 26>(preAllocRadical), "Exchange", 2.0);

    // 6. 设置运行参数
    AMHBOptimizer<OpType, OpResult, 3, 26>::AMHBParameters param{
        100000000, // maxIterations: AMHB 是以固定迭代次数驱动的
        40         // TempStart
    };

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

    // 8. 启动优化
    optimizer.solve<nextTempType>(param, nextTemp);

    return 0;
}