#include "../include/branch_and_bound.h"
#include "../include/data.h"
#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

void LNS(encodingData<3> &encoding, mappingData<3> &mp, int preAllocRadical, std::mt19937 &gen) {
    auto startTime = std::chrono::steady_clock::now();

    std::vector<int> candidates;
    for (int i = preAllocRadical + 1; i <= mp.numRadical; ++i) {
        candidates.push_back(i);
    }

    int k = 3;
    int round = 0;
    double T = 0.5;
    int best = 120;
    std::uniform_real_distribution<double> probGen(0.0, 1.0);
    std::uniform_int_distribution<int> radicalGen(preAllocRadical + 1, encoding.numRadical);
    std::uniform_int_distribution<int> keyGen(1, 26);
    int minorReheat = 0;
    int majorReheat = 0;

    auto nextTemp = [&minorReheat, &majorReheat](double &x) {
        if(x>0.275||x<0.1) {
            x *= 0.9999995;
        }
        else {
            x*=0.999999725;
        }
        if (x <= 0.05) {
            if (minorReheat <= 5) {
                x += 0.25;
                minorReheat++;
            }
            else {
                x += 0.4;
                minorReheat = 0;
                majorReheat++;
            }
        }
    };

    while (true) {
        round++;
        int dupBefore = encoding.dupCnt;

        // SA
        bool flag = false;
        for (int numItr = 1; numItr <= encoding.numRadical * 20; numItr++) {
            uint16_t radIdx = radicalGen(gen);
            uint16_t oldKey = encoding.mapping.radicalToKey[radIdx];
            uint16_t newKey = keyGen(gen);
            while (oldKey == newKey)
                newKey = keyGen(gen);
            encoding.modifyMapping(radIdx, newKey);
            int deltaE = encoding.dupCnt - dupBefore;
            if (deltaE < 0 || probGen(gen) < std::exp(-deltaE / T)) {
                dupBefore = encoding.dupCnt;
                flag = true;
                if (encoding.dupCnt < best)
                    best = encoding.dupCnt;
            }
            else {
                encoding.modifyMapping(radIdx, oldKey);
            }
            nextTemp(T);
        }
        for (int numItr = 1; numItr <= encoding.numRadical * 20; numItr++) {
            uint16_t radIdx1,radIdx2;
            auto getRandomPair=[&radIdx1,&radIdx2,&radicalGen,&gen](){
                radIdx1 = radicalGen(gen),radIdx2=radicalGen(gen);
                while(radIdx2==radIdx1)
                    radIdx2=radicalGen(gen);
            };
            while(encoding.mapping.radicalToKey[radIdx1]==encoding.mapping.radicalToKey[radIdx2])
                getRandomPair();
            auto exchangeKey=[&encoding](uint16_t a,uint16_t b) {
                auto tmp=encoding.mapping.radicalToKey[b];
                encoding.modifyMapping(b, encoding.mapping.radicalToKey[a]);
                encoding.modifyMapping(a, tmp);
            };
            exchangeKey(radIdx1,radIdx2);
            int deltaE = encoding.dupCnt - dupBefore;
            if (deltaE < 0 || probGen(gen) < std::exp(-deltaE / T)) {
                dupBefore = encoding.dupCnt;
                flag = true;
                if (encoding.dupCnt < best)
                    best = encoding.dupCnt;
            }
            else {
                exchangeKey(radIdx1,radIdx2);
            }
            nextTemp(T);
        }
        std::shuffle(candidates.begin(), candidates.end(), gen);
        for (int numItr = 0; flag && numItr < 6; numItr++) {
            // 备份旧按键，并将这些字根清零
            std::vector<int> useRadical(candidates.begin() + numItr * k, candidates.begin() + (numItr + 1) * k);
            std::vector<int> backupKeys;
            backupKeys.reserve(k);
            for (int r : useRadical) {
                backupKeys.push_back(mp.radicalToKey[r]);
                encoding.modifyMapping(r, 0);
            }

            // B&B
            BnbOptimizer<3> bnb(encoding, 26, dupBefore, false);
            bool success = bnb.solve(useRadical);

            // 根据结果更新或回滚
            if (success) {
                for (int r : useRadical) {
                    encoding.modifyMapping(r, bnb.bestMapping[r]);
                }
                dupBefore = encoding.dupCnt;
                if (encoding.dupCnt < best) {
                    best = encoding.dupCnt;
                }
            }
            else {
                for (size_t i = 0; i < useRadical.size(); ++i) {
                    encoding.modifyMapping(useRadical[i], backupKeys[i]);
                }
            }
        }
        if (round % 10 == 0) {
            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();
            std::cout
                << std::format(
                       "\r[T: {:.4f} | Round {} | Best: {} | Current: {} | Minor Reheat: {} | Major Reheat: {} | Time: {}s]", T,
                       round, best, encoding.dupCnt, minorReheat, majorReheat, elapsed)
                << std::flush;
        }
        if (encoding.dupCnt <= 100) {
            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();
            std::cout << "\nOptimization Finished" << std::endl;
            std::cout << std::format("\r[T: {:.4f} | Round {} | Current: {} | Time: {}s]", T, round, encoding.dupCnt, elapsed)
                      << std::flush;
            break;
        }
        
    }
}

int main() {
    mappingData<3> mp;
    mp.readRadicalEncoding("testdata-triple/3d_data.txt");
    mp.readMapping("testdata-triple/26-02-17-06-26-01.4573490-mapping.txt");

    // std::mt19937 gen(std::random_device{}());
    std::mt19937 gen(0);
    std::uniform_int_distribution<int> dist(1, 26);

    int preAllocRadical = 26;
    encodingData<3> encoding(mp);
    encoding.createEncodingFromMapping();
    encoding.buildHash();

    LNS(encoding, mp, 26, gen);

    return 0;
}
