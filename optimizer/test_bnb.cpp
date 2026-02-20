#include "include/data.h"
#include "include/branch_and_bound.h"

int main() {
    mappingData<3> mp;
    mp.readRadicalEncoding("testdata-triple/3d_data.txt");

    int preAllocRadical = mp.numRadical-8;
    std::random_device rd;
    // std::mt19937 gen(rd());
    std::mt19937 gen(0);
    std::uniform_int_distribution<int> dist(1, 26);
    std::vector<int> useRadical;
    for (int i = 1; i <= mp.numRadical; i++) {
        if (i > preAllocRadical) {
            useRadical.push_back(i);
        }
        else{
            mp.radicalToKey[i] = dist(gen);
        }
    }
    encodingData<3> encoding(mp);
    encoding.createEncodingFromMapping();
    encoding.buildHash();

    BnbOptimizer<3> optimizer(encoding, 26, 1000,true);
    optimizer.solve(useRadical);
    encoding.writeKeyEncoding("testdata-triple/output");
}