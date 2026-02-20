#include "include/data.h"
#include "include/tabu.h"

int main() {
    mappingData<3> mp;
    mp.readRadicalEncoding("testdata-triple/3d_data.txt");
    // mp.readMapping("testdata-triple/26-02-17-06-26-01.4573490-mapping.txt");
    std::mt19937 gen(0);
    std::uniform_int_distribution<int> dist(1, 26);

    int preAllocRadical = 26;
    for (int i = 1; i <= mp.numRadical; i++) {
        if (i <= preAllocRadical)
            mp.radicalToKey[i] = i;
        else
            mp.radicalToKey[i] = dist(gen);
    }

    encodingData<3> encoding(mp);
    encoding.createEncodingFromMapping();
    encoding.buildHash();

    TabuOptimizer<3> ts(encoding, 26, 26, true);
    ts.solve(50000, 15, 200, 100);

    return 0;
}