#include "data_numerical_hashonly.h"
#include "tabu.h"

int main() {
    mapping<3,26> encoding;
    encoding.readRadicalEncoding("test_data/niwang_sanma/radical_data.txt");
    std::mt19937 gen(0);
    std::uniform_int_distribution<int> dist(1, 26);

    int preAllocRadical = 26;
    for (int i = 1; i <= encoding.numRadical; i++) {
        if (i <= preAllocRadical)
            encoding.radicalToKey[i] = i;
        else
            encoding.radicalToKey[i] = dist(gen);
    }

    encoding.createEncodingFromMapping();
    encoding.buildHash();

    TabuOptimizer<3,26>::TabuParameters param{
        500000, 
        200, 
        500, 
        0
    };

    TabuOptimizer<3,26> ts(encoding, preAllocRadical, true);
    ts.solve(param);

    return 0;
}