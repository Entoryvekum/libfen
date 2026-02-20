#include "include/data.h"
#include "include/sa.h"

int main() {
    mappingData<3> mp;
    mp.readRadicalEncoding("testdata-triple/3d_data.txt");

    std::mt19937 gen(std::random_device{}());
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

    SAOptimizer<3> optimizer(encoding, 26, preAllocRadical, true);
    optimizer.solve(20.0,  95, [](double T, int step, int best) {
        static int cnt = 0;
        if (T > 1)
            return T * 0.99999925;
        else if (T > 0.5) {
            return T * 0.9999999;
        }
        else if (T > 0.2) {
            return T * 0.9999999;
        }
        else {
            if (cnt < 10) {
                T += 0.25;
                cnt++;
            }
            else {
                cnt = 0;
                T += 0.5;
            }
            return T * 0.9999999;
        }
    },0);
    return 0;
}
