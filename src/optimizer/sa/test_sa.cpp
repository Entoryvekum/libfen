#include "data_numerical_hashonly.h"
#include "sa.h"

int main() {
    mapping<3, 26> encoding;
    encoding.readRadicalEncoding("test_data/niwang_sanma/radical_data.txt");

    std::mt19937 gen(std::random_device{}());
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
    SAOptimizer<3, 26> optimizer(encoding, preAllocRadical, true);
    SAOptimizer<3, 26>::SAParameters param{
        22.5,
        0,
    };
    
    struct nextTempType
    {
        double operator()(double T, int numItr, double currentEnergy) {
            if (currentEnergy < 120 && !below120) {
                below120 = true;
                return 0.75;
            }
            if(below120==true&&cntBelow120<=1000000) {
                cntBelow120++;
                return T;
            }
            if (T > 7.5)
                return T * 0.999999;
            else if (T > 1.5)
                return T * 0.9999999;
            else if (T > 0.25)
                return T * 0.999999975;
            else if (T > 0.1) {
                return T * 0.99999995;
            }
            else if (T > 0.005) {
                return T * 0.999999;
            }
            else {
                cntBelow120=0;
                below120=false;
                return -1;
            }
        }
        int cntBelow120=0;
        bool below120 = false;
    } nextTemp;
    optimizer.solve<nextTempType>(param, nextTemp);
    return 0;
}
