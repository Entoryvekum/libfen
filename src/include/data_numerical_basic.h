#pragma once
#include <array>
#include <chrono>
#include <cstdint>
#include <format>
#include <fstream>
#include <set>
#include <sstream>
#include <vector>

template <int maxLength, int numKeys> class mapping
{
public:
    mapping() {
        int maxEncodingSpace = 1;
        for (int i = 0; i < maxLength; i++)
            maxEncodingSpace *= (numKeys + 1);
        hashTable.assign(maxEncodingSpace, {});
    }
    void readRadicalEncoding(const std::string &filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        std::string line;
        radicalEncoding.clear();
        numRadical = 0;
        while (std::getline(file, line)) {
            if (line.empty())
                continue;
            std::istringstream iss(line);
            std::array<uint16_t, maxLength> radicalRow{};
            int value;
            for (int count = 0; iss >> value && count < maxLength; count++) {
                radicalRow[count] = static_cast<uint16_t>(value);
                if (value > numRadical) {
                    numRadical = value;
                }
            }
            radicalEncoding.push_back(radicalRow);
        }
        keyEncoding.assign(radicalEncoding.size(), {});
        radicalToKey.assign(numRadical + 1, {});
        radicalToHanzi.assign(numRadical + 1, {});
        std::set<uint16_t> radicalSet;
        for (int i = 0; i < radicalEncoding.size(); i++) {
            radicalSet.clear();
            for (int j = 0; j < maxLength; j++) {
                if (radicalEncoding[i][j] != 0) {
                    if (radicalSet.find(radicalEncoding[i][j]) == radicalSet.end()) {
                        radicalToHanzi[radicalEncoding[i][j]].push_back(i);
                        radicalSet.insert(radicalEncoding[i][j]);
                    }
                }
                else
                    break;
            }
        }
    }
    void readMapping(const std::string &filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        int i = 1;
        for (; i <= numRadical && (file >> radicalToKey[i]); i++)
            ;
        if (i - 1 != numRadical)
            throw std::runtime_error("Number of Encodings Do not Match: " + filename);
    }
    void writeMapping(const std::string &directory) {
        auto now = std::chrono::system_clock::now();
        std::string filename = std::format("{}/{:%y-%m-%d-%H-%M-%S}-d{}-mapping.txt", directory, now, dupCnt);

        std::ofstream file(filename);
        if (!file.is_open())
            throw std::runtime_error("Cannot create file: " + filename);

        for (int i = 1; i < radicalToKey.size(); i++) {
            file << radicalToKey[i] << "\n";
        }
    }
    void createEncodingFromMapping() {
        for (int i = 0; i < radicalEncoding.size(); i++) {
            for (int j = 0; j < maxLength; j++) {
                keyEncoding[i][j] = radicalToKey[radicalEncoding[i][j]];
            }
        }
    }
    uint32_t hashIndex(const std::array<uint16_t, maxLength> &encoding) {
        uint32_t rst = encoding[0];
        for (int i = 1; i < maxLength; i++) {
            rst = rst * (numKeys + 1) + encoding[i];
        }
        return rst;
    }
    void buildHash() {
        dupCnt = 0;
        std::fill(hashTable.begin(), hashTable.end(), 0);
        for (int i = 0; i < keyEncoding.size(); i++) {
            auto &v = hashTable[hashIndex(keyEncoding[i])];
            if (v > 0)
                dupCnt++;
            v++;
        }
    }
    void modifyRadical(int radicalIdx, uint16_t newKey) {
        radicalToKey[radicalIdx] = newKey;
        const auto &affectedHanzi = radicalToHanzi[radicalIdx];
        for (int hanziIdx : affectedHanzi) {
            auto tmp = hashIndex(keyEncoding[hanziIdx]);
            dupCnt -= (--hashTable[tmp] > 0);
            for (int j = 0; j < maxLength; ++j) {
                // 无分支修改
                const bool flag = (radicalEncoding[hanziIdx][j] == radicalIdx);
                keyEncoding[hanziIdx][j] = (newKey * flag) | (keyEncoding[hanziIdx][j] * (!flag));
            }
            tmp = hashIndex(keyEncoding[hanziIdx]);
            dupCnt += (++hashTable[tmp] > 1);
        }
    }
    int numRadical;
    int dupCnt;
    std::vector<uint16_t> radicalToKey;                           // 字根对应的字母
    std::vector<std::array<uint16_t, maxLength>> radicalEncoding; // 汉字的字根编码
    std::vector<std::array<uint16_t, maxLength>> keyEncoding;     // 汉字的字母编码
    std::vector<std::vector<int>> radicalToHanzi;                 // 包合某个字根的汉字
    std::vector<int> hashTable;                                   // 汉字哈希桶
};