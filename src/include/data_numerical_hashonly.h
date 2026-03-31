#pragma once
#include <array>
#include <chrono>
#include <cstdint>
#include <format>
#include <fstream>
#include <map>
#include <sstream>
#include <utility>
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
        numberEncoding.assign(radicalEncoding.size(), {});
        radicalToKey.assign(numRadical + 1, {});
        radicalToHanzi.assign(numRadical + 1, {});
        radicalMask.assign(radicalEncoding.size(), {});
        std::map<uint16_t, int> radicalMap;
        for (int i = 0; i < radicalEncoding.size(); i++) {
            radicalMap.clear();
            std::array<uint32_t, maxLength> mask{};
            for (int j = 0; j < maxLength; j++) {
                for (int k = 0; k < maxLength; k++)
                    mask[k] *= (numKeys + 1);
                if (radicalEncoding[i][j] != 0) {
                    if (!radicalMap.contains(radicalEncoding[i][j])) {
                        radicalToHanzi[radicalEncoding[i][j]].push_back(std::make_pair(i, j));
                        radicalMap[radicalEncoding[i][j]] = j;
                    }
                    mask[radicalMap[radicalEncoding[i][j]]]++;
                }
            }
            radicalMask[i] = mask;
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
            auto tmp = hashIndex(keyEncoding[i]);
            numberEncoding[i] = tmp;
            auto &v = hashTable[tmp];
            if (v > 0)
                dupCnt++;
            v++;
        }
    }
    void modifyRadical(int radicalIdx, uint16_t newKey) {
        uint16_t oldKey = radicalToKey[radicalIdx];
        radicalToKey[radicalIdx] = newKey;
        const auto &affectedHanzi = radicalToHanzi[radicalIdx];
        for (auto hanzi : affectedHanzi) {
            auto hanziIdx = hanzi.first;
            auto pos = hanzi.second;
            dupCnt -= (--hashTable[numberEncoding[hanziIdx]] > 0);
            numberEncoding[hanziIdx] =
                numberEncoding[hanziIdx] - radicalMask[hanziIdx][pos] * (static_cast<int>(oldKey) - newKey);
            dupCnt += (++hashTable[numberEncoding[hanziIdx]] > 1);
        }
    }
    int numRadical;
    int dupCnt;
    std::vector<uint16_t> radicalToKey;                           // 字根对应的字母
    std::vector<std::array<uint16_t, maxLength>> radicalEncoding; // 汉字的字根编码
    std::vector<std::array<uint16_t, maxLength>> keyEncoding;     // 汉字的字母编码
    std::vector<std::vector<std::pair<int, int>>> radicalToHanzi; // 包合某个字根的汉字
    std::vector<uint16_t> hashTable;                              // 汉字哈希桶
    std::vector<uint32_t> numberEncoding;                         // 汉字的哈希编码
    std::vector<std::array<uint32_t, maxLength>> radicalMask;     // 字根在汉字中是掩码
};
