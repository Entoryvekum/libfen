#pragma once
#include "hash_table6.hpp"
#include <array>
#include <chrono>
#include <cstdint>
#include <format>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>

// mappingData: 字根数据、字根-按键映射数据、字根-汉字反查数据。
// encodingData: 汉字编码数据，维护重码数。若编码某一位为0，则表示该位待定。

constexpr uint16_t isRadical = 1 << 15;

template <int max_length> class mappingData
{
public:
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
            std::array<uint16_t, max_length> radicalRow{};
            int value;
            for (int count = 0; iss >> value && count < max_length; count++) {
                radicalRow[count] = static_cast<uint16_t>(value);
                if (value > numRadical) {
                    numRadical = value;
                }
            }
            radicalEncoding.push_back(radicalRow);
        }
        radicalToKey.assign(numRadical + 1, {});
        radicalToHanzi.assign(numRadical + 1, {});
        radicalWeight.assign(numRadical + 1, {});
        std::set<uint16_t> radicalSet;
        for (int i = 0; i < radicalEncoding.size();i++) {
            radicalSet.clear();
            for (int j = 0; j < max_length; j++) {
                if (radicalEncoding[i][j] != 0) {
                    if(radicalSet.find(radicalEncoding[i][j])==radicalSet.end()) {
                        radicalToHanzi[radicalEncoding[i][j]].push_back(i);
                        radicalSet.insert(radicalEncoding[i][j]);
                    }
                    radicalWeight[radicalEncoding[i][j]]++;
                }
                else
                    break;
            }
        }
    }
    void writeMapping(const std::string &directory) {
        auto now = std::chrono::system_clock::now();
        std::string filename = std::format("{}/{:%y-%m-%d-%H-%M-%S-mapping}.txt", directory, now);

        std::ofstream file(filename);
        if (!file.is_open())
            throw std::runtime_error("Cannot create file: " + filename);

        for (int i=1;i< radicalToKey.size();i++) {
            file <<radicalToKey[i] << "\n";
        }
    }
    void readMapping(const std::string &filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        int i=1;
        for(;i<=numRadical&&(file>>radicalToKey[i]);i++);
        if(i-1!=numRadical)
            throw std::runtime_error("Read mapping error: " + filename);
    }
    int numRadical;
    std::vector<uint16_t> radicalToKey;
    std::vector<double> radicalWeight;
    std::vector<std::array<uint16_t, max_length>> radicalEncoding;
    std::vector<std::vector<int>> radicalToHanzi;
};

template <int len> uint64_t index(const std::array<uint16_t, len> &radical, const std::array<uint16_t, len> &encoding) {
    uint64_t rst = 0;
    for (int i = 0; i < len; i++) {
        uint64_t v = 0;
        if (encoding[i] != 0)
            v = encoding[i];
        else
            v = static_cast<uint64_t>(radical[i] | isRadical);
        rst |= (v << (i * 16));
    }
    return rst;
}

template <int max_length> class encodingData
{
public:
    encodingData(mappingData<max_length> &mp)
        : mapping(mp), keyEncoding(mp.radicalEncoding.size()), numRadical(mp.numRadical) {}
    void buildHash() {
        dupCnt=0;
        hash.clear();
        hash.reserve(keyEncoding.size());
        for (int i = 0; i < keyEncoding.size(); i++) {
            auto &v = hash[index<max_length>(mapping.radicalEncoding[i], keyEncoding[i])];
            if (v > 0)
                dupCnt++;
            v++;
        }
    }
    void modify(int pos, const std::array<uint16_t, max_length> &aft) {
        const auto srcIndex = index<max_length>(mapping.radicalEncoding[pos], keyEncoding[pos]);
        auto &src = hash.at(srcIndex);
        if (src == 0)
            throw std::runtime_error("Bucket not found");
        if (src == 1)
            hash.erase(srcIndex);
        else {
            src--;
            dupCnt--;
        }
        const auto trgIndex = index<max_length>(mapping.radicalEncoding[pos], aft);
        auto &trg = hash[trgIndex];
        if (trg == 0)
            trg = 1;
        else {
            trg++;
            dupCnt++;
        }
    }
    void createEncodingFromMapping() {
        for (int i = 0; i < mapping.radicalEncoding.size(); i++) {
            for (int j = 0; j < max_length; j++) {
                keyEncoding[i][j] = mapping.radicalToKey[mapping.radicalEncoding[i][j]];
            }
        }
    }
    void writeKeyEncoding(const std::string &directory) {
        auto now = std::chrono::system_clock::now();
        std::string filename = std::format("{}/{:%y-%m-%d-%H-%M-%S-encoding}.txt", directory, now);

        std::ofstream file(filename);
        if (!file.is_open())
            throw std::runtime_error("Cannot create file: " + filename);

        file<<dupCnt<<std::endl;
        for (const auto &row : keyEncoding) {
            for (int i = 0; i < max_length; ++i) {
                file << row[i] << (i == max_length - 1 ? "" : " ");
            }
            file << "\n";
        }
    }
    void modifyMapping(int radicalIdx, int newKey) {
        mapping.radicalToKey[radicalIdx] = newKey;
        const auto &affectedHanzi = mapping.radicalToHanzi[radicalIdx];
        for (int hanziIdx : affectedHanzi) {
            std::array<uint16_t, max_length> nextEncoding = keyEncoding[hanziIdx];
            for (int j = 0; j < max_length; ++j) {
                if (mapping.radicalEncoding[hanziIdx][j] == radicalIdx)
                    nextEncoding[j] = newKey;
                else
                    nextEncoding[j] = keyEncoding[hanziIdx][j];
            }
            modify(hanziIdx, nextEncoding);
            keyEncoding[hanziIdx] = nextEncoding;
        }
    }
    int dupCnt;
    const int numRadical;
    mappingData<max_length> &mapping;
    std::vector<std::array<uint16_t, max_length>> keyEncoding;
    emhash6::HashMap<uint64_t, uint16_t> hash;
};