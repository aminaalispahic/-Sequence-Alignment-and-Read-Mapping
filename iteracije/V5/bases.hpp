#pragma once
#include <cstdint>
#include <string>
#include <string_view>

namespace sw {

using Base = uint8_t;
constexpr Base BASE_N = 'N';

// Zamijeni A/C/G/T s bajtovima 0..3. Ostali znakovi zadržavaju identitet.
inline char remapBaseByte(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return c;
    }
}

inline Base encodeBase(char c) { return static_cast<Base>(remapBaseByte(c)); }
inline char decodeBase(Base b) { return remapBaseByte(static_cast<char>(b)); }
inline char decodeBase(char c) { return c; }
inline bool isCanonical(Base b) { return b < 4; }

struct EncodedView {
    std::string_view bytes;

    size_t size() const { return bytes.size(); }
    Base operator[](size_t i) const { return static_cast<Base>(bytes[i]); }
    EncodedView substr(size_t start, size_t length) const {
        return {bytes.substr(start, length)};
    }
};

inline void appendEncoded(std::string_view sequence, std::string& out) {
    for (char c : sequence) out.push_back(static_cast<char>(encodeBase(c)));
}

inline std::string encodeSequence(std::string_view sequence) {
    std::string out;
    out.reserve(sequence.size());
    appendEncoded(sequence, out);
    return out;
}

inline void appendReverseComplement(EncodedView sequence, std::string& out) {
    for (size_t i = sequence.size(); i > 0; --i) {
        const Base b = sequence[i - 1];
        out.push_back(static_cast<char>(isCanonical(b) ? 3 - b : b));
    }
}

}  // namespace sw
