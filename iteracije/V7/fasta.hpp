// fasta.hpp: čitanje FASTA/FASTQ (automatski prepoznaje format) i pisanje
#pragma once
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sw {

struct Record { std::string id, seq, qual; };

inline std::string firstToken(const std::string& s) {
    size_t e = s.find_first_of(" \t");
    return e == std::string::npos ? s : s.substr(0, e);
}
inline void rstrip(std::string& s) { while (!s.empty() && (s.back() == '\r' || s.back() == '\n')) s.pop_back(); }
inline void upper(std::string& s) { for (char& c : s) c = (char)std::toupper((unsigned char)c); }

inline std::vector<Record> readSequences(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Ne mogu otvoriti datoteku: " + path);
    std::vector<Record> v; std::string line;
    int c = f.peek();
    if (c == '>') {
        while (std::getline(f, line)) {
            rstrip(line);
            if (line.empty()) continue;
            if (line[0] == '>') v.push_back({firstToken(line.substr(1)), "", ""});
            else if (!v.empty()) { upper(line); v.back().seq += line; }
        }
    } else if (c == '@') {
        std::string h, s, p, q;
        while (std::getline(f, h)) {
            rstrip(h); if (h.empty()) continue;
            std::getline(f, s); std::getline(f, p); std::getline(f, q);
            rstrip(s); rstrip(q); upper(s);
            v.push_back({firstToken(h.substr(1)), s, q});
        }
    } else throw std::runtime_error("Nepoznat format (očekujem FASTA ili FASTQ): " + path);
    return v;
}

inline void writeFasta(const std::string& path, const std::vector<Record>& v, size_t wrap = 80) {
    std::ofstream o(path);
    for (const auto& r : v) {
        o << '>' << r.id << '\n';
        for (size_t i = 0; i < r.seq.size(); i += wrap) o << r.seq.substr(i, wrap) << '\n';
    }
}
inline void writeFastq(const std::string& path, const std::vector<Record>& v) {
    std::ofstream o(path);
    for (const auto& r : v)
        o << '@' << r.id << '\n' << r.seq << "\n+\n"
          << (r.qual.empty() ? std::string(r.seq.size(), 'I') : r.qual) << '\n';
}

}  // namespace sw
