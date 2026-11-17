$ErrorActionPreference = 'Stop'

$originalPath = $env:Path
$compilerDir = 'C:\msys64\ucrt64\bin'
if (Test-Path -LiteralPath (Join-Path $compilerDir 'g++.exe')) {
    $env:Path = "$compilerDir;$env:Path"
}

$scratch = Join-Path ([IO.Path]::GetTempPath()) ('prs-sw-bafers-' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $scratch | Out-Null
$probe = Join-Path $scratch 'probe.cpp'
$header = Join-Path $scratch 'align.hpp'
$exe = Join-Path $scratch 'probe.exe'
$mapperHeader = Join-Path $scratch 'mapper.hpp'
$fastaHeader = Join-Path $scratch 'fasta.hpp'
$basesHeader = Join-Path $scratch 'bases.hpp'

try {
    @'
#include <cstdio>
#include <random>
#include <type_traits>
#include "fasta.hpp"
#include "mapper.hpp"

bool checkAlignment() {
    std::mt19937 rng(7);
    std::string ref(250, 'A');
    for (char& base : ref) base = "ACGT"[rng() & 3];
    const std::string read = ref.substr(40, 150);
    sw::KmerIndex index;
    sw::MapParams params;
    sw::Stats stats;
    index.build(ref, params.k);
    const auto result = sw::mapRead(read, index, ref, params, stats);
    return stats.extensions == 1 && result.pos == 40 && result.score == 450;
}

void printInputLayout() {
    std::printf("\nBig ulaz nakon pripreme\n");
#if PRS_PROBE_VERSION == 3
    std::string ref;
    {
        const auto refs = sw::readSequences("data/big_ref.fa");
        ref = sw::encodeSequence(refs.at(0).seq);
    }
    const auto reads = sw::makeBatch(sw::readSequences("data/big_reads.fq"));
    size_t idBytes = 0;
    for (const auto& id : reads.ids) idBytes += id.size();
    std::printf("Readovi: %zu | referenca: %zu B\n", reads.size(), ref.size());
    std::printf("Baze readova: %zu B | FASTQ kvalitete: 0 B\n", reads.bases.size());
    std::printf("Reverzni komplementi: %zu B | offseti: %zu elemenata, %zu B\n",
                reads.reverseBases.size(), reads.offsets.size(),
                reads.offsets.size() * sizeof(reads.offsets[0]));
    std::printf("Tekst ID-eva: %zu B | objekti za ID-eve: %zu B\n",
                idBytes, reads.ids.size() * sizeof(reads.ids[0]));
#else
    const auto refs = sw::readSequences("data/big_ref.fa");
    const auto reads = sw::readSequences("data/big_reads.fq");
    size_t bases = 0, qualities = 0, idBytes = 0;
    for (const auto& record : reads) {
        bases += record.seq.size();
        qualities += record.qual.size();
        idBytes += record.id.size();
    }
    std::printf("Readovi: %zu | referenca: %zu B\n", reads.size(), refs.at(0).seq.size());
    std::printf("Baze readova: %zu B | FASTQ kvalitete: %zu B\n", bases, qualities);
    std::printf("Reverzni komplementi u skupu readova: 0 B | offseti: 0 B\n");
    std::printf("Tekst ID-eva: %zu B | objekti Record: %zu B\n",
                idBytes, reads.size() * sizeof(reads[0]));
#endif
}

int main() {
    if (!checkAlignment()) return 1;
    printInputLayout();
    return 0;
}
'@ | Set-Content -LiteralPath $probe -Encoding ascii

    $fullAnchor = 'std::vector<uint8_t> T((n + 1) * W, STOP);'
    $fullSizes = @'
    const size_t scoreCount = H.size(), directionCount = T.size();
    const size_t scoreBytes = scoreCount * sizeof(H[0]);
    const size_t directionBytes = directionCount * sizeof(T[0]);
'@
    $rollingAnchor = 'directions.resize(((size_t)n + 1) * width);'
    $rollingSizes = @'
    const size_t scoreCount = previous.size() + current.size();
    const size_t directionCount = directions.size();
    const size_t scoreBytes = scoreCount * sizeof(previous[0]);
    const size_t directionBytes = directionCount * sizeof(directions[0]);
'@
    $printSizes = @'
    const size_t regionCopyBytes =
        std::is_same_v<std::decay_t<decltype(Y)>, std::string> ? Y.size() : 0;
    std::printf("Referentni prozor: %zu baza | kopija znakova: %zu B\n",
                Y.size(), regionCopyBytes);
    std::printf("Skorovi: %zu elemenata, %zu B\n", scoreCount, scoreBytes);
    std::printf("Smjerovi: %zu elemenata, %zu B\n", directionCount, directionBytes);
    std::printf("Ukupno SW baferi: %zu B\n", scoreBytes + directionBytes);
'@

    Push-Location -LiteralPath $PSScriptRoot
    try {
        foreach ($number in 0..3) {
            $version = "V$number"
            $directory = Join-Path $PSScriptRoot $version
            Copy-Item -LiteralPath (Join-Path $directory 'mapper.hpp') -Destination $mapperHeader
            Copy-Item -LiteralPath (Join-Path $directory 'fasta.hpp') -Destination $fastaHeader
            if ($number -eq 3) {
                Copy-Item -LiteralPath (Join-Path $directory 'bases.hpp') -Destination $basesHeader
            }
            $code = [IO.File]::ReadAllText((Join-Path $directory 'align.hpp'))
            if ($number -lt 2) {
                $anchor = $fullAnchor
                $sizes = $fullSizes
            } else {
                $anchor = $rollingAnchor
                $sizes = $rollingSizes
            }
            if ([regex]::Matches($code, [regex]::Escape($anchor)).Count -ne 1) {
                throw "Nije pronadjena jedinstvena linija za $version."
            }
            $code = $code.Replace($anchor, "$anchor`n$sizes`n$printSizes")
            [IO.File]::WriteAllText($header, $code, [Text.UTF8Encoding]::new($false))

            & g++ -O2 -std=c++17 "-DPRS_PROBE_VERSION=$number" -o $exe $probe
            if ($LASTEXITCODE -ne 0) { throw "Kompajliranje provjere nije uspjelo za $version." }

            Write-Host "$version - read 150 baza, prozor 190 baza"
            & $exe
            if ($LASTEXITCODE -ne 0) { throw "Provjera poravnanja nije uspjela za $version." }
            Write-Host ''
        }
    } finally {
        Pop-Location
    }
} finally {
    $env:Path = $originalPath
    foreach ($file in $probe, $header, $exe, $mapperHeader, $fastaHeader, $basesHeader) {
        if (Test-Path -LiteralPath $file) { Remove-Item -LiteralPath $file }
    }
    Remove-Item -LiteralPath $scratch
}
