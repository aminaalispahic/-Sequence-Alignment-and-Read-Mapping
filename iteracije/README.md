# Iteracije Smith–Waterman mappera

Direktoriji V0–V8 čuvaju izvorni kod po iteracijama. V0 je referentni mapper s opcijom profiliranja, a svaki sljedeći direktorij dodaje izmjenu opisanu u pripadajućem diffu. Za pregled izmjerenih vremena i provjeru memorije pogledati sažeti izvještaj.

## Ulazni podaci

Lokalno napraviti direktorij data unutar iteracije i u njega kopirati zajedničke fajlove big_ref.fa i big_reads.fq koji su korišteni u mjerenju. Radni direktorij u VTuneu je direktorij aplikacije, npr. iteracije/V1.

## Kompajliranje na Windowsu

Otvoriti PowerShell u direktoriju verzije. Primjer za V1:

```powershell
$env:Path = 'C:\msys64\ucrt64\bin;' + $env:Path
g++ -O2 -fno-reorder-functions -gcodeview -std=c++17 -o mapper_V1_pdb.exe mapper.cpp "-Wl,--pdb=mapper_V1_pdb.pdb"
```

Koristiti isti skup zastavica za sve verzije u jednom poređenju. Za drugu verziju promijeniti naziv izvršnog fajla i PDB-a. Oba fajla ostaviti zajedno. Generator, demo i test poravnanja kompajliraju se zasebno jer imaju vlastiti main.

## Profiliranje u Intel VTuneu

U VTuneu odabrati Microarchitecture Exploration, Detailed i interval uzorkovanja 1 ms. Po verzijama koristiti sljedeće argumente:

- V0–V5: --ref ../data/big_ref.fa --reads ../data/big_reads.fq --profile 1
- V6: --ref ../data/big_ref.fa --reads ../data/big_reads.fq --profile 1 --band 12
- V7–V8: --ref ../data/big_ref.fa --reads ../data/big_reads.fq --profile 1 --band 12 --rel-votes 0.5

Za izlaz u TSV koristiti --profile 0 i --out s putanjom izlaza. PowerShell skriptu provjeri_sw_bafere.ps1 pokretati iz direktorija iteracije; koristi privremene fajlove i lokalne big ulaze.

## Izvorne izmjene i rezultati

original-V0-profile.diff prikazuje uvođenje profilnog režima u početnu verziju. Ostali diff fajlovi porede susjedne verzije. Kompajlirani fajlovi, sirovi VTune rezultati i lokalni ulazi čuvaju se izvan Git historije.
