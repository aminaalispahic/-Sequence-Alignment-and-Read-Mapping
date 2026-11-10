# Sequence Alignment i Read Mapping

Projekat iz predmeta **Paralelni računarski sistemi**, namijenjen proučavanju algoritama i njihovih performansi. <b></br>
Cilj je implementirati Smith–Waterman algoritam i jednostavan read mapper, a zatim uporediti serijsko izvršavanje s paralelnim verzijama pomoću OpenMP-a i CUDA-e.


## Cilj projekta

Program treba da pronađe gdje se kratki DNK readovi nalaze u dužoj referentnoj sekvenci. Pošto readovi mogu sadržavati zamjene, umetanja ili brisanja, tražena lokacija ne mora biti potpuno identična readu.

Smith–Waterman pronalazi najbolje lokalno poravnanje između reada i kandidatske regije genoma. Serijska verzija služi kao osnovna verzija za provjeru rezultata i poređenje performansi.

## Kako radi

### Smith–Waterman

Smith–Waterman je algoritam za lokalno poravnanje dvije sekvence. Za svaku kombinaciju pozicija računa se najbolji skor u matrici dinamičkog programiranja. Poklapanja dobijaju pozitivan broj bodova, a nepoklapanja i praznine se kažnjavaju.

Nakon popunjavanja matrice, najveći skor označava najbolje lokalno poravnanje. Traceback od te ćelije do nule određuje koje baze su poravnate i omogućava ispis poravnanja i CIGAR zapisa. Linearna kazna za prazninu obračunava se po znaku praznine.

### Read mapping

Mapper koristi strategiju **seed-and-extend**:

1. **Indeksiranje:** referentni genom se dijeli na k-mere, kratke podsekvence fiksne dužine. Planirana početna vrijednost je **k = 12–15**. Indeks je sortirani niz parova (k-mer, pozicija); prečesti k-meri mogu se preskočiti.
2. **Seeding:** k-meri iz reada traže se u indeksu. Njihove pozicije služe za glasanje o mogućim lokacijama reada u genomu.
3. **Extension:** Smith–Waterman detaljno poravna read s manjim regijama oko odabranih kandidata.
4. **Rezultat:** mapper bira najbolju kandidatsku lokaciju i vraća poziciju, skor i CIGAR zapis. Read bez odgovarajućeg kandidata označava se kao nemapiran.

Ovim pristupom se Smith–Waterman ne pokreće nad cijelim genomom za svaki read, što bi bilo presporo.

## Verzije i paralelizacija

- **Serijska verzija:** obrađuje readove jedan po jedan i služi kao referenca za poređenje.
- **OpenMP verzija:** raspoređuje nezavisne readove na CPU niti. Verzije se pokreću na istim readovima i s istim parametrima; provjerava se da vraćaju iste skorove i pozicije.
- **CUDA verzija:** obrađuje grupe nezavisnih parova (read, kandidatska regija) na GPU-u. Rezultati se porede sa serijskom verzijom. Veličina grupe i vrijeme prenosa podataka bilježe se tokom mjerenja.

## Tehnologije

| Dio | Tehnologija |
| --- | --- |
| Jezik | C++17 |
| CPU paralelizacija | OpenMP |
| GPU paralelizacija | CUDA |
| Provjera poravnanja | Biopython PairwiseAligner ili EMBOSS water |
| Testni podaci | Vlastiti generator ili wgsim |
| Mjerenje vremena | std::chrono |
| Profiliranje | perf ili gprof |

## Testni podaci i provjera

Planirani skupovi podataka:

- **Mali primjer:** TGTTACGG i GGTTGACTA za ručnu provjeru poravnanja.
- **Sintetički podaci:** nasumična referentna sekvenca od približno 1–5 miliona baza i readovi dužine 100–150 baza s kontrolisanim zamjenama i indelima. Poznate početne pozicije omogućavaju provjeru tačnosti mapiranja.
- **Stvarni podaci:** bakterijski genom, na primjer E. coli, i simulirani readovi za realističniju evaluaciju.

Serijska implementacija provjerava se na malim, ručno pregledanim primjerima i poređenjem s Biopython PairwiseAlignerom ili EMBOSS waterom. Parametri bodovanja moraju biti isti u implementaciji i referentnom alatu. Rezultati OpenMP-a i CUDA-e provjeravaju se prema serijskoj verziji.

## Mjerenje performansi

Sva poređenja rade se na istom skupu readova i uz iste parametre. Mjerenja se ponavljaju više puta, a prikazuje se medijana.

Za OpenMP se mijenja broj CPU niti p:

- **Ubrzanje:** S_OpenMP(p) = T_serijski / T_OpenMP(p)
- **Efikasnost:** E_OpenMP(p) = S_OpenMP(p) / p

Za CUDA se mijenja veličina grupe b:

- **Ubrzanje:** S_CUDA(b) = T_serijski / T_CUDA(b)

Vrijeme prenosa podataka na GPU bilježi se posebno ili se uključuje u ukupno CUDA vrijeme; izvještaj treba navesti koja je mjera korištena.

**GCUPS** predstavlja broj stvarno izračunatih ćelija Smith–Waterman matrica podijeljen vremenom poravnanja i sa 10^9. Tačnost mapiranja je udio readova pronađenih na pravoj poziciji, uz unaprijed određenu toleranciju od nekoliko baza.

## Status i pokretanje

Ovaj README opisuje plan projekta. Komande za kompajliranje i pokretanje biće dodane kada se utvrde struktura izvornog koda, ulazni argumenti programa i dostupna CUDA konfiguracija.
