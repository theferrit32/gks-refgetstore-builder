# Some residues do not round-trip through `RefgetStore`'s encoded storage mode

When a sequence is added to a gtars `RefgetStore` in the default
`StorageMode.Encoded`, a small set of residues comes back from the store as a
different character. The sequence digest is computed from the input bytes and
is correct; the bytes returned by `stream_sequence` differ from the input at
those positions, so they no longer hash to that digest.

This document describes the cases observed, where they arise in the source, and
a minimal script that reproduces them.

- **gtars versions tested:** 0.10.0 (the latest release on PyPI at the time of
  writing) and 0.9.2. Both behave identically.
- **Source references:** tag `gtars-python-v0.10.0` / `gtars-refget-v0.10.0`,
  commit [`90141b68`](https://github.com/databio/gtars/tree/90141b683aaa1c5bae4d452718fc440825a5f711).
  `digest/alphabet.rs` and `digest/encoder.rs` are unchanged between the 0.9.2
  tag and the `gtars-python-v0.11.0` tag (`9777791b`).
- **Entry point used:** `RefgetStore.add_sequence_collection_from_fasta`.

## Summary

| input residue | alphabet selected | returned as | where it changes | original recoverable from stored bytes |
|---|---|---|---|---|
| `D` (A/G/T) | `DnaIupac` | `H` | decoding | yes |
| `H` (A/C/T) | `DnaIupac` | `V` | decoding | yes |
| `U` | `DnaIupac` | `T` | encoding | no |
| `U` (selenocysteine) | `Protein` | `A` | encoding | no |
| `B` (Asx) | `Protein` | `A` | encoding | no |

(`DnaIupac` is reported as `dnaio` in sequence metadata and in the output below.)

There are two underlying causes:

1. In the `DnaIupac` alphabet, the decoding table maps two codes to different
   letters than the encoding table assigned to them. The stored bits are
   distinct for `D`, `H` and `V`, so the original sequence can be recovered by
   correcting the decoding table.
2. Some residues are accepted into an alphabet that has no code of their own
   for them, so they are stored using another residue's code. This applies to
   `U` in `DnaIupac` (stored as `T`'s code) and to `U` and `B` in `Protein`
   (stored as `A`'s code).

In `StorageMode.Raw` all of these residues round-trip unchanged.

## Observed occurrences

These cases were found while building a `RefgetStore` of human reference
sequences with [gks-refgetstore-builder](https://github.com/theferrit32/gks-refgetstore-builder).
Its verification step reads every stored sequence back with `stream_sequence`,
recomputes the sha512t24u digest of the returned bytes, and compares it with
the digest the sequence is stored under.

The store holds 1,779,052 sequences from NCBI, Ensembl and EBI sources, built
with gtars 0.9.2. Reading it with 0.10.0 gives the same results. 133 sequences
(0.0075%) return bytes that differ from their input:

| alphabet | sequences | returned bytes differ from input |
|---|---:|---:|
| dna2bit | 1,426,421 | 0 |
| dna3bit | 643 | 0 |
| dnaio | 107 | 27 |
| protein | 351,881 | 106 |

Re-ingesting the same source FASTA files produces the same results each time.

The `dnaio` group is small by count but holds 1.82 Gbp: a single ambiguity code
anywhere in a sequence selects this alphabet for the whole sequence, and several
primary-assembly chromosomes (for example `NC_000001.11`) contain one. The
chromosomes in this store happen not to contain `D` or `H`, and round-trip
correctly.

## Reproducing

`repro_minimal.py`, in the same directory as this document, depends only on
gtars. It needs no existing store, downloaded data or network access:

    uv run --with gtars==0.10.0 python repro_minimal.py

It runs three stages:

1. **Detection:** calls `digest_sequence` on each test sequence and prints the
   alphabet selected and the digest.
2. **Encoded round trip:** adds the sequences to `RefgetStore.in_memory()` in
   the default `StorageMode.Encoded`, reads each one back, and reports the first
   differing position.
3. **Raw round trip:** repeats stage 2 after `store.disable_encoding()`.

It exits 0 when the results match the expectations in this document. It exits
non-zero if any result changes, including a case that now round-trips
correctly, so it can also confirm a fix.

Output with gtars 0.10.0 (abridged):

    stage 1 -- detection (pure functions, no store)
      case               len  alphabet   sha512t24u                        expect
      selenoprotein      180  protein    z_3gTL7__q3R6SR1-8NLSFz7-H1RbDdV  protein
      protein-as-dnaio    12  dnaio      13xx4yn7JousnAXLxTBX-bbWg4TitdN0  dnaio
      nucleotide-iupac    19  dnaio      Vu6wOvjbNQLTfbmxUCCcf9tBB9KZJYDu  dnaio
      pyrrolysine         19  ASCII      p8m0yov-U6dljEuwih4z-J4D_bHJo1eJ  ASCII
      plain-protein       18  protein    Yh8OEOTbVpE7N7GDVoGQS-x77N0Wyuyj  protein
      plain-dna           12  dna2bit    tlkjFUbZBvLI4IYBwYgqhjZI_Aetdq5B  dna2bit

    stage 2 -- round trip through RefgetStore (StorageMode.Encoded, the default)
      case              expect   actual   substitution
      selenoprotein     differs  differs  index 109: U -> A
      protein-as-dnaio  differs  differs  index 5: H -> V
      nucleotide-iupac  differs  differs  index 4: H -> V
      pyrrolysine       matches  matches
      plain-protein     matches  matches
      plain-dna         matches  matches

    stage 3 -- the same round trip at StorageMode.Raw (control)
      ... all six match ...

A final per-residue check places individual residues in a fixed nucleotide or
protein context and reports which ones round-trip. Every "returned as" result in
this document comes from that check.

## How sequences are stored

On import, gtars selects an alphabet from the characters in the sequence, then
packs the sequence through that alphabet's lookup table
([`store/import.rs:277-296`](https://github.com/databio/gtars/blob/90141b683aaa1c5bae4d452718fc440825a5f711/gtars-refget/src/store/import.rs#L277-L296);
the batched path at `:344-346` selects the alphabet the same way):

```rust
let mut guesser = crate::digest::AlphabetGuesser::new();
guesser.update(&raw_bytes);
let alphabet = guesser.guess();
// …
let sequence_data = match mode {
    StorageMode::Encoded => {
        let mut encoder = SequenceEncoder::new(alphabet, length);
        encoder.update(&raw_bytes);
        drop(raw_bytes);
        encoder.finalize()
    }
    StorageMode::Raw => raw_bytes,
};
```

The encoder converts each byte with a table lookup
([`digest/encoder.rs:50-64`](https://github.com/databio/gtars/blob/90141b683aaa1c5bae4d452718fc440825a5f711/gtars-refget/src/digest/encoder.rs#L50-L64)).
A byte with no entry in the table takes the table's default value, `0`:

```rust
pub fn update(&mut self, sequence: &[u8]) {
    for &byte in sequence {
        let code = self.alphabet.encoding_array[byte as usize] as u64;
        self.buffer = (self.buffer << self.alphabet.bits_per_symbol) | code;
        // … flush whole bytes
    }
}
```

In `StorageMode.Raw` the bytes are stored as given, and all of the sequences in
this document round-trip. This places the changes in the encoding tables and in
alphabet selection, rather than in FASTA parsing.

## 1. `D` and `H` in the `DnaIupac` alphabet

### Encoding and decoding tables

The encoding table gives `D`, `H` and `V` three distinct 4-bit codes
([`digest/alphabet.rs:194-230`](https://github.com/databio/gtars/blob/90141b683aaa1c5bae4d452718fc440825a5f711/gtars-refget/src/digest/alphabet.rs#L194-L230)):

```rust
const DNA_IUPAC_ENCODING_ARRAY: [u8; 256] = {
    let mut arr = [0u8; 256];
    // …
    arr[b'B' as usize] = 0b1100; // C or G or T
    arr[b'D' as usize] = 0b1101; // A or G or T
    arr[b'H' as usize] = 0b1110; // A or C or T
    arr[b'V' as usize] = 0b1111; // A or C or G
```

The decoding table maps codes `0b1101` and `0b1110` to `H` and `V`
([`digest/alphabet.rs:234-253`](https://github.com/databio/gtars/blob/90141b683aaa1c5bae4d452718fc440825a5f711/gtars-refget/src/digest/alphabet.rs#L234-L253)):

```rust
const DNA_IUPAC_DECODING_ARRAY: [u8; 256] = {
    let mut arr = [b'N'; 256]; // Default to 'N' for all
    // …
    arr[0b1011] = b'D'; // A or G or T
    arr[0b1100] = b'B'; // C or G or T
    arr[0b1101] = b'H'; // A or C or T
    arr[0b1110] = b'V'; // A or C or G
    arr[0b1111] = b'V'; // A or C or G (or use 'N' if you'd rather treat 0b1111 as invalid)
```

Taken together:

| input | encoded as | decoded as |
|---|---|---|
| `B` | `0b1100` | `B` |
| `D` | `0b1101` | `H` |
| `H` | `0b1110` | `V` |
| `V` | `0b1111` | `V` |

No input is encoded as `0b1011`, so the decoding entry `0b1011 → D` is not
reached. Because `0b1101` and `0b1110` are each produced by only one input
residue, the stored bits still identify the original residue.

Testing each IUPAC code on its own in an `ACGT` context:

| code | alphabet selected | result |
|---|---|---|
| `R` `Y` `N` | `dna3bit` | round-trips |
| `S` `W` `K` `M` `B` `V` | `dnaio` | round-trips |
| `D` (A/G/T) | `dnaio` | returned as `H` (A/C/T) |
| `H` (A/C/T) | `dnaio` | returned as `V` (A/C/G) |

### Example: a RefSeq non-coding RNA

`NR_103745.1` (PTENP1-AS, a human long non-coding RNA, 873 nt, from
`GCF_000001405.26_GRCh38_rna.fna.gz`) contains seven ambiguity codes:

```
ambiguity codes at  542 K   599 K   620 K   668 Y   744 Y   814 Y   852 D
index 852           input D  ->  returned H

sha512t24u(input bytes)     o_Juegh-7Rliguu4uJ72pColqqJQt4mz   <- stored digest
sha512t24u(returned bytes)  a7NEDaBDsfrHPxhhYH9hPovgecVJOIhj
```

The `K` and `Y` codes round-trip. The `D` at index 852 is returned as `H`.

### Example: a short protein

Alphabet selection is based on which characters occur, not on sequence length
or composition. A protein made up entirely of letters that are also IUPAC
nucleotide codes is assigned the `DnaIupac` alphabet.

`ENSP00000499040.1` (NOTCH2, 12 aa, from Ensembl release 113
`Homo_sapiens.GRCh38.pep.all.fa.gz`):

```
input     MCVTYHNGTGYC
returned  MCVTYVNGTGYC
               ^ index 5: H -> V

sha512t24u(input bytes)     13xx4yn7JousnAXLxTBX-bbWg4TitdN0   <- stored digest
sha512t24u(returned bytes)  fQh3MU3P7sOQ7EihvRecrb08e86fBG0G
```

The alphabet selection alone would not change the returned sequence, since
`DnaIupac` has a code for every letter here. The change at index 5 comes from
the `H` decoding entry.

### Occurrences in the test store

All 27 `dnaio` sequences whose returned bytes differ are explained by these two
substitutions. Each was confirmed by reversing the substitutions and checking
that the result hashes to the stored digest. 22 contain at least one `D`, 9
contain at least one `H`, and 4 contain both.

26 of the 27 are Ensembl `ENSP*` protein records of 3 to 20 residues. The
remaining one is `NR_103745.1`, the RNA record above.

### Recovering existing data

Because the stored bits are distinct, correcting the two decoding entries
(`0b1101 → D`, `0b1110 → H`) recovers the original sequences without
re-importing them. `check_decoder_fix.py` decodes every `dnaio` payload in a
store twice: once with the current table, which reproduces gtars' output
exactly, and once with the two entries changed:

| round-trips with current table | round-trips with corrected table | sequences |
|---|---|---:|
| no | yes | 27 |
| yes | yes | 80 |

All 107 `dnaio` sequences (1.82 Gbp, including the chromosomes) round-trip with
the corrected table, and none of the 80 that round-trip today changes. This is
also expected from the tables: a sequence that round-trips today contains
neither of the two affected codes.

## 2. `U` in the `DnaIupac` alphabet

The `DnaIupac` encoding table gives `U` the same code as `T`
([`digest/alphabet.rs:199-200`](https://github.com/databio/gtars/blob/90141b683aaa1c5bae4d452718fc440825a5f711/gtars-refget/src/digest/alphabet.rs#L199-L200)):

```rust
    arr[b'T' as usize] = 0b1000; // T
    arr[b'U' as usize] = 0b1000; // U (common in RNA)
```

`U` is accepted by the `DnaIupac` alphabet, so an RNA sequence written with `U`
is assigned `DnaIupac` even when it has no ambiguity codes, and is returned
with `T` in place of `U`:

```
input           returned
ACGUACGUACGU    ACGTACGTACGT
MCVTYUNGTGYC    MCVTYTNGTGYC
```

The digest is computed from the input, so it identifies the `U` form while the
store returns the `T` form. Since both letters share one code, the original
cannot be recovered from the stored bytes. The human sources used for the test
store write RNA with `T`, so it contains no instances; the reproducer's
per-residue check includes this case.

## 3. `U` and `B` in the `Protein` alphabet

Selenocysteine is written `U` in the IUPAC one-letter amino-acid code, and
appears in both Ensembl and RefSeq protein files. A sequence containing `U` is
assigned the `Protein` alphabet, and each `U` is returned as `A`.

`ENSP00000473614` (GPX4, glutathione peroxidase 4, 180 aa, from Ensembl release
76 `Homo_sapiens.GRCh38.pep.all.fa.gz`):

```
input charset     ACDEFGHIKLMNPQRSTUVWY
returned charset  ACDEFGHIKLMNPQRSTVWY
index 109         input U  ->  returned A

sha512t24u(input bytes)     z_3gTL7__q3R6SR1-8NLSFz7-H1RbDdV   <- stored digest
sha512t24u(returned bytes)  4otSS66T9mEsIeu_jmydFFDB_x-LFqsb
```

```python
from gtars.refget import RefgetStore
store = RefgetStore.open_local("store")
store.stream_sequence("z_3gTL7__q3R6SR1-8NLSFz7-H1RbDdV").read_all()[109]
# 'A'; the input FASTA has 'U'
```

`B` (Asx, asparagine or aspartic acid) behaves the same way.

### Why `U` and `B` are returned as `A`

The protein encoding table is initialized to zero, and `0` is also the code for
alanine
([`digest/alphabet.rs:255-259`](https://github.com/databio/gtars/blob/90141b683aaa1c5bae4d452718fc440825a5f711/gtars-refget/src/digest/alphabet.rs#L255-L259)):

```rust
const PROTEIN_ENCODING_ARRAY: [u8; 256] = {
    let mut arr = [0u8; 256];
    // Standard amino acids (20) + special characters
    arr[b'A' as usize] = 0b00000;
    arr[b'a' as usize] = 0b00000; // Alanine
```

The table assigns codes `0b00000` through `0b10111` (the 20 standard residues,
`*`, `X`, `-` and `.`), leaving 8 of the 32 five-bit codes unused. `U` and `B`
have no entry, so they are encoded as `0` and decoded as `A`. The original
residue cannot be recovered from the stored bytes.

### Why the `Protein` alphabet is selected

Alphabet selection checks each character against the tables in a fixed order
and records the first alphabet that accepts it; the IUPAC table is checked
before the protein table
([`digest/alphabet.rs:488-511`](https://github.com/databio/gtars/blob/90141b683aaa1c5bae4d452718fc440825a5f711/gtars-refget/src/digest/alphabet.rs#L488-L511)):

```rust
fn get_minimum_alphabet_for_char(byte: u8) -> AlphabetType {
    // Check 2-bit DNA first (most restrictive)
    if matches!(byte, b'A' | b'C' | b'G' | b'T') {
        return AlphabetType::Dna2bit;
    }

    // Check 3-bit DNA
    if matches!(byte, b'N' | b'R' | b'Y') {
        return AlphabetType::Dna3bit;
    }

    // Check IUPAC DNA
    if DNA_IUPAC_ENCODING_ARRAY[byte as usize] != 0 || byte == b'N' {
        return AlphabetType::DnaIupac;
    }

    // Check Protein
    if PROTEIN_ENCODING_ARRAY[byte as usize] != 0 || byte == b'-' || byte == b'*' {
        return AlphabetType::Protein;
    }

    // Default to ASCII
    AlphabetType::Ascii
}
```

The sequence-level alphabet is then the most general of the per-character
results, using a fixed ordering
([`digest/alphabet.rs:514-527`](https://github.com/databio/gtars/blob/90141b683aaa1c5bae4d452718fc440825a5f711/gtars-refget/src/digest/alphabet.rs#L514-L527)):

```rust
fn is_more_general_alphabet(a: AlphabetType, b: AlphabetType) -> bool {
    let alphabet_hierarchy = [
        AlphabetType::Dna2bit,
        AlphabetType::Dna3bit,
        AlphabetType::DnaIupac,
        AlphabetType::Protein,
        AlphabetType::Ascii,
    ];
    // …
    pos_a > pos_b
}
```

This ordering treats each alphabet as able to encode every character that the
alphabets before it can. That holds for most characters, but `U` and `B` are
accepted by `DnaIupac` and have no code in `Protein`. When a sequence already at
`Protein` contains `U`, the per-character result is `DnaIupac`, which ranks
below `Protein`, so the sequence stays at `Protein`. Every other IUPAC
nucleotide letter is also an amino-acid code, so `U` and `B` are the only
characters affected.

For comparison, `J`, `O` (pyrrolysine) and `Z` (Glx) appear in neither table,
so they select the `ASCII` alphabet and round-trip byte for byte. `X` has a
protein code and round-trips.

### Occurrences in the test store

All 106 `protein` sequences whose returned bytes differ contain selenocysteine.
Each was confirmed by substituting `U` back in at `A` positions until the
result hashed to the stored digest:

| `U` residues restored | sequences |
|---|---:|
| 1 | 92 |
| 2 | 8 |
| 3 | 3 |
| 4 | 1 |
| 10 (selenoprotein P) | 2 |

54 are Ensembl `ENSP*` records and 52 are RefSeq `NP_*` records, 41 to 698
residues long. No `B` occurs in these sources; the `B` behaviour was confirmed
with the reproducer only.

## Test coverage

The existing round-trip tests use inputs that avoid the affected residues. The
`DnaIupac` encoding test uses `ACGTRYMK`
([`digest/encoder.rs:510-528`](https://github.com/databio/gtars/blob/90141b683aaa1c5bae4d452718fc440825a5f711/gtars-refget/src/digest/encoder.rs#L510-L528)),
and the protein test (`:531-542`) uses the 20 standard residues plus `*X-`:

```rust
fn test_dna_iupac_encoding() {
    let sequence = b"ACGTRYMK";
```

The streaming-decoder tests include `D` and `H` in their fixture
(`ACGTRYMKSWBDHVN-`, `:246`). Those tests compare the streaming decoder's output
with `decode_substring_from_bytes` rather than with the original input, so both
decoders produce the same output and the comparison passes
([`digest/streaming_decoder.rs:209-238`](https://github.com/databio/gtars/blob/90141b683aaa1c5bae4d452718fc440825a5f711/gtars-refget/src/digest/streaming_decoder.rs#L209-L238)):

```rust
let encoded = encode_sequence(sequence, alphabet);
// …
let expected = decode_substring_from_bytes(&encoded, start, end, alphabet);
// …
assert_eq!(out, expected, /* … */);
```

## Impact

- The returned sequence has the expected length and consists of valid residue
  codes, so the difference is not visible without comparing against the input.
- Store metadata, indexes and digests are all internally consistent. The
  difference shows up only when the returned bytes are re-hashed and compared
  with the stored digest.
- Consumers that look up a sequence by digest receive a sequence that hashes to
  a different digest.
- For `D` and `H`, the stored data is intact and a decoder change restores it.
  For `U` and `B`, the affected sequences would need to be imported again after
  an encoder change.

## Possible fixes

Listed roughly in order of impact.

1. **Correct the two `DnaIupac` decoding entries** to `0b1101 → D` and
   `0b1110 → H`. This restores all existing `DnaIupac` payloads without
   re-importing (27 of 27 recovered in the test store, 80 of 80 unchanged).
   Changing the encoding side instead would change how existing payloads
   decode, so the decoding side is the lower-risk place to correct.
2. **Compare round trips against the original input** in the tests. An
   optional check at encode time, falling back to a lossless alphabet or
   returning an error when a sequence does not round-trip, would also prevent
   similar cases in future.
3. **Give `U` its own `DnaIupac` code**, or select `ASCII` for sequences
   containing `U`. No input currently encodes to `0b1011`, so assigning
   `U → 0b1011` would not affect existing payloads.
4. **Add `U` and `B` to the `Protein` alphabet** (8 five-bit codes are unused),
   or have alphabet selection check each alphabet's own table so these residues
   select `ASCII`, as `J`, `O` and `Z` already do. Existing affected sequences
   would need to be imported again.
5. **Take sequence length or composition into account** when choosing a
   nucleotide alphabet, for example by requiring a minimum length or checking
   for residues that are not nucleotide codes. This only changes which
   sequences reach `DnaIupac`; with (1) in place it has no effect on
   round-tripping.

## Supporting files

These files are in the same directory as this document.

- `repro_minimal.py` — the reproducer described above. Depends only on gtars.
- `check_decoder_fix.py` — decodes every `dnaio` payload in a store with the
  current and the corrected decoding table, and reports which round-trip.
  Requires [gks-refgetstore-builder](https://github.com/theferrit32/gks-refgetstore-builder)
  and a built store.
- `find_examples.py` — lists stored sequences whose returned bytes differ from
  their digest, stopping after a set number per alphabet. Requires
  gks-refgetstore-builder and a built store. `examples.tsv` is a sample of its
  output from the test store.

The complete list of the 133 sequences is
[`seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv`](https://github.com/theferrit32/gks-refgetstore-builder/blob/main/seqrepo_equivalence/known_divergence/gtars_encoding_roundtrip.tsv)
in the gks-refgetstore-builder repository.
