# Chimeric alignment tags in subread / subjunc BAM output

## Scope

When `subjunc` (or `subread-align` with fusion / long-deletion detection enabled)
reports a read that is split across more than one reference location, the
**primary** alignment is written in the standard SAM columns
(`RNAME`, `POS`, `CIGAR`, `SEQ`, `FLAG`, ...). Every **additional** section of
the read is described by a set of four auxiliary tags:

| Tag    | Type | Meaning |
|--------|------|---------|
| `CG:Z` | Z    | CIGAR of the chimeric section, written in **read orientation** |
| `CP:i` | i    | 1-based reference coordinate of the 5' end of the chimeric section |
| `CT:Z` | Z    | Strand of the chimeric section (`+` or `-`) |
| `CC:Z` | Z    | Chromosome of the chimeric section |

A read with *n* chimeric sections carries *n* copies of each tag, in the same
order. The primary section is **not** repeated in `CG:Z`; it is the `CIGAR`
column.

---

## `CG:Z` — chimeric CIGAR

`CG:Z` is a CIGAR string describing the chimeric section **in read orientation**,
i.e. 5' → 3' of the `SEQ` field. It uses the same operators as the main CIGAR:

* `M` — bases of the read that align to the chimeric chromosome.
* `S` — bases of the read that belong to *other* sections (soft-clipped here).
* `I`, `D`, `N` — as usual.

The `M` block is the part of the read that is placed on `CC` at `CP`.

Examples:

| `CG:Z` | Meaning |
|--------|---------|
| `64M86S` | first 64 read bases align to the chimeric chromosome; last 86 belong elsewhere |
| `77S73M` | first 77 read bases belong elsewhere; last 73 align to the chimeric chromosome |
| `24M126S` | first 24 read bases align to the chimeric chromosome |
| `138S12M` | last 12 read bases align to the chimeric chromosome |

---

## `CP:i` — chimeric position

`CP:i` is the **1-based** reference coordinate of the **5' end of the chimeric
section**, expressed in the orientation of the chimeric alignment:

* If the chimeric `M` block aligns to the **forward** strand, `CP:i` is the
  **leftmost** coordinate of the block.
* If the chimeric `M` block aligns to the **reverse** strand, `CP:i` is the
  **rightmost** coordinate of the block.

In other words, `CP:i` is the reference base that the **first base of the
chimeric read segment** (the 5' base of that segment) is aligned to.

### Formula in the source

`src/core.c`, around line 1486:

```c
SUBreadSprintf(r->additional_information + strlen(r->additional_information),
    CORE_ADDITIONAL_INFO_LENGTH + 1 - strlen(r->additional_information),
    "\tCG:Z:%s\tCP:i:%u\tCT:Z:%c\tCC:Z:%s",
    r->out_cigars[xk1],
    max(1, chimeric_pos + soft_clipping_movement + 1),
    strand_xor ? '-' : '+',
    chimaric_chr);
```

where

```c
soft_clipping_movement = get_soft_clipping_length(r->out_cigars[xk1]);
```

and `get_soft_clipping_length()` returns the **leading** soft-clip length of the
chimeric CIGAR (0 if the CIGAR starts with `M`). So:

```
CP:i = max(1, chimeric_pos + leading_soft_clip(CG:Z) + 1)
```

`chimeric_pos` is the 0-based reference position of the read head in the
chimeric alignment; adding the leading soft-clip and 1 converts it to the
1-based coordinate of the first aligned base of the chimeric `M` block.

---

## `CT:Z` — chimeric strand

`CT:Z` is the strand of the chimeric section, `+` or `-`.

It is reported in the **orientation of the original sequenced read**, not
necessarily the reference forward strand. Concretely:

```
CT = reference_strand_of_chimeric_M_block  XOR  (read is reverse-strand)
```

So for a **forward-strand read** (`FLAG & 0x10 == 0`), `CT` equals the reference
strand of the chimeric block. For a **reverse-strand read** (`FLAG & 0x10 != 0`),
`CT` is the opposite of the reference strand of the chimeric block.

This is why a chimeric block that matches the reference forward strand can still
carry `CT:Z:-` when the read itself is reverse-strand (see Example 2).

---

## `CC:Z` — chimeric chromosome

`CC:Z` is the name of the chromosome (or contig) that the chimeric `M` block
aligns to. It may be the same as `RNAME` (intra-chromosomal chimeric alignment)
or different (inter-chromosomal fusion).

---

## Worked examples

All examples are from the toy fusion dataset
(`UncommonCases/subjuncAlljunc/ToyExample`, reference `reference.fa`,
4 chromosomes of 3,000,000 bp). Coordinates below are 1-based.

### Example 1 — different chromosome, different strand

```
toy000006  FLAG=73  RNAME=chr4  POS=870784  CIGAR=64S86M
           CG:Z:64M86S  CP:i:1392337  CT:Z:-  CC:Z:chr2
```

* Primary: `read[64:150]` matches `chr4:870784-870869` on the **forward** strand.
* Chimeric: `read[0:64]` = `CGCACGTTTTCGAACGGAATATCTCTTGATACACATCGCATATCTCCGTTACTGTGAGCTTACC`
  whose reverse complement is
  `GGTAAGCTCACAGTAACGGAGATATGCGATGTGTATCAAGAGATATTCCGTTCGAAAACGTGCG`,
  matching `chr2:1392274-1392337` on the **reverse** strand.
* `CP:i:1392337` is the **rightmost** coordinate of the chimeric block
  (reverse-strand block).
* `CT:Z:-` because the read is forward and the chimeric block is reverse.

### Example 2 — different chromosome, same strand

```
toy000005  FLAG=83  RNAME=chr1  POS=1949809  CIGAR=77M73S
           CG:Z:77S73M  CP:i:2539751  CT:Z:-  CC:Z:chr2
```

* Primary: `read[0:77]` matches `chr1:1949809-1949885` on the **forward** strand.
* Chimeric: `read[77:150]` =
  `ATGAACACCGTTGTCTGGCTATGCCCTGATACTCGCGTAACCGGCTACTGCTCAAGGTCGTGTAGTAAGAAAC`
  matches `chr2:2539751-2539823` on the **forward** strand.
* `CP:i:2539751` is the **leftmost** coordinate of the chimeric block
  (forward-strand block).
* `CT:Z:-` even though the chimeric block is forward, because the read itself is
  reverse-strand (`FLAG & 0x10`), so `CT` is flipped relative to the reference.

### Example 3 — same chromosome, different strand

```
toy000001  FLAG=113  RNAME=chr4  POS=2223048  CIGAR=24S46M80S
           CG:Z:24M126S  CP:i:631078  CT:Z:+  CC:Z:chr4
```

* Primary: `read[24:70]` matches `chr4:2223048-2223093` on the **forward** strand.
* Chimeric: `read[0:24]` = `CTACTCAGGTGTAACGGAATATAG`
  whose reverse complement is `CTATATTCCGTTACACCTGAGTAG`,
  matching `chr4:631055-631078` on the **reverse** strand.
* `CP:i:631078` is the **rightmost** coordinate of the chimeric block
  (reverse-strand block).
* `CT:Z:+` because the read is reverse-strand and the chimeric block is reverse,
  so the two cancel out.

### Example 4 — same chromosome, same strand

```
toy009797  FLAG=65  RNAME=chr1  POS=2112515  CIGAR=138M12S
           CG:Z:138S12M  CP:i:2680568  CT:Z:+  CC:Z:chr1
```

* Primary: `read[0:138]` matches `chr1:2112515-2112652` on the **forward** strand.
* Chimeric: `read[138:150]` = `GCGCCCTACATG`
  matches `chr1:2680568-2680579` on the **forward** strand.
* `CP:i:2680568` is the **leftmost** coordinate of the chimeric block
  (forward-strand block).
* `CT:Z:+` because the read is forward and the chimeric block is forward.

---

## Quick reference

| Case | `CG:Z` | `CP:i` points to | `CT:Z` |
|------|--------|------------------|--------|
| chimeric block forward | `...S<M>M` | leftmost base of `M` | `+` if read forward, `-` if read reverse |
| chimeric block reverse | `<M>M...S` | rightmost base of `M` | `-` if read forward, `+` if read reverse |

---

## Notes and caveats

1. **`CG:Z` is in read orientation, not reference orientation.** The `M` block
   may be at the start or the end of the string depending on which part of the
   read is chimeric.
2. **`CP:i` is the 5' end of the chimeric segment**, so it is the leftmost
   coordinate for a forward block and the rightmost coordinate for a reverse
   block. Do not assume it is always the leftmost coordinate.
3. **`CT:Z` is relative to the original read orientation**, not the reference
   forward strand. For reverse-strand reads it is the opposite of the reference
   strand of the chimeric block.
4. **Multiple chimeric sections.** A read spanning more than one fusion carries
   several `CG`/`CP`/`CT`/`CC` tag sets, in the order the sections appear in the
   read.
5. **Primary section.** The primary alignment is in the standard SAM columns and
   is not duplicated in `CG:Z`.
6. **Known issue (v2.2.0).** For chimeric reads whose primary alignment is the
   3' segment of the read, the BAM writer can emit the primary `CIGAR` in read
   orientation (e.g. `78M72S`) with `POS` at the rightmost base, instead of the
   reference orientation (`72S78M`) with `POS` at the leftmost base. In that
   case `CG:Z` and `CP:i` are also affected. See
   `FUSION_NEIGHBOUR_REMOVAL_PROBLEM.md` for the related investigation.
