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

---

## Extended CIGAR encoding: `N`, `B`, `n`, `b`

### 1. What the extended CIGAR is

Before any tag is written, `subjunc` realigns the read and stores the result in
`current_result->cigar_string`. This is **not** a SAM CIGAR. Besides `M`, `I`,
`D` and `S` it may contain four *jump* operators that describe how the read
moves from one reference location to another:

| Operator | Case | Direction | Strand |
|----------|------|-----------|--------|
| `N` | upper | forward (increasing coordinate) | unchanged |
| `B` | upper | backward (decreasing coordinate) | unchanged |
| `n` | lower | forward | **flipped** |
| `b` | lower | backward | **flipped** |

The number in front of a jump operator is the **jump length in reference
bases**, not a read length. `M`, `I`, `D`, `S` keep their SAM meaning.

`N` and `B` are the ordinary splice / long-deletion operators. `n` and `b` are
the fusion operators: the read continues on the opposite strand.

### 2. The walk in `chimeric_cigar_parts()`

`src/core.c`:

```c
int chimeric_cigar_parts(global_context_t * global_context, unsigned int sel_pos,
                         int is_first_section_negative_strand,
                         int is_first_section_reversed,
                         char * in_cigar, unsigned int * out_poses,
                         char ** out_cigars, char * out_strands,
                         int read_len, short * perfect_lens, char * read_name)
```

The string is walked **left to right in reference order**, while a separate
counter tracks how many read bases have been consumed. The state is:

| Variable | Meaning |
|----------|---------|
| `read_cursor` | number of read bases consumed so far |
| `current_perfect_cursor` | reference cursor, used when `is_reversed == 0` |
| `current_perfect_map_start` | reference position of the current section's start |
| `is_reversed` | 1 = the reference cursor moves right-to-left |
| `is_negative` | 1 = the section is on the reference **reverse** strand |
| `tmpi` | the number currently being accumulated |

`is_reversed` and `is_negative` are **independent**. `is_reversed` is a
bookkeeping direction; `is_negative` is the biological strand. Both are flipped
by a lowercase jump, but they start from different values:

```c
chimeric_cigar_parts(global_context, r->linear_position,
                     is_first_section_jumped ^ current_strand,   /* is_negative */
                     is_first_section_jumped,                    /* is_reversed */
                     r->current_cigar_decompress, ...);
```

Initialisation of section 0:

```c
out_poses[0]   = current_perfect_map_start - (is_reversed ? 1 : 0);
out_strands[0] = is_negative ? '-' : '+';
```

### 3. Jump target

For a jump operator with length `len`:

```c
if(is_reversed)
    target = current_perfect_map_start - 1 + (toupper(ncch)=='N' ? +len : -len);
else
    target = current_perfect_cursor      + (toupper(ncch)=='N' ? +len : -len);
```

The `-1` in the `is_reversed` branch is needed because the cursor is already one
past the last consumed base.

### 4. When a jump actually starts a new section

```c
if(is_chro_jump || islower(ncch) || ncch == 'B' || is_long_jump)
```

Consequences:

* **Lowercase `n` / `b` always start a new section** (they are fusion jumps).
* **Uppercase `B` always starts a new section.**
* **Uppercase `N` starts a new section only if it changes chromosome or is
  longer than `maximum_intron_length`.** A short intra-chromosomal `N` is an
  ordinary intron: it is handled by the `D`/`N` branch, which advances the
  reference cursor but **not** `read_cursor`, and no new section is created.

### 5. Section bookkeeping

When a jump does start a new section:

```c
current_perfect_cursor = target;

if(islower(ncch)){
    is_reversed = !is_reversed;
    is_negative = !is_negative;
    if(is_reversed) current_perfect_cursor--;      /* one more step back */
}

current_perfect_map_start = current_perfect_cursor;
tmpi = 0;

/* close the old section with its trailing soft clip */
out_cigars[sec] += "<read_len - read_cursor>S";
perfect_lens[sec] = perfect_len;
perfect_len = 0;
sec++;

/* open the new section */
out_poses[sec]   = current_perfect_map_start - read_cursor;
out_strands[sec] = is_negative ? '-' : '+';
out_cigars[sec]  = "<read_cursor>S";               /* leading soft clip */
```

The critical line is

```c
out_poses[sec] = current_perfect_map_start - read_cursor;
```

This is **not** the block's leftmost base. It is the reference position that
read base 0 *would* have if the section's alignment were extended backwards to
the start of the read.

Then, for every `M`, `S`, `D` or `N` inside the section:

```c
if(ncch == 'M' || ncch == 'S'){
    read_cursor += tmpi;
    if(ncch == 'M') perfect_len += tmpi;
    if(is_reversed) out_poses[sec] += tmpi;        /* accumulate */
    else            current_perfect_cursor += tmpi;/* out_poses[sec] untouched */
}
else if(ncch == 'D' || ncch == 'N'){
    if(is_reversed) out_poses[sec] += tmpi;
    else            current_perfect_cursor += tmpi;
}
else if(ncch == 'I'){
    read_cursor += tmpi;                           /* reference cursor unchanged */
}
```

So `out_poses[sec]` means two different things:

* `is_reversed == 1` -> it accumulates over the section and ends at the block's
  **rightmost** base.
* `is_reversed == 0` -> it stays at `target - read_cursor`, i.e. the **virtual
  position of read base 0**.

Note that `S` also advances `out_poses[sec]` when `is_reversed == 1`, so a
trailing soft clip shifts the stored value.

### 6. Recovering the block's leftmost base

Two cases, depending on `is_reversed` **at the start of the section**:

```
is_reversed == 1:   leftmost = out_poses[sec] - (mlen - 1)
is_reversed == 0:   leftmost = out_poses[sec] + read_cursor_at_section_start
                                - (is_negative ? mlen - 1 : 0)
```

where

* `read_cursor_at_section_start` is the leading soft clip of the section's
  CIGAR (`out_cigars[sec]`), and
* `mlen` is the sum of `M`/`D`/`N` in that section.

The `- (mlen - 1)` term is needed because a reverse-strand block extends
**backwards** from the jump point.

### 7. Worked examples

All examples use the toy reference in
`UncommonCases/subjuncAlljunc/ToyExample` (4 chromosomes of 3,000,000 bp).
Coordinates are **0-based global** coordinates; `sel_pos` is the value passed
into `chimeric_cigar_parts()`.

#### Example 1 - `n` (forward jump, strand flip), chimeric block is the 2nd section

```
in_cigar = 80M3519224n70M      sel_pos = 647148      read_len = 150
entry: is_reversed = 1, is_negative = 0
```

| Operator | `read_cursor` | `current_perfect_cursor` | `out_poses` | `out_strands` | What happened |
|----------|---------------|--------------------------|-------------|---------------|---------------|
| init | 0 | 647148 | `[647147]` | `+` | `sel_pos - 1` because `is_reversed` |
| `80M` | 80 | 647148 | `[647227]` | `+` | `is_reversed` -> `out_poses[0] += 80` |
| `3519224n` | 80 | 4166371 | `[647227, 4166291]` | `+`, `-` | target `= 647148 - 1 + 3519224`; flip -> `is_reversed=0`, `is_negative=1`; `out_poses[1] = 4166371 - 80` |
| `70M` | 150 | 4166441 | `[647227, 4166291]` | `+`, `-` | `is_reversed=0` -> cursor advances, `out_poses[1]` unchanged |

Sections:

| sec | read interval | `CG` | strand | global ref interval | leftmost |
|-----|---------------|------|--------|---------------------|----------|
| 0 | `read[0:80]` | `80M70S` | `+` | 647148 - 647227 | 647148 |
| 1 | `read[80:150]` | `80S70M` | `-` | 4166302 - 4166371 | 4166302 |

For section 1: `leftmost = 4166291 + 80 - (70-1) = 4166302`. The 5' base of a
reverse block is its rightmost base, `4166371`.

#### Example 2 - `n` (forward jump, strand flip), chimeric block is the 1st section

```
in_cigar = 64M5483319n86M      sel_pos = 4395887      read_len = 150
entry: is_reversed = 1, is_negative = 1
```

| Operator | `read_cursor` | `current_perfect_cursor` | `out_poses` | `out_strands` |
|----------|---------------|--------------------------|-------------|---------------|
| init | 0 | 4395887 | `[4395886]` | `-` |
| `64M` | 64 | 4395887 | `[4395950]` | `-` |
| `5483319n` | 64 | 9879205 | `[4395950, 9879141]` | `-`, `+` |
| `86M` | 150 | 9879291 | `[4395950, 9879141]` | `-`, `+` |

Sections:

| sec | read interval | `CG` | strand | global ref interval | leftmost |
|-----|---------------|------|--------|---------------------|----------|
| 0 | `read[0:64]` | `64M86S` | `-` | 4395887 - 4395950 | 4395887 |
| 1 | `read[64:150]` | `64S86M` | `+` | 9879141 - 9879226 | 9879141 |

For section 0: `leftmost = 4395950 - (64-1) = 4395887`.

#### Example 3 - `B` (backward jump, same strand)

```
in_cigar = 83M3015162B52M15S   sel_pos = 9879241      read_len = 150
entry: is_reversed = 0, is_negative = 0
```

| Operator | `read_cursor` | `current_perfect_cursor` | `out_poses` | `out_strands` |
|----------|---------------|--------------------------|-------------|---------------|
| init | 0 | 9879241 | `[9879241]` | `+` |
| `83M` | 83 | 9879324 | `[9879241]` | `+` |
| `3015162B` | 83 | 6864162 | `[9879241, 6864079]` | `+`, `+` |
| `52M` | 135 | 6864214 | `[9879241, 6864079]` | `+`, `+` |
| `15S` | 150 | 6864229 | `[9879241, 6864079]` | `+`, `+` |

Sections:

| sec | read interval | `CG` | strand | global ref interval | leftmost |
|-----|---------------|------|--------|---------------------|----------|
| 0 | `read[0:83]` | `83M67S` | `+` | 9879241 - 9879323 | 9879241 |
| 1 | `read[83:135]` | `83S52M15S` | `+` | 6864162 - 6864213 | 6864162 |

For section 1: `leftmost = 6864079 + 83 - 0 = 6864162`. The trailing `15S` does
**not** move `out_poses[1]` because `is_reversed == 0`.

#### Example 4 - `N` (forward jump, same strand, chromosome change)

```
in_cigar = 64M9616328N86M      sel_pos = 378455       read_len = 150
entry: is_reversed = 0, is_negative = 0
```

| Operator | `read_cursor` | `current_perfect_cursor` | `out_poses` | `out_strands` |
|----------|---------------|--------------------------|-------------|---------------|
| init | 0 | 378455 | `[378455]` | `+` |
| `64M` | 64 | 378519 | `[378455]` | `+` |
| `9616328N` | 64 | 9994847 | `[378455, 9994783]` | `+`, `+` |
| `86M` | 150 | 9994933 | `[378455, 9994783]` | `+`, `+` |

Here `N` **does** start a new section because the target is on a different
chromosome (`is_chro_jump`). If the same `N` had stayed on the same chromosome
and been shorter than `maximum_intron_length`, no new section would have been
created and the `N` would simply have advanced `current_perfect_cursor`.

#### Example 5 - `n` with a trailing soft clip

```
in_cigar = 75M3642234n66M9S    sel_pos = 8284954      read_len = 150
entry: is_reversed = 1, is_negative = 1
```

| Operator | `read_cursor` | `current_perfect_cursor` | `out_poses` | `out_strands` |
|----------|---------------|--------------------------|-------------|---------------|
| init | 0 | 8284954 | `[8284953]` | `-` |
| `75M` | 75 | 8284954 | `[8285028]` | `-` |
| `3642234n` | 75 | 11927262 | `[8285028, 11927112]` | `-`, `+` |
| `66M` | 141 | 11927262 | `[8285028, 11927112]` | `-`, `+` |
| `9S` | 150 | 11927262 | `[8285028, 11927112]` | `-`, `+` |

Sections:

| sec | read interval | `CG` | strand | global ref interval | leftmost |
|-----|---------------|------|--------|---------------------|----------|
| 0 | `read[0:75]` | `75M75S` | `-` | 8284954 - 8285028 | 8284954 |
| 1 | `read[75:141]` | `75S66M9S` | `+` | 11927187 - 11927252 | 11927187 |

For section 1: `leftmost = 11927112 + 75 - 0 = 11927187`.

#### Example 6 - leading soft clip before an `n`

```
in_cigar = 1S112M526736n37M    sel_pos = 8390944      read_len = 150
entry: is_reversed = 0, is_negative = 0
```

| Operator | `read_cursor` | `current_perfect_cursor` | `out_poses` | `out_strands` |
|----------|---------------|--------------------------|-------------|---------------|
| init | 0 | 8390944 | `[8390944]` | `+` |
| `1S` | 1 | 8390945 | `[8390944]` | `+` |
| `112M` | 113 | 8391057 | `[8390944]` | `+` |
| `526736n` | 113 | 8917792 | `[8390944, 8917679]` | `+`, `-` |
| `37M` | 150 | 8917792 | `[8390944, 8917716]` | `+`, `-` |

Sections:

| sec | read interval | `CG` | strand | global ref interval | leftmost |
|-----|---------------|------|--------|---------------------|----------|
| 0 | `read[1:113]` | `1S112M37S` | `+` | 8390945 - 8391056 | 8390945 |
| 1 | `read[113:150]` | `113S37M` | `-` | 8917680 - 8917716 | 8917680 |

For section 1, `is_reversed == 1` at section start, so
`leftmost = out_poses[1] - (37-1) = 8917716 - 36 = 8917680`.

### 8. Primary selection and the section swap

After the walk:

```c
char main_piece_strand = (is_first_section_negative_strand == is_first_section_reversed) ? '+' : '-';

for(xk1 = 0; xk1 < section_count; xk1++)
    if(best_match < 0 ||
       (main_piece_strand == out_strands[xk1] && perfect_lens[xk1] > perfect_lens[best_match]))
        best_match = xk1;

if(best_match > 0) { /* swap section 0 and section best_match */ }
```

The primary section is therefore the **longest section on the main strand**, and
the four arrays (`out_poses`, `out_cigars`, `out_strands`, `perfect_lens`) are
**reordered**. After the swap, `out_poses[0]` / `out_cigars[0]` are no longer the
first section in reference order, and `out_poses[1]` is no longer necessarily the
second.

### 9. How the tags are derived

```c
unsigned int vitual_head_pos = r->out_poses[xk1];
char strand_xor = (r->out_strands[xk1] == '-') != is_second_read;

locate_gene_position_max(vitual_head_pos, &chromosome_table,
                         &chimaric_chr, &chimeric_pos, NULL, NULL,
                         0 + r->out_lens[xk1]);

int soft_clipping_movement = get_soft_clipping_length(r->out_cigars[xk1]);

SUBreadSprintf(..., "\tCG:Z:%s\tCP:i:%u\tCT:Z:%c\tCC:Z:%s",
               r->out_cigars[xk1],
               max(1, chimeric_pos + soft_clipping_movement + 1),
               strand_xor ? '-' : '+',
               chimaric_chr);
```

* `CG:Z` is `out_cigars[xk1]`, already in read orientation.
* `CC:Z` is the chromosome of `out_poses[xk1]`.
* `CT:Z` is `out_strands[xk1]` XOR the read's own strand (`is_second_read`).
* `CP:i` is `chimeric_pos + leading_soft_clip + 1`, where `chimeric_pos` is the
  chromosome offset of `out_poses[xk1]`.

Because `out_poses[xk1]` is a *virtual read-head* coordinate rather than the
block's leftmost base (see section 5), `CP:i` is only correct when the section
happens to satisfy `out_poses[xk1] == leftmost - leading_soft_clip`. This is the
root cause of the known `CP:i` issue described in the Notes section.

### 10. Summary of operator effects

| Operator | `read_cursor` | reference cursor | `out_poses[sec]` | `is_reversed` | `is_negative` | new section |
|----------|---------------|------------------|------------------|---------------|---------------|-------------|
| `M` | `+len` | `+len` if `!is_reversed` | `+len` if `is_reversed` | - | - | no |
| `S` | `+len` | `+len` if `!is_reversed` | `+len` if `is_reversed` | - | - | no |
| `I` | `+len` | unchanged | unchanged | - | - | no |
| `D` | unchanged | `+len` if `!is_reversed` | `+len` if `is_reversed` | - | - | no |
| `N` (short, same chr) | unchanged | `+len` if `!is_reversed` | `+len` if `is_reversed` | - | - | no |
| `N` (chr change or long) | unchanged | `= target` | reset | - | - | **yes** |
| `B` | unchanged | `= target` | reset | - | - | **yes** |
| `n` | unchanged | `= target` | reset | **flip** | **flip** | **yes** |
| `b` | unchanged | `= target` | reset | **flip** | **flip** | **yes** |
