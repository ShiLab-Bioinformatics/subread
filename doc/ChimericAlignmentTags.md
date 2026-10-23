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
| `CP:i` | i    | 1-based coordinate calculated from the stored section position and leading soft clip; it is not reliably the section's 5' base |
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

`CP:i` is calculated from the section position stored by `chimeric_cigar_parts()`.
That position is not consistently the leftmost or rightmost aligned base, so
`CP:i` must not be interpreted generally as the 5' base of the segment.

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

`chimeric_pos` is the chromosome offset returned for `out_poses[xk1]` by
`locate_gene_position_max()`. Adding the leading soft clip does not always
recover the first aligned base; see the extended CIGAR section below.

---

## `CT:Z` — chimeric strand

`CT:Z` is the strand of the chimeric section, `+` or `-`.

The writer computes it as follows (`is_second_read` identifies mate 2, not
the read's strand):

```
CT_is_negative = (out_strands[xk1] == '-') XOR is_second_read
```

Thus mate 1 uses `out_strands` directly and mate 2 inverts it. `FLAG & 0x10`
is not the operand in this calculation. Interpreting `CT` as a biological
strand requires accounting for that mate-dependent convention.

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
* The recorded `CT:Z:-` is the writer's mate-adjusted `out_strands` value.

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
* The recorded `CT:Z:-` reflects the writer's mate-dependent inversion; the
  code does not use `FLAG & 0x10` to compute this tag.

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
* The recorded `CT:Z:+` is the writer's mate-adjusted `out_strands` value.

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
* The recorded `CT:Z:+` is the writer's mate-adjusted `out_strands` value.

---

## Quick reference

| Field | Writer value |
|-------|--------------|
| `CG:Z` | selected section's generated CIGAR |
| `CP:i` | `max(1, chromosome_offset(out_poses) + leading_soft_clip + 1)` |
| `CT:Z` | `'-'` when `(out_strands == '-') != is_second_read`, otherwise `'+'` |

---

## Notes and caveats

1. **`CG:Z` is in read orientation, not reference orientation.** The `M` block
   may be at the start or the end of the string depending on which part of the
   read is chimeric.
2. **`CP:i` is a computed coordinate**, not a guaranteed 5' or leftmost base.
3. **`CT:Z` uses mate identity** to invert the stored section strand for mate 2.
4. **Multiple chimeric sections.** A read spanning more than one fusion carries
   several `CG`/`CP`/`CT`/`CC` tag sets. The primary-selection swap can change
   their order from the original CIGAR walk.
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

The string is walked **left to right in CIGAR-string order**, while a separate
counter tracks how many read bases have been consumed. The state is:

| Variable | Meaning |
|----------|---------|
| `read_cursor` | number of read bases consumed so far |
| `current_perfect_cursor` | reference cursor, used when `is_reversed == 0` |
| `current_perfect_map_start` | reference position of the current section's start |
| `is_reversed` | 1 = `M`/`S`/`D`/unsplit `N` accumulate in `out_poses` instead of `current_perfect_cursor` |
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

In the `is_reversed` branch the anchor is the fixed section start minus one;
`current_perfect_cursor` is not advanced by `M` in that branch.

### 4. When a jump actually starts a new section

```c
if(is_chro_jump || islower(ncch) || ncch == 'B' || is_long_jump)
```

Consequences:

* **Lowercase `n` / `b` always start a new section** (they are fusion jumps).
* **Uppercase `B` always starts a new section.**
* **Uppercase `N` starts a new section only if it changes chromosome or
  `abs(current_perfect_cursor - target) >= maximum_intron_length`.** This
  distance need not equal the encoded `N` length when `is_reversed` is true.
  A shorter intra-chromosomal `N` is an
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

* `is_reversed == 1` -> it accumulates the lengths of processed `M`, `S`, `D`
  and unsplit `N`, starting from `map_start - read_cursor_at_entry`. It is not
  generally the block's rightmost base.
* `is_reversed == 0` -> it stays at `target - read_cursor`, i.e. the **virtual
  position of read base 0**.

Note that `S` also advances `out_poses[sec]` when `is_reversed == 1`, so a
trailing soft clip shifts the stored value.

### 6. Position bookkeeping and its limits

For a later section, let `entry_read_cursor` be the value when the section
opens and `map_start` be its adjusted jump target. The code initializes
`out_poses[sec] = map_start - entry_read_cursor`. If `is_reversed` is false,
that value remains unchanged. If true, each subsequent `M`, `S`, `D`, and
unsplit `N` adds its length to `out_poses[sec]`; `I` does not. The initial
section instead starts at `sel_pos - (is_reversed ? 1 : 0)`.

For either lowercase jump (`n` or `b`), the code flips `is_reversed` and
subtracts one from the jump target **only if the new value is true**. The
likely purpose is an inclusive/exclusive boundary adjustment: the calculated
target acts as the boundary after the first base of a section traversed in
reverse, and `target - 1` selects that base. This is an inference from the
arithmetic, not a coordinate convention explicitly documented by the source.
The adjustment depends on the new traversal direction, not on whether the
jump was forward (`n`) or backward (`b`).

These are the actual stored values, not a general formula for the leftmost
aligned base. In particular, `is_negative` affects `out_strands`, not position
arithmetic. A caller needs the section's CIGAR, entry state and orientation
to infer its genomic interval. The tag writer does not perform that recovery.

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

Before primary selection, the generated CIGARs are `80M70S` and `80S70M`.
Their stored positions are `647227` and `4166291`; the second section's
non-reversed cursor walk begins at `4166371`. The section strand flags are
`+` and `-`. These stored positions are not a general leftmost-base formula.

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

Before primary selection, the generated CIGARs are `64M86S` and `64S86M`.
The stored positions are `4395950` and `9879141`; the second section's
non-reversed cursor walk starts at `9879205`. The strand flags are `-` and
`+`. The code does not convert the first stored position to `4395887`.

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

Before primary selection, the generated CIGARs are `83M67S` and
`83S52M15S`, with stored positions `9879241` and `6864079` and strand flags
`+` and `+`. The second cursor starts at `6864162`; its trailing `15S` moves
that cursor but leaves `out_poses[1]` unchanged.

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
| `3642234n` | 75 | 11927187 | `[8285028, 11927112]` | `-`, `+` |
| `66M` | 141 | 11927253 | `[8285028, 11927112]` | `-`, `+` |
| `9S` | 150 | 11927262 | `[8285028, 11927112]` | `-`, `+` |

Before primary selection, the generated CIGARs are `75M75S` and
`75S66M9S`, with stored positions `8285028` and `11927112` and strand flags
`-` and `+`. The second cursor starts at `11927187`; `66M9S` advances it to
`11927262`, while its stored position remains `11927112`.

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

Before primary selection, the generated CIGARs are `1S112M37S` and
`113S37M`, with stored positions `8390944` and `8917716` and strand flags
`+` and `-`. The second section starts with `out_poses[1] = 8917679`;
because it is reversed, `37M` adds 37 to that stored value.

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

Section 0 is the initial choice. A later section replaces it only when that
section is on `main_piece_strand` and has a strictly greater `M` count than
the current choice. Therefore the primary is not necessarily the longest
section on the main strand. After selection,
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
* `CT:Z` is `out_strands[xk1]` XOR mate identity (`is_second_read`).
* `CP:i` is `chimeric_pos + leading_soft_clip + 1`, where `chimeric_pos` is the
  chromosome offset of `out_poses[xk1]`.

`out_poses[xk1]` is not generally the block's leftmost base (see section 6),
so this calculation does not generally yield the first aligned base. The
chromosome lookup also uses `out_lens[xk1]` as its final argument.

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
| `n` | unchanged | `= target - 1` if newly reversed, else `target` | reset | **flip** | **flip** | **yes** |
| `b` | unchanged | `= target - 1` if newly reversed, else `target` | reset | **flip** | **flip** | **yes** |
