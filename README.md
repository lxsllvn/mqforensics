# mqforensics

Fast C/htslib tool to extract per-site and per-interval mapping statistics from BAMs. Designed for RADseq/GBS/amplicon panels (short loci) and scalable to thousands of BAMs.

---

# Overview of per-site statistics

**Basic coverage and error metrics**

|Column|Meaning|
|---|---|
|depth|Total reads overlapping the site (includes deletions, like mpileup).|
|mismatch_rate|Fraction of covered bases that disagree with the reference (mismatch_bases / depth).|
|ins_rate|Inserted bases per covered base (sum of insertion lengths / depth).|
|del_rate|Deleted bases per covered base (sum of deletion lengths / depth).|
|clip_rate|Fraction of clipped (soft + hard) bases relative to coverage.|
|strand_bias_z| Binomial-approximation Z-score for read depth imbalance, calculated from all aligned reads.|

---

**Mapping quality (MQ) distributions**

|Column|Meaning|
|---|---|
|mq_mean, mq_sd|Mean ± SD of MAPQ emitted by the aligner.|
|capmq60_mean, capmq60_sd|Mean ± SD of the capped MAPQ after invoking the -C mapping-quality cap, rescaled to 0–60. |
|effmq_mean, effmq_sd|Mean ± SD of effective MAPQ used after applying the -C mapping-quality cap (0-60). |
|frac_capped|Fraction of reads with a capped MAPQ. |
|mean_delta, sd_delta|Mean ± SD of (raw − eff) MAPQ reductions. |
|mq_median_hist, effmq_median_hist, clipfrac_median_hist|Histogram-based medians for MAPQ, eff MQ, and clip-fraction respectively|

---

**Mismatch- and clipping-derived features features**

|Column|Meaning|
|---|---|
|subQ_mean, subQ_sd|Mean ± SD of per-read substitution quality sum (capped base-quality scores ≥ 13).|
|clipQ_mean, clipQ_sd|Mean ± SD of per-read clip quality sum (soft-clip qsum + 13×hardclip length). |
|clipfrac_mean, clipfrac_sd|Mean ± SD of the fraction of clipped bases per read (soft + hard / read length).|

**MAPQ-distribution divergence metrics**

|Column|Meaning|
|---|---|
|ks_mq_eff|Kolmogorov–Smirnov distance between raw MQ and eff MQ histograms. |
|w1_mq_eff|Wasserstein-1 (Earth-Mover’s) distance between raw MQ and eff MQ histograms. |
|js_mq_eff|Jensen–Shannon divergence between raw MQ and eff MQ histograms.|

---

 **Sequence-entropy and base-composition metrics**

These metrics are based on only bases (BQ ≥ 13) matching the reference in order to describe local sequence complexity rather than population genetic variation. 

| Column                                          | Meaning |
|------------------------------------------------|---------|
| entropy_ref_pooled, entropy_ref_fwd, entropy_ref_rev | Shannon entropy (in nats) of base composition among reference-matching A/C/G/T bases with BQ ≥ 13, pooled across strands, and separately for forward and reverse strands. |
| alph_eff_ref_pooled, alph_eff_ref_fwd, alph_eff_ref_rev | Effective alphabet size for the same ref-matching bases. Defined as exp(entropy), interpreted as the equivalent number of equally frequent bases (range 1–4). |
| gc_frac_pooled, gc_frac_fwd, gc_frac_rev       | Fraction of G + C among ref-matching A/C/G/T bases (BQ ≥ 13), pooled across strands, and separately for forward and reverse strands. |

---

**Flanking-region context**

Flanking-region context computes per-site windowed means of coverage and clip-fraction in the surrounding region, optionally restricted to reference-matching bases only.

At present, the `--ref-only` flag used to invoke the entropy-related metrics also applies to the flanking-region context.  This means that you can't get the entropy/GC metrics *and* all-bases flanking-region context in the same run.

|Column|Meaning|
|---|---|
|flank_cov_mean, flank_cov_sd|Mean ± SD of flanking coverage (within ± --flank bp window) — gives local coverage homogeneity.|
|flank_cf_mean, flank_cf_sd|Mean ± SD of flanking clip-fraction — whether reads around the site are consistently clipped. Sensitive to local repeats.|

---

# Usage

In `analyze` mode, `MQforensics` emits sufficient triplets and histograms for a single BAM alignment.  Per-site statistics are obtained by first running `analyze` on each BAM alignment and then supplying the output to `summarize` mode. Per-sample statistics can also be obtained from the `analyze` output using the supplied `awk` script. 

This two-step framework makes `MQforensics` fast and fairly memory efficient (≤ 6 Gb of RAM on my BAMs), although very high read depths (*ca.* 500x over intervals summing to 1Mbp) will require more. 

At present, `mq_forensics` uses a maximum of 500 reads by default, which can be turned off completely with `--max-reads-per-site 0`  or changed by passing any other integer to `--max-reads-per-site`.

## Step 1. Extract metrics from one BAM

```bash
mq_forensics analyze -b <in.bam> -r <regions.bed> -f <ref.fa> \
-C 100 -d 3 --emit-suffstats --emit-hist --flank 10 --ref-only \
-o per_site.tsv -O per_interval.tsv 
```

Required:
- `b` : input BAM (indexed .bai required)
- `r` : BED (0-based, half-open)
- `C` : cap threshold (samtools-style -C), e.g. 100
- `d` : minimum depth for site to emit numeric outputs

Technically optional but not really:
`--emit-suffstats` : emit n, sum, sumsq triplets instead of direct means/medians. 

Without this flag, `mqforensics` will attempt to store every read-level QC value at every site. This might not work,and if it does, it will be very slow and memory-intensive for non-trivial BAMs. In the future, `--emit-suffstats` will be the default mode. 

Optional: 
`--emit-hist`      :  emit histograms for approximating divergences and medians in `summarize` mode
`-f/--fasta ref.fa`:  FASTA reference; required for `--ref-only`  (entropy, GC) and flanking-context features (--flank)
`--flank N`         : Enable flanking-context metrics using a ±N bp window inside each BED interval.
`--ref-only`       :  Enables the GC and entropy-based metrics. If `--flank` is also invoked, it will limit the flanking-region context to reference-matching bases but has no effect on any other metric (i.e.,  MAPQ, mismatch_rate, clipfrac_mean, clip_rate and so on)

## Step 2. Per-site summary statistics across samples

To get the per-site statistics, feed the `analyze` outputs to `summarize` mode: 

```bash
{ head -n1 *.site.tsv | head -n1
  tail -q -n +2 *.site.tsv \
  | LC_ALL=C sort -S 3G --parallel=1 -k1,1 -k2,2n
} \
| mq_forensics summarize -i - -o per_site_summary.tsv
```

Required:
- `-i` : sorted suffstats TSVs;
- `-o` : outfile  name 

---

# Build

Requires:
- GCC or Clang with C11 support
- htslib v1.18+ installed and discoverable via pkg-config

```bash
make
```

or manual:
```bash
gcc -O3 -std=c11 -o mq_forensics src/*.c \
  $(pkg-config --cflags --libs htslib) -lm
```

---

# Important caveats 

`mqforensics` calculates everything I could come up with; these metrics may not be sensible, robust, or useful. 

I have only tested BWA-MEM alignments. 

`--direct-mode` is likely to be very slow and/or memory-intensive, especially if you do not limit the number of reads. Use `--emit-suffstats`.

---
