#!/usr/bin/awk -f
#
# Summarize mq_forensics suffstats (*.advanced_MQs.tsv) to one row per BAM.
# Run as:
#   awk -f summarize_by_bam.awk *_advanced_MQs.tsv > mq_bam_summary.tsv
#
# Columns:
# sample n_sites_depth_ge1 depth_mean mismatch_rate ins_rate del_rate clip_rate \
# mq_mean mq_sd effmq_mean effmq_sd subQ_mean subQ_sd clipQ_mean clipQ_sd \
# clipfrac_mean clipfrac_sd frac_capped mean_delta sd_delta \
# flank_cov_mean flank_cov_sd flank_cf_mean flank_cf_sd \
# gc_frac_fwd gc_frac_rev frac_N_pooled strand_bias_z prop_reads_clipfrac_gt0.5
#

BEGIN {
    FS = OFS = "\t"

    # Print header once
    print "sample",
          "n_sites_depth_ge1",
          "depth_mean",
          "mismatch_rate","ins_rate","del_rate","clip_rate",
          "mq_mean","mq_sd",
          "effmq_mean","effmq_sd",
          "subQ_mean","subQ_sd",
          "clipQ_mean","clipQ_sd",
          "clipfrac_mean","clipfrac_sd",
          "frac_capped","mean_delta","sd_delta",
          "flank_cov_mean","flank_cov_sd",
          "flank_cf_mean","flank_cf_sd",
          "gc_frac_fwd","gc_frac_rev",
          "frac_N_pooled",
          "strand_bias_z",
          "prop_reads_clipfrac_gt0.5"

    SAMPLE_ENV = SAMPLE   # keep original env value if any
    reset_accumulators()
}

# ---------- helpers ----------

# derive sample name from filename if SAMPLE not set
function derive_sample_name(raw,   s) {
    s = raw
    sub(/^.*\//, "", s)                     # strip path
    sub(/_advanced_MQs\.tsv$/, "", s)       # common mqforensics suffix
    sub(/\.tsv$/, "", s)                    # fallback
    return s
}

# reset all per-sample accumulators
function reset_accumulators(  i) {
    n_sites = 0
    total_depth  = 0
    total_depth2 = 0

    tot_mm = tot_ins = tot_del = tot_clip = 0

    # pooled moments for MQ, effMQ, subQ, clipQ, clipfrac
    N_mq = S_mq = SS_mq = 0
    N_eff = S_eff = SS_eff = 0
    N_subQ = S_subQ = SS_subQ = 0
    N_clipQ = S_clipQ = SS_clipQ = 0
    N_cf = S_cf = SS_cf = 0

    # capping/delta pooled across reads
    N_capped = 0
    S_delta = SS_delta = 0

    # flanking suffstats
    N_flank_cov = S_flank_cov = SS_flank_cov = 0
    N_flank_cf  = S_flank_cf  = SS_flank_cf  = 0

    # per-strand base counts (all reads)
    A_fwd = C_fwd = G_fwd = T_fwd = N_fwd = 0
    A_rev = C_rev = G_rev = T_rev = N_rev = 0

    # clipfrac histogram over reads (from hist_clipfrac)
    for (i = 0; i < 10; i++) clipHist[i] = 0

    current_sample = ""
}

# pooled mean/sd helpers
function pooled_mean(N,S) {
    return (N > 0 ? S / N : "NA")
}
function pooled_sd(N,S,SS,   m,v) {
    if (N <= 0) return "NA"
    m = S / N
    v = SS / N - m * m
    if (v < 0) v = 0
    return sqrt(v)
}

# emit summary line for current sample (if any)
function flush_sample(   depth_mean,depth_var, \
                         mismatch_rate,ins_rate,del_rate,clip_rate, \
                         mq_mean,mq_sd,effmq_mean,effmq_sd, \
                         subQ_mean,subQ_sd,clipQ_mean,clipQ_sd, \
                         cf_mean,cf_sd, \
                         frac_capped,mean_delta,sd_delta, \
                         flank_cov_mean,flank_cov_sd, \
                         flank_cf_mean,flank_cf_sd, \
                         num_gc_fwd,den_gc_fwd,gc_frac_fwd, \
                         num_gc_rev,den_gc_rev,gc_frac_rev, \
                         total_N,total_all,frac_N_pooled, \
                         depth_fwd,depth_rev,total_depth_reads,strand_bias_z, \
                         total_clip_reads,reads_highclip,prop_clip_gt05, \
                         vfc,vff,var_delta,i) {

    if (current_sample == "") return

    if (n_sites == 0) {
        # No sites with depth>=1; still emit a row of NAs
        print current_sample,
              0,
              "NA",
              "NA","NA","NA","NA",
              "NA","NA",
              "NA","NA",
              "NA","NA",
              "NA","NA",
              "NA","NA",
              "NA","NA","NA",
              "NA","NA",
              "NA","NA",
              "NA","NA",
              "NA",
              "NA",
              "NA",
              "NA",
              "NA"
        return
    }

    # depth stats across sites
    depth_mean = total_depth / n_sites
    depth_var  = (total_depth2 / n_sites) - depth_mean * depth_mean
    if (depth_var < 0) depth_var = 0

    # bp-level rates
    mismatch_rate = (total_depth > 0 ? tot_mm   / total_depth : "NA")
    ins_rate      = (total_depth > 0 ? tot_ins  / total_depth : "NA")
    del_rate      = (total_depth > 0 ? tot_del  / total_depth : "NA")
    clip_rate     = (total_depth > 0 ? tot_clip / total_depth : "NA")

    mq_mean      = pooled_mean(N_mq, S_mq)
    mq_sd        = pooled_sd(N_mq, S_mq, SS_mq)

    effmq_mean   = pooled_mean(N_eff, S_eff)
    effmq_sd     = pooled_sd(N_eff, S_eff, SS_eff)

    subQ_mean    = pooled_mean(N_subQ, S_subQ)
    subQ_sd      = pooled_sd(N_subQ, S_subQ, SS_subQ)

    clipQ_mean   = pooled_mean(N_clipQ, S_clipQ)
    clipQ_sd     = pooled_sd(N_clipQ, S_clipQ, SS_clipQ)

    cf_mean      = pooled_mean(N_cf, S_cf)
    cf_sd        = pooled_sd(N_cf, S_cf, SS_cf)

    if (total_depth > 0) {
        frac_capped = N_capped / total_depth

        mean_delta  = S_delta / total_depth
        var_delta   = SS_delta / total_depth - mean_delta * mean_delta
        if (var_delta < 0) var_delta = 0
        sd_delta    = sqrt(var_delta)
    } else {
        frac_capped = "NA"
        mean_delta  = "NA"
        sd_delta    = "NA"
    }

    # flanking stats
    if (N_flank_cov > 0) {
        flank_cov_mean = S_flank_cov / N_flank_cov
        vfc = SS_flank_cov / N_flank_cov - flank_cov_mean * flank_cov_mean
        if (vfc < 0) vfc = 0
        flank_cov_sd = sqrt(vfc)
    } else {
        flank_cov_mean = "NA"
        flank_cov_sd   = "NA"
    }

    if (N_flank_cf > 0) {
        flank_cf_mean = S_flank_cf / N_flank_cf
        vff = SS_flank_cf / N_flank_cf - flank_cf_mean * flank_cf_mean
        if (vff < 0) vff = 0
        flank_cf_sd = sqrt(vff)
    } else {
        flank_cf_mean = "NA"
        flank_cf_sd   = "NA"
    }

    # GC and N fractions from base counts (all reads)
    num_gc_fwd = C_fwd + G_fwd
    den_gc_fwd = A_fwd + C_fwd + G_fwd + T_fwd
    gc_frac_fwd = (den_gc_fwd > 0 ? num_gc_fwd / den_gc_fwd : "NA")

    num_gc_rev = C_rev + G_rev
    den_gc_rev = A_rev + C_rev + G_rev + T_rev
    gc_frac_rev = (den_gc_rev > 0 ? num_gc_rev / den_gc_rev : "NA")

    total_N = N_fwd + N_rev
    total_all = A_fwd + C_fwd + G_fwd + T_fwd + N_fwd \
              + A_rev + C_rev + G_rev + T_rev + N_rev
    frac_N_pooled = (total_all > 0 ? total_N / total_all : "NA")

    depth_fwd = A_fwd + C_fwd + G_fwd + T_fwd + N_fwd
    depth_rev = A_rev + C_rev + G_rev + T_rev + N_rev
    total_depth_reads = depth_fwd + depth_rev
    if (total_depth_reads > 0) {
        # simple symmetric z for strand imbalance
        strand_bias_z = (depth_fwd - depth_rev) / sqrt(total_depth_reads)
    } else {
        strand_bias_z = "NA"
    }

    # prop_reads_clipfrac_gt0.5 from hist_clipfrac (bins 0.5–1.0)
    total_clip_reads = 0
    reads_highclip = 0
    for (i = 0; i < 10; i++) {
        total_clip_reads += clipHist[i]
        if (i >= 5) reads_highclip += clipHist[i]
    }
    if (total_clip_reads > 0) {
        prop_clip_gt05 = reads_highclip / total_clip_reads
    } else {
        prop_clip_gt05 = "NA"
    }

    print current_sample,
          n_sites,
          depth_mean,
          mismatch_rate, ins_rate, del_rate, clip_rate,
          mq_mean, mq_sd,
          effmq_mean, effmq_sd,
          subQ_mean, subQ_sd,
          clipQ_mean, clipQ_sd,
          cf_mean, cf_sd,
          frac_capped, mean_delta, sd_delta,
          flank_cov_mean, flank_cov_sd,
          flank_cf_mean, flank_cf_sd,
          gc_frac_fwd, gc_frac_rev,
          frac_N_pooled,
          strand_bias_z,
          prop_clip_gt05
}

# ---------- main logic ----------

# New file: flush previous sample, then reset
FNR == 1 {
    if (NR > 1) {
        flush_sample()
        reset_accumulators()
    }

    if (SAMPLE_ENV != "") {
        current_sample = SAMPLE_ENV
    } else {
        current_sample = derive_sample_name(FILENAME)
    }

    next
}

# skip header in each file
$1 == "chrom" && $2 == "pos" { next }

# Per-line accumulation from suffstats
{
    depth = $3 + 0
    if (depth < 1) next

    n_sites++
    total_depth  += depth
    total_depth2 += depth * depth

    mm    = $4 + 0
    ins   = $5 + 0
    del   = $6 + 0
    clips = $7 + 0

    tot_mm   += mm
    tot_ins  += ins
    tot_del  += del
    tot_clip += clips

    # six blocks of (n, sum, sumsq):
    # mq:      8–10
    # cap:     11–13
    # eff:     14–16
    # subQ:    17–19
    # clipQ:   20–22
    # clipfrac:23–25

    n_mq   = $8 + 0
    s_mq   = $9 + 0
    ss_mq  = $10 + 0
    N_mq  += n_mq
    S_mq  += s_mq
    SS_mq += ss_mq

    n_eff   = $14 + 0
    s_eff   = $15 + 0
    ss_eff  = $16 + 0
    N_eff  += n_eff
    S_eff  += s_eff
    SS_eff += ss_eff

    n_subQ   = $17 + 0
    s_subQ   = $18 + 0
    ss_subQ  = $19 + 0
    N_subQ  += n_subQ
    S_subQ  += s_subQ
    SS_subQ += ss_subQ

    n_clipQ   = $20 + 0
    s_clipQ   = $21 + 0
    ss_clipQ  = $22 + 0
    N_clipQ  += n_clipQ
    S_clipQ  += s_clipQ
    SS_clipQ += ss_clipQ

    n_cf   = $23 + 0
    s_cf   = $24 + 0
    ss_cf  = $25 + 0
    N_cf  += n_cf
    S_cf  += s_cf
    SS_cf += ss_cf

    # capping / delta
    N_capped += ($26 + 0)
    S_delta  += ($27 + 0)
    SS_delta += ($28 + 0)

    # flanking suffstats
    n_fc   = $29 + 0
    s_fc   = $30 + 0
    ss_fc  = $31 + 0
    N_flank_cov  += n_fc
    S_flank_cov  += s_fc
    SS_flank_cov += ss_fc

    n_ff   = $32 + 0
    s_ff   = $33 + 0
    ss_ff  = $34 + 0
    N_flank_cf  += n_ff
    S_flank_cf  += s_ff
    SS_flank_cf += ss_ff

    # per-strand base counts (all reads)
    A_fwd += ($35 + 0)
    C_fwd += ($36 + 0)
    G_fwd += ($37 + 0)
    T_fwd += ($38 + 0)
    N_fwd += ($39 + 0)

    A_rev += ($40 + 0)
    C_rev += ($41 + 0)
    G_rev += ($42 + 0)
    T_rev += ($43 + 0)
    N_rev += ($44 + 0)

    # hist_clipfrac: last 10 columns if --emit-hist
    if (NF >= 10) {
        start = NF - 9
        b = 0
        for (i = start; i <= NF; i++) {
            clipHist[b] += ($i + 0)
            b++
        }
    }
}

END {
    flush_sample()
}
