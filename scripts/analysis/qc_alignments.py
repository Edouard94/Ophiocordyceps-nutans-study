#!/usr/bin/env python3
"""
QC filter for DNA alignments.

- Computes a simple majority-rule consensus per column (A/C/G/T only, ignoring gaps and N)
- For each sequence, computes:
  * non-gap A/C/G/T sites (sites_considered)
  * proportion of Ns among non-gap characters (n_prop)
  * mismatch proportion to consensus at considered sites (mm_prop)
- Removes sequences exceeding thresholds.

Usage:
  python3 scripts/qc_alignments.py \
    --input path/to/aln.fasta \
    --output path/to/filtered.fasta \
    --log logs/qc_gene.tsv \
    --max-n 0.5 \
    --max-mismatch 0.3 \
    --min-sites 100
"""

from __future__ import annotations
import argparse
from collections import Counter


def read_fasta(path: str) -> list[tuple[str, str]]:
    recs: list[tuple[str, str]] = []
    with open(path, "r", encoding="utf-8") as fh:
        hdr = None
        seq_parts = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(seq_parts)))
                hdr = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if hdr is not None:
            recs.append((hdr, "".join(seq_parts)))
    return recs


def write_fasta(path: str, recs: list[tuple[str, str]]):
    with open(path, "w", encoding="utf-8") as out:
        for h, s in recs:
            out.write(f">{h}\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i+60] + "\n")


def consensus(seq_list: list[str]) -> str:
    if not seq_list:
        return ""
    L = max((len(s) for s in seq_list), default=0)
    cols = []
    for i in range(L):
        cnt = Counter()
        for s in seq_list:
            if i >= len(s):
                continue
            c = s[i].upper()
            if c in ("A", "C", "G", "T"):
                cnt[c] += 1
        if not cnt:
            cols.append("N")
        else:
            cols.append(cnt.most_common(1)[0][0])
    return "".join(cols)


def per_seq_stats(seq: str, cons: str) -> tuple[int, int, int]:
    """Return (sites_considered, n_count, mismatches).
    sites_considered = positions where seq has A/C/G/T (ignore gaps and N in seq)
    n_count = positions where seq has N
    mismatches = positions where seq is A/C/G/T and differs from consensus A/C/G/T at same site
    """
    L = min(len(seq), len(cons))
    sites = 0
    ncount = 0
    mm = 0
    for i in range(L):
        q = seq[i].upper()
        if q == "N":
            ncount += 1
            continue
        if q not in ("A", "C", "G", "T"):
            continue  # ignore gaps and others
        sites += 1
        c = cons[i].upper() if i < len(cons) else "N"
        if c in ("A", "C", "G", "T") and q != c:
            mm += 1
    return sites, ncount, mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--log", required=True)
    ap.add_argument("--max-n", type=float, default=0.5, help="Max proportion of N among non-gap characters")
    ap.add_argument("--max-mismatch", type=float, default=0.3, help="Max mismatch proportion to consensus at considered sites")
    ap.add_argument("--min-sites", type=int, default=100, help="Minimum considered sites (A/C/G/T) required to evaluate a sequence")
    args = ap.parse_args()

    recs = read_fasta(args.input)
    if not recs:
        # nothing to do
        write_fasta(args.output, recs)
        return

    cons = consensus([s for _, s in recs])

    kept: list[tuple[str, str]] = []
    removed: list[tuple[str, float, float, int]] = []  # (hdr, n_prop, mm_prop, sites)

    for hdr, seq in recs:
        sites, ncount, mm = per_seq_stats(seq, cons)
        if sites < args.min_sites:
            # too little signal; keep but note that QC is inconclusive
            kept.append((hdr, seq))
            continue
        n_prop = ncount / (ncount + sites) if (ncount + sites) > 0 else 0.0
        mm_prop = mm / sites if sites > 0 else 0.0
        if n_prop > args.max_n or mm_prop > args.max_mismatch:
            removed.append((hdr, n_prop, mm_prop, sites))
        else:
            kept.append((hdr, seq))

    # Write outputs
    write_fasta(args.output, kept)

    # Log
    with open(args.log, "w", encoding="utf-8") as lf:
        lf.write("Header\tSitesConsidered\tN_Prop\tMismatch_Prop\tDecision\n")
        for hdr, seq in recs:
            sites, ncount, mm = per_seq_stats(seq, cons)
            n_prop = ncount / (ncount + sites) if (ncount + sites) > 0 else 0.0
            mm_prop = mm / sites if sites > 0 else 0.0
            decision = "REMOVE" if any(hdr == h for h, _, _, _ in removed) else ("KEEP_LOW_SITES" if sites < args.min_sites else "KEEP")
            lf.write(f"{hdr}\t{sites}\t{n_prop:.3f}\t{mm_prop:.3f}\t{decision}\n")

    print(f"QC: kept {len(kept)} / {len(recs)} sequences; removed {len(removed)} (see {args.log})")


if __name__ == "__main__":
    main()
