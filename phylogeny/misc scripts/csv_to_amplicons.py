#!/usr/bin/env python3
import sys
import csv
import argparse
from pathlib import Path
import re

GENES = {"ITS", "SSU", "LSU", "TEF", "RPB1"}

def clean_header(s: str) -> str:
    # Keep simple, remove spaces/semicolons, enforce ASCII safe characters
    s = s.strip()
    s = s.replace(" ", "_").replace(";", "_")
    s = re.sub(r"[^A-Za-z0-9_\-\|\.]+", "_", s)
    return s

def clean_seq(seq: str) -> str:
    seq = re.sub(r"[^ACGTNacgtn]", "N", seq)
    return seq.upper()

def autodetect_delimiter(sample: str) -> str:
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        return dialect.delimiter
    except Exception:
        # Fall back to semicolon if it looks like it, else comma
        if sample.count(";") > sample.count(","):
            return ";"
        return ","

def main():
    p = argparse.ArgumentParser(description="Extract amplicon FASTAs per marker from a CSV.")
    p.add_argument("csv", help="Input CSV path (Illumina/Flongle/Sanger combined)")
    p.add_argument("outdir", help="Output directory for amplicon FASTAs")
    args = p.parse_args()

    in_path = Path(args.csv)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_data = in_path.read_text(encoding="utf-8", errors="replace")
    delim = autodetect_delimiter(sample_data[:10000])

    # Open file and iterate rows
    with open(in_path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        # Normalize fieldnames
        field_map = {name.lower(): name for name in reader.fieldnames or []}
        required = ["marker", "phylo_id", "sequence"]
        for r in required:
            if r not in field_map:
                sys.stderr.write(f"Error: required column '{r}' not found in CSV header. Columns: {reader.fieldnames}\n")
                sys.exit(1)

        # Open writers per gene lazily
        handles = {}
        counts = {g: 0 for g in GENES}
        total = 0

        for row in reader:
            marker = (row[field_map["marker"]] or "").strip().upper()
            if marker not in GENES:
                continue
            seq = (row[field_map["sequence"]] or "").strip()
            if not seq:
                continue
            phylo_id = (row[field_map["phylo_id"]] or "").strip()
            if not phylo_id:
                # Fallback to SeqID or Sample if missing
                sid = row.get(field_map.get("seqid", ""), "").strip() or \
                      (row.get(field_map.get("sample", ""), "").strip() + "_" + marker)
                phylo_id = sid

            header = clean_header(phylo_id)
            seq = clean_seq(seq)

            if marker not in handles:
                fp = out_dir / f"{marker}_amplicons.fasta"
                handles[marker] = open(fp, "w", encoding="utf-8")
            handles[marker].write(f">{header}\n{seq}\n")
            counts[marker] += 1
            total += 1

        for h in handles.values():
            h.close()

    sys.stderr.write("Wrote amplicons per marker:\n")
    for g in GENES:
        sys.stderr.write(f"  {g}: {counts[g]}\n")
    sys.stderr.write(f"Total sequences: {total}\n")

if __name__ == "__main__":
    main()
