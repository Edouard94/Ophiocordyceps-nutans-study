#!/usr/bin/env python3
# Example usage:
# python3 scripts/compare_platforms.py --amplicon-dir results/amplicons --outdir results/comparisons
import argparse
from pathlib import Path
from collections import defaultdict
import re
import sys
import subprocess
import tempfile
import os

try:
    import numpy as np
    import pandas as pd
    import seaborn as sns  # type: ignore
    import matplotlib.pyplot as plt  # type: ignore
    HAS_PLOTTING = True
except Exception:
    HAS_PLOTTING = False
    import numpy as np
    import pandas as pd


def read_fasta(p: Path):
    name = None
    seq = []
    with open(p) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, ''.join(seq)


def id_parts(h):
    # Expect Phylo_ID like Platform_Sample_Marker or similar
    # We'll detect marker as one of ITS/LSU/SSU/TEF at the end
    marker_match = re.search(r"_(ITS|LSU|SSU|TEF)$", h, re.IGNORECASE)
    marker = marker_match.group(1).upper() if marker_match else None
    plat = h.split('_', 1)[0]
    # Sample: take everything between first underscore and last underscore before marker
    sample = None
    if marker_match:
        sample = h.split('_' + marker)[0]
        sample = sample.split('_', 1)[-1]
    return plat, sample, marker


def pident(a, b):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as qf, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as sf:
        qf.write(f">query\n{a}\n")
        sf.write(f">subject\n{b}\n")
        qf.flush()
        sf.flush()
        try:
            # Try standard blastn first
            result = subprocess.run(['blastn', '-query', qf.name, '-subject', sf.name, 
                                   '-outfmt', '6 pident', '-max_hsps', '1'], 
                                    capture_output=True, text=True, timeout=30)
            if result.returncode == 0 and result.stdout.strip():
                # Take the first line only
                first_line = result.stdout.strip().split('\n')[0]
                return float(first_line)
            # If no alignment, try with relaxed parameters
            result = subprocess.run(['blastn', '-query', qf.name, '-subject', sf.name, 
                                   '-outfmt', '6 pident', '-evalue', '10', '-word_size', '7', '-max_hsps', '1'], 
                                    capture_output=True, text=True, timeout=30)
            if result.returncode == 0 and result.stdout.strip():
                first_line = result.stdout.strip().split('\n')[0]
                return float(first_line)
        except (subprocess.TimeoutExpired, ValueError, FileNotFoundError) as e:
            print(f"Error in pident: {e}")
        finally:
            os.unlink(qf.name)
            os.unlink(sf.name)
    return np.nan


def main():
    ap = argparse.ArgumentParser(description="Compute platform agreement heatmaps from amplicon FASTAs.")
    ap.add_argument('--amplicon-dir', required=True, help='Directory with amplicon FASTAs per gene')
    ap.add_argument('--outdir', required=True, help='Output directory for CSVs and figures')
    args = ap.parse_args()

    amplicon_dir = Path(args.amplicon_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pairs = [('Illumina', 'Flongle'), ('Sanger', 'Flongle'), ('Sanger', 'Illumina')]

    all_rows = []
    presence_rows = []

    for gene in ['ITS', 'LSU', 'SSU', 'TEF']:
        f = amplicon_dir / f"{gene}_amplicons.fasta"
        if not f.exists():
            continue
        seqs = dict(read_fasta(f))

        # index by sample and platform
        by_sample = defaultdict(dict)
        for h, s in seqs.items():
            plat, sample, marker = id_parts(h)
            if marker != gene:
                continue
            if sample:
                by_sample[sample][plat] = s

        # compute identities for requested pairs and collect presence
        for sample, d in by_sample.items():
            for plat in ['Illumina', 'Sanger', 'Flongle']:
                if plat in d:
                    presence_rows.append({'Sample': sample, 'Marker': gene, 'Platform': plat})
            for p1, p2 in pairs:
                if p1 in d and p2 in d:
                    pid = pident(d[p1], d[p2])
                    all_rows.append({'Sample': sample, 'Marker': gene, 'Pair': f'{p1}_vs_{p2}', 'PercentIdentity': pid})

    if not all_rows:
        print("No comparable pairs found.")
        sys.exit(0)

    df = pd.DataFrame(all_rows)
    df.to_csv(outdir / 'platform_identity_long.csv', index=False)

    # Process presence data
    df_presence = pd.DataFrame(presence_rows)
    df_presence_agg = df_presence.groupby(['Sample', 'Marker'])['Platform'].apply(lambda x: '\n'.join(sorted(x))).reset_index()
    presence_pivot = df_presence_agg.pivot(index='Sample', columns='Marker', values='Platform').fillna('')
    presence_pivot.to_csv(outdir / 'amplicon_recovery.csv')

    # Heatmaps per pair: rows samples, columns markers
    if HAS_PLOTTING:
        for pair in set(df['Pair']):
            sub = df[df['Pair'] == pair]
            pivot = sub.pivot_table(index='Sample', columns='Marker', values='PercentIdentity')
            # Skip if insufficient data
            if pivot.shape[0] < 2 or pivot.shape[1] < 2:
                print(f"Skipping heatmap for {pair}: insufficient data (shape {pivot.shape})")
                pivot.to_csv(outdir / f'matrix_{pair}.csv')
                continue
            plt.figure(figsize=(6, max(3, 0.3 * len(pivot))))
            sns.heatmap(pivot, annot=True, fmt='.1f', cmap='viridis', vmin=0, vmax=100)
            plt.title(f'Percent identity: {pair}')
            plt.tight_layout()
            plt.savefig(outdir / f'heatmap_{pair}.png', dpi=200)
            plt.close()
            # Also save the matrix as CSV for R plotting
            pivot.to_csv(outdir / f'matrix_{pair}.csv')

        # Plot amplicon recovery
        color_pivot = presence_pivot.map(lambda x: len(x.split('\n')) if x else 0)
        plt.figure(figsize=(8, max(4, 0.5 * len(presence_pivot))))
        ax = sns.heatmap(color_pivot, annot=presence_pivot, fmt='', cmap='Blues', vmin=0, vmax=3, 
                         cbar_kws={'label': 'Number of platforms'}, annot_kws={'fontsize': 8, 'va': 'center'})
        plt.title('Amplicon recovery by platform')
        plt.tight_layout()
        plt.savefig(outdir / 'amplicon_recovery.png', dpi=200, bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    main()
