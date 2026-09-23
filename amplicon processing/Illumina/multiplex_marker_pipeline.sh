#!/bin/bash

# Multiplex Marker Extraction Pipeline for Consensus Barcode Recovery
# Workflow: QC -> marker-specific primer trim -> merge/filter per marker -> per-sample ASV consensus (UNOISE3) -> taxonomy
# Tools: fastp, cutadapt, vsearch, seqkit, blastn, taxonkit

# Multiplex primer sequences (same for all samples):
# Forward primers:
# ITS2f_350bp	GCATCGATGAAGAACGCAGC
# SSUf_305bp	CCGTGGTAATTCTAGAGCTAATACATGC
# LSUf_308bp	AGCGCACAAGTAGAGTGATCG
# TEFf_268bp	TGGTACAAGGGYTGGGAGAAGG
#
# Reverse primers:
# ITS2r_350bp	YTTTTCCTCCGCTTATTGATATGC
# SSUr_305bp	RTCGGGATTGGGTAATTTGCGC
# LSUr_308bp	GGTCCGTGTTTCAAGACGG
# TEFr_268bp	GCTGCTCGTGGTGCATYTCV

set -euo pipefail

# ---------------------------
# Config and helpers
# ---------------------------
THREADS=${THREADS:-20}
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
DEBUG=${DEBUG:-0}

if [[ "$DEBUG" == "1" ]]; then
  set -x
fi

# ASV consensus settings (UNOISE3-like)
# Keep top N sequences per sample (by derep size)
TOP_N_ASVS=${TOP_N_ASVS:-2}
# Quality filter thresholds for FASTQ -> FASTA step
MAX_EE=${MAX_EE:-3}
# Enforce minimum read length early; adjust via MIN_LEN or DEREP_MIN_LEN as needed
MIN_LEN=${MIN_LEN:-150}
# Minimum length for dereplication input (per-sample, per-marker)
DEREP_MIN_LEN=${DEREP_MIN_LEN:-$MIN_LEN}
# NCBI nt BLAST database and TaxonKit taxonomy database (optional taxonomy assignment)
NTDB=${NTDB:-/data2/nt_NCBI_database/nt}
# BLAST tuning: task and taxonomic scope (default fungi taxid=4751)
BLAST_TASK=${BLAST_TASK:-megablast}
BLAST_TAXIDS=${BLAST_TAXIDS:-4751}
TAXONKIT_DB=${TAXONKIT_DB:-$HOME/taxdump}

ahelp() {
    python3 -c "import os,sys; print(os.path.abspath(sys.argv[1]))" "$1"
}

# Reverse complement helper with IUPAC support
reverse_complement() {
  # IUPAC complements: A<->T, C<->G, R<->Y, W<->W, S<->S, K<->M, B<->V, D<->H, N<->N
  # Example: RYWSKMBVDH -> DHBVKMSWRY
  echo "$1" | tr 'ACGTRYWSKMBVDHNacgtrywskmbvdhn' 'TGCAYRWSMKVBHDNtgcayrwsmkvhdbn' | rev
}

# Print key tool versions
version_log() {
  { command -v fastp >/dev/null 2>&1 && fastp --version; } || true
  { command -v cutadapt >/dev/null 2>&1 && cutadapt --version; } || true
  { command -v vsearch >/dev/null 2>&1 && vsearch --version | head -n1; } || true
  { command -v seqkit >/dev/null 2>&1 && seqkit version; } || true
  { command -v blastn >/dev/null 2>&1 && blastn -version | head -n1; } || true
}

# Usage
if [ $# -ne 3 ]; then
  echo "Usage: $0 <input_dir> <output_dir> <csv_file>"
  exit 1
fi
input_dir="$1"
output_dir="$2"
csv_file="$3"

# ---------------------------
# Dependency checks
# ---------------------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "Error: $1 not found"; exit 1; }; }
need fastp; need cutadapt; need vsearch; need python3; need seqkit

# ---------------------------
# Setup
# ---------------------------
mkdir -p "$output_dir"/{qc,trimmed,merged,markers/{tef,ssu,lsu,its},taxonomy,final_results}

version_log || true

echo "=== Multiplex Marker Extraction Pipeline Starting ==="
echo "Input: $input_dir"
echo "Output: $output_dir"
echo "CSV: $csv_file"

# ---------------------------
# Parse samples from CSV (primers are fixed for multiplex)
# ---------------------------
PARSER="$SCRIPT_DIR/parse_csv.py"
declare -a samples
while IFS=$'\t' read -r sample fwd rev fwdrc revrc; do
  [[ -z "$sample" ]] && continue
  samples+=("$sample")
done < <(python3 "$PARSER" "$csv_file")

echo "Samples: ${samples[*]}"

# ---------------------------
# Define markers and primers
# ---------------------------
declare -A fwd_primers rev_primers
fwd_primers=(
  ["its"]="GCATCGATGAAGAACGCAGC"
  ["ssu"]="CCGTGGTAATTCTAGAGCTAATACATGC"
  ["lsu"]="AGCGCACAAGTAGAGTGATCG"
  ["tef"]="TGGTACAAGGGYTGGGAGAAGG"
)
rev_primers=(
  ["its"]="YTTTTCCTCCGCTTATTGATATGC"
  ["ssu"]="RTCGGGATTGGGTAATTTGCGC"
  ["lsu"]="GGTCCGTGTTTCAAGACGG"
  ["tef"]="GCTGCTCGTGGTGCATYTCV"
)
markers=("tef" "ssu" "lsu" "its")

# ---------------------------
# Step 1: QC with fastp
# ---------------------------
echo "Step 1: QC (fastp)"
for sample in "${samples[@]}"; do
  r1=$(find "$input_dir" -type f \( -name "*${sample}*_R1*.fastq.gz" -o -name "*${sample}*_1*.fastq.gz" \) | grep -E '(_R1|_1)[^/]*\.fastq\.gz$' | head -1 || true)
  r2=$(find "$input_dir" -type f \( -name "*${sample}*_R2*.fastq.gz" -o -name "*${sample}*_2*.fastq.gz" \) | grep -E '(_R2|_2)[^/]*\.fastq\.gz$' | head -1 || true)
  if [[ -z "$r1" || -z "$r2" ]]; then echo "Warning: Files not found for $sample, skipping"; continue; fi

  fastp \
    -i "$r1" -I "$r2" \
    -o "$output_dir/qc/${sample}_R1_qc.fastq.gz" \
    -O "$output_dir/qc/${sample}_R2_qc.fastq.gz" \
    --detect_adapter_for_pe \
    --correction \
    --trim_poly_g \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --length_required 15 \
    --n_base_limit 5 \
    --thread $THREADS \
    --json "$output_dir/qc/${sample}_fastp.json" \
    --html "$output_dir/qc/${sample}_fastp.html" \
    >/dev/null 2>&1 || echo "Warning: fastp failed for $sample, skipping"
done

if command -v multiqc >/dev/null 2>&1; then
  multiqc "$output_dir/qc" -o "$output_dir/qc" -n multiqc_report.html --title "Multiplex Marker Pipeline QC Report" || echo "Warning: MultiQC failed"
fi

# ---------------------------
# Step 2: Marker-specific primer trimming (cutadapt)
# ---------------------------
echo "Step 2: Marker-specific primer trimming (cutadapt)"
for marker in "${markers[@]}"; do
  echo "Trimming for $marker"
  mkdir -p "$output_dir/trimmed/$marker"
  for sample in "${samples[@]}"; do
    [[ -f "$output_dir/qc/${sample}_R1_qc.fastq.gz" ]] || continue

    # Compute reverse complements for paired-end trimming
    rev_rc=$(reverse_complement "${rev_primers[$marker]}")
    fwd_rc=$(reverse_complement "${fwd_primers[$marker]}")

    # Comprehensive primer trimming: both 5' and 3' adapters
    # R1: -g FWD (5' end) and -a REV_RC (3' end, read-through)
    # R2: -G REV (5' end) and -A FWD_RC (3' end, read-through)
    # This handles both normal amplicons and short amplicons with primer read-through
    cutadapt \
      -g "${fwd_primers[$marker]}" \
      -a "${rev_rc}" \
      -G "${rev_primers[$marker]}" \
      -A "${fwd_rc}" \
      -o "$output_dir/trimmed/$marker/${sample}_R1_trimmed.fastq.gz" \
      -p "$output_dir/trimmed/$marker/${sample}_R2_trimmed.fastq.gz" \
      "$output_dir/qc/${sample}_R1_qc.fastq.gz" \
      "$output_dir/qc/${sample}_R2_qc.fastq.gz" \
      --discard-untrimmed \
      --minimum-length 10 \
      --error-rate 0.15 \
      --overlap 8 \
      --cores $THREADS \
      >"$output_dir/trimmed/$marker/${sample}_cutadapt.log" 2>&1 || echo "Warning: cutadapt failed for $sample $marker"
  done
done

# ---------------------------
# Step 3: Merge pairs and filter per marker (vsearch)
# ---------------------------
echo "Step 3: Merge pairs and filter per marker (vsearch)"
for marker in "${markers[@]}"; do
  echo "Merging for $marker"
  mkdir -p "$output_dir/merged/$marker"
  for sample in "${samples[@]}"; do
    [[ -f "$output_dir/trimmed/$marker/${sample}_R1_trimmed.fastq.gz" ]] || continue

    vsearch \
      --fastq_mergepairs "$output_dir/trimmed/$marker/${sample}_R1_trimmed.fastq.gz" \
      --reverse "$output_dir/trimmed/$marker/${sample}_R2_trimmed.fastq.gz" \
      --fastqout "$output_dir/merged/$marker/${sample}_merged.fastq" \
      --fastq_maxee 5 \
      --fastq_minmergelen 50 \
      --fastq_maxdiffs 30 \
      --fastq_maxdiffpct 100 \
      --threads $THREADS \
      || echo "Warning: vsearch merge failed for $sample $marker"
  done
done

# ---------------------------
# Step 4: Convert to FASTA per marker
# ---------------------------
echo "Step 4: Convert to FASTA per marker (quality filter)"
for marker in "${markers[@]}"; do
  for sample in "${samples[@]}"; do
    [[ -f "$output_dir/merged/$marker/${sample}_merged.fastq" ]] || continue
    vsearch \
      --fastq_filter "$output_dir/merged/$marker/${sample}_merged.fastq" \
      --fastq_maxee $MAX_EE \
      --fastq_minlen $MIN_LEN \
      --fastaout "$output_dir/markers/$marker/${sample}_$marker.fasta"
  done
done

# Prepare consensus summary table
CONSENSUS_SUMMARY="$output_dir/final_results/consensus_sequences.tsv"
echo -e "Sample\tMarker\tRank\tSize\tSeqID\tQueryID\tSequence" > "$CONSENSUS_SUMMARY"

# ---------------------------
# Step 5: Per-sample dereplication, chimera removal, sort by abundance, extract top-N
# ---------------------------
echo "Step 5: Per-sample dereplication, chimera removal, and top-N extraction"

for marker in "${markers[@]}"; do
  marker_dir="$output_dir/markers/$marker"
  mkdir -p "$marker_dir"

  for sample in "${samples[@]}"; do
    sample_fasta="$marker_dir/${sample}_${marker}.fasta"
    if [[ ! -s "$sample_fasta" ]]; then continue; fi

    # Enforce minimum length prior to dereplication
    sample_input="$marker_dir/${sample}_lenfilt.fasta"
    seqkit seq -m "$DEREP_MIN_LEN" "$sample_fasta" > "$sample_input" || true
    if ! grep -q '^>' "$sample_input" 2>/dev/null; then
      echo "Info: $marker $sample has no sequences >= ${DEREP_MIN_LEN} bp after filtering; skipping"
      continue
    fi

    # Dereplicate full-length sequences and track abundance
    vsearch --fastx_uniques "$sample_input" \
      --fastaout "$marker_dir/${sample}_uniques.fasta" \
      --sizeout \
      --relabel Uniq

    # Chimera removal on dereplicated sequences (uses size annotations)
    vsearch --uchime3_denovo "$marker_dir/${sample}_uniques.fasta" \
      --nonchimeras "$marker_dir/${sample}_nochim.fasta" \
      --sizein \
      --threads $THREADS || true

    # Choose input for sorting: prefer non-chimeras if available, otherwise use uniques
    input_for_sort="$marker_dir/${sample}_uniques.fasta"
    if [[ -s "$marker_dir/${sample}_nochim.fasta" ]]; then
      input_for_sort="$marker_dir/${sample}_nochim.fasta"
    fi

    # Sort by size (abundance)
    vsearch --sortbysize "$input_for_sort" \
      --output "$marker_dir/${sample}_sorted.fasta" \
      --sizein --sizeout

    # Take top N most abundant sequences as consensuses
    seqkit head -n "$TOP_N_ASVS" "$marker_dir/${sample}_sorted.fasta" > "$marker_dir/${sample}_consensus.fasta" || true

    # Fallback: if empty for any reason, take the very first record
    if ! grep -q '^>' "$marker_dir/${sample}_consensus.fasta" 2>/dev/null; then
      awk 'BEGIN{RS=">"; ORS=""} NR==2{print ">"$0}' "$marker_dir/${sample}_sorted.fasta" > "$marker_dir/${sample}_consensus.fasta" || true
    fi
  done

  # Combine consensuses for taxonomy
  : > "$marker_dir/consensuses.fasta"
  for sample in "${samples[@]}"; do
    if [[ -f "$marker_dir/${sample}_consensus.fasta" ]]; then
      # Ensure sample name in headers for downstream tracking
      awk -v s="$sample" 'BEGIN{OFS=""} /^>/{print ">" s "|" substr($0,2); next} {print}' "$marker_dir/${sample}_consensus.fasta" \
        >> "$marker_dir/consensuses.fasta"
    fi
  done

  if [[ -s "$marker_dir/consensuses.fasta" ]]; then
    seqkit fx2tab -n -s "$marker_dir/consensuses.fasta" \
      | awk -v m="$marker" 'BEGIN{FS="\t"; OFS="\t"}
        {
          query=$1;
          seq=$2;
          split(query, parts, "|");
          sample=parts[1];
          rest=parts[2];
          size="";
          if (match(rest, /;size=([0-9]+)/, arr)) {
            size=arr[1];
          }
          seqid=rest;
          sub(/;size=.*/, "", seqid);
          rank[sample] += 1;
          print sample, m, rank[sample], size, seqid, query, seq;
        }' >> "$CONSENSUS_SUMMARY"
  fi
done

# ---------------------------
# Step 6: BLAST taxonomy assignment (NCBI nt)
# ---------------------------
echo "Step 6: BLAST taxonomy assignment (NCBI nt)"

mkdir -p "$output_dir/taxonomy"

# Aggregated BLAST top hits and taxon IDs for all markers
BLAST_TOP_ALL="$output_dir/taxonomy/blast_top_all.tsv"
TAXIDS_ALL="$output_dir/taxonomy/taxids_all.txt"
: > "$BLAST_TOP_ALL"
: > "$TAXIDS_ALL"

# Ensure BLAST can locate taxonomy files (taxdb.bti/btd or taxonomy4blast.sqlite3)
TAXDB_DIR=$(dirname "$NTDB")
# Prepend DB dir to BLASTDB search path so blastn finds taxonomy files alongside the DB
if [[ -n "${BLASTDB:-}" ]]; then
  export BLASTDB="$TAXDB_DIR:$BLASTDB"
else
  export BLASTDB="$TAXDB_DIR"
fi
echo "BLASTDB set to: $BLASTDB"

for marker in "${markers[@]}"; do
  mdir="$output_dir/markers/$marker"
  CENTROIDS_FASTA="$output_dir/taxonomy/${marker}_centroids.fasta"
  cp "$mdir/consensuses.fasta" "$CENTROIDS_FASTA" 2>/dev/null || > "$CENTROIDS_FASTA"

  if [[ -s "$CENTROIDS_FASTA" ]]; then
    BLAST_OUT="$output_dir/taxonomy/${marker}_blast_results.tsv"
    # Determine if taxonomy files are available; if so, enable -taxids, else fall back gracefully
    TAXIDS_OPT=""
    if [[ -f "$TAXDB_DIR/taxdb.bti" && -f "$TAXDB_DIR/taxdb.btd" ]] || [[ -f "$TAXDB_DIR/taxonomy4blast.sqlite3" ]]; then
      TAXIDS_OPT=( -taxids "$BLAST_TAXIDS" )
      echo "Running blastn for $marker against: $NTDB (task=$BLAST_TASK; taxids=$BLAST_TAXIDS)"
    else
      echo "Warning: taxonomy files not found in $TAXDB_DIR; running BLAST without -taxids filtering for $marker"
      TAXIDS_OPT=()
      echo "Running blastn for $marker against: $NTDB (task=$BLAST_TASK; no taxid filter)"
    fi

    # Attempt BLAST with taxonomic filter if available; on failure, retry without -taxids
    if ! blastn -task "$BLAST_TASK" "${TAXIDS_OPT[@]}" -query "$CENTROIDS_FASTA" -db "$NTDB" -out "$BLAST_OUT" \
      -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle' \
      -evalue 1e-20 -perc_identity 80 -qcov_hsp_perc 80 -max_target_seqs 10 -num_threads $THREADS; then
      if [[ -n "${TAXIDS_OPT[*]}" ]]; then
        echo "Warning: blastn with -taxids failed for $marker; retrying without taxonomic filtering"
        blastn -task "$BLAST_TASK" -query "$CENTROIDS_FASTA" -db "$NTDB" -out "$BLAST_OUT" \
          -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle' \
          -evalue 1e-20 -perc_identity 80 -qcov_hsp_perc 80 -max_target_seqs 10 -num_threads $THREADS || echo "Warning: blastn failed for $marker (no-taxids)"
      else
        echo "Warning: blastn failed for $marker"
      fi
    fi

    if [[ -f "$BLAST_OUT" && -s "$BLAST_OUT" ]]; then
      BLAST_TOP="$output_dir/taxonomy/${marker}_blast_top.tsv"
      sort -k1,1 -k12,12nr "$BLAST_OUT" | awk '!seen[$1]++' > "$BLAST_TOP"

      # Aggregate BLAST top hit rows and taxids for global processing
      cat "$BLAST_TOP" >> "$BLAST_TOP_ALL"
      cut -f13 "$BLAST_TOP" | tr ';' '\n' | grep -v '^$' >> "$TAXIDS_ALL"
    else
      echo "BLAST produced no results for $marker; skipping taxonomy mapping"
    fi
  else
    echo "No centroids for $marker; skipping BLAST"
  fi
done

# ---------------------------
# Step 6b: Run taxonkit once and merge consensus, BLAST, and taxonomy summaries
# ---------------------------
FINAL_TABLE="$output_dir/final_results/barcode_summary.tsv"

# Prepare taxonomy mapping for all taxids (if taxonkit is available)
TAXONKIT_REFORMAT_ALL="$output_dir/taxonomy/taxids_reformat_all.tsv"
if command -v taxonkit >/dev/null 2>&1; then
  sort -u "$TAXIDS_ALL" > "$output_dir/taxonomy/taxids_all_unique.txt" || true
  if [[ -s "$output_dir/taxonomy/taxids_all_unique.txt" ]]; then
  taxonkit lineage --data-dir "$TAXONKIT_DB" "$output_dir/taxonomy/taxids_all_unique.txt" > "$output_dir/taxonomy/taxids_lineage_all.tsv" || true
  taxonkit reformat -i 2 -f '{k};{p};{c};{o};{f};{g};{s}' --data-dir "$TAXONKIT_DB" "$output_dir/taxonomy/taxids_lineage_all.tsv" > "$TAXONKIT_REFORMAT_ALL" || true
  fi
fi

python3 - "$CONSENSUS_SUMMARY" "$BLAST_TOP_ALL" "$TAXONKIT_REFORMAT_ALL" "$FINAL_TABLE" <<'PY'
import csv
import sys
from pathlib import Path

cons_path = Path(sys.argv[1])
blast_path = Path(sys.argv[2])
tax_path = Path(sys.argv[3])
final_path = Path(sys.argv[4])

headers = [
  "Sample", "Marker", "Rank", "Size", "SeqID", "QueryID", "Sequence",
  "FoundMarkerBLAST", "SubjectID", "SubjectTitle", "Pident", "AlignLength",
  "Bitscore", "Evalue", "TaxID", "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species"
]

def read_consensus(path: Path):
  if not path.exists() or path.stat().st_size == 0:
    return []
  with path.open() as handle:
    return list(csv.DictReader(handle, delimiter='\t'))

def read_blast_top(path: Path):
  # BLAST top file has no header; columns per outfmt:
  # 0 qseqid,1 sseqid,2 pident,3 length,4 mismatch,5 gapopen,6 qstart,7 qend,8 sstart,9 send,10 evalue,11 bitscore,12 staxids,13 stitle
  rows = []
  if not path.exists() or path.stat().st_size == 0:
    return rows
  with path.open() as handle:
    for line in handle:
      parts = line.rstrip('\n').split('\t')
      if len(parts) < 14:
        continue
      rows.append({
        'Query': parts[0],
        'SubjectID': parts[1],
        'Pident': parts[2],
        'AlignLength': parts[3],
        'Evalue': parts[10],
        'Bitscore': parts[11],
        'TaxID': parts[12].split(';')[0] if parts[12] else '',
        'SubjectTitle': parts[13]
      })
  return rows

def read_taxonkit_reformat(path: Path):
  # Expect: taxid \t <name?> \t k;p;c;o;f;g;s
  mapping = {}
  if not path.exists() or path.stat().st_size == 0:
    return mapping
  with path.open() as handle:
    for row in csv.reader(handle, delimiter='\t'):
      if not row:
        continue
      taxid = row[0]
      ranks = row[2] if len(row) > 2 else ''
      mapping[taxid] = ranks.split(';') if ranks else ['']*7
  return mapping

consensus_rows = read_consensus(cons_path)
blast_rows = read_blast_top(blast_path)
tax_map = read_taxonkit_reformat(tax_path)

blast_lookup = {row['Query']: row for row in blast_rows}

final_path.parent.mkdir(parents=True, exist_ok=True)
with final_path.open('w', newline='') as handle:
  writer = csv.writer(handle, delimiter='\t')
  writer.writerow(headers)

  for row in consensus_rows:
    query = row.get("QueryID", "")
    blast = blast_lookup.get(query, {})
    ranks = tax_map.get(blast.get('TaxID',''), ['']*7)
    title = blast.get('SubjectTitle', '')
    lower = title.lower()
    found = ''
    if 'internal transcribed spacer' in lower or 'its' in lower:
      found = 'ITS'
    elif 'small subunit' in lower or '18s' in lower or 'ssu' in lower:
      found = 'SSU'
    elif 'large subunit' in lower or '28s' in lower or 'lsu' in lower:
      found = 'LSU'
    elif 'elongation factor' in lower or 'tef' in lower:
      found = 'TEF'

    writer.writerow([
      row.get("Sample", ""),
      row.get("Marker", ""),
      row.get("Rank", ""),
      row.get("Size", ""),
      row.get("SeqID", ""),
      query,
      row.get("Sequence", ""),
      found,
      blast.get("SubjectID", ""),
      title,
      blast.get("Pident", ""),
      blast.get("AlignLength", ""),
      blast.get("Bitscore", ""),
      blast.get("Evalue", ""),
      blast.get("TaxID", ""),
      ranks[0] if len(ranks)>0 else '',
      ranks[1] if len(ranks)>1 else '',
      ranks[2] if len(ranks)>2 else '',
      ranks[3] if len(ranks)>3 else '',
      ranks[4] if len(ranks)>4 else '',
      ranks[5] if len(ranks)>5 else '',
      ranks[6] if len(ranks)>6 else '',
    ])
PY

echo "Final barcode summary table generated: $FINAL_TABLE"

# ---------------------------
# Step 7: Generate pipeline statistics report
# ---------------------------
generate_report() {
  if ! command -v seqkit >/dev/null 2>&1; then
    echo "Warning: seqkit not found, skipping pipeline stats report"
    return
  fi

  REPORT="$output_dir/final_results/pipeline_stats.tsv"
  mkdir -p "$output_dir/final_results"

  echo -e "Sample\tStep\tNum_Seqs\tMin_Len\tAvg_Len\tMax_Len" > "$REPORT"

  for sample in "${samples[@]}"; do
    # After QC
    if [[ -f "$output_dir/qc/${sample}_R1_qc.fastq.gz" ]]; then
      stats=$(seqkit stats "$output_dir/qc/${sample}_R1_qc.fastq.gz" --tabular | tail -n1 | awk '{print $4"\t"$6"\t"$7"\t"$8}')
      echo -e "$sample\tQC\t$stats" >> "$REPORT"
    fi

    # After trimming
    if [[ -f "$output_dir/trimmed/${sample}_R1_trimmed.fastq.gz" ]]; then
      stats=$(seqkit stats "$output_dir/trimmed/${sample}_R1_trimmed.fastq.gz" --tabular | tail -n1 | awk '{print $4"\t"$6"\t"$7"\t"$8}')
      echo -e "$sample\tTrimmed\t$stats" >> "$REPORT"
    fi

    # After merging
    if [[ -f "$output_dir/merged/${sample}_merged.fasta" ]]; then
      stats=$(seqkit stats "$output_dir/merged/${sample}_merged.fasta" --tabular | tail -n1 | awk '{print $4"\t"$6"\t"$7"\t"$8}')
      echo -e "$sample\tMerged\t$stats" >> "$REPORT"
    fi

    # Per marker
    for marker in "${markers[@]}"; do
      if [[ -f "$output_dir/markers/$marker/${sample}_${marker}.fasta" ]]; then
        stats=$(seqkit stats "$output_dir/markers/$marker/${sample}_${marker}.fasta" --tabular | tail -n1 | awk '{print $4"\t"$6"\t"$7"\t"$8}')
        echo -e "$sample\t${marker^^}\t$stats" >> "$REPORT"
      fi
      if [[ -f "$output_dir/markers/$marker/${sample}_filtered.fasta" ]]; then
        stats=$(seqkit stats "$output_dir/markers/$marker/${sample}_filtered.fasta" --tabular | tail -n1 | awk '{print $4"\t"$6"\t"$7"\t"$8}')
        echo -e "$sample\t${marker^^}_Filtered\t$stats" >> "$REPORT"
      fi
    done
  done

  # Total consensuses per marker
  for marker in "${markers[@]}"; do
    if [[ -f "$output_dir/markers/$marker/consensuses.fasta" ]]; then
      stats=$(seqkit stats "$output_dir/markers/$marker/consensuses.fasta" --tabular | tail -n1 | awk '{print $4"\t"$6"\t"$7"\t"$8}')
      echo -e "Total\t${marker^^}_Consensuses\t$stats" >> "$REPORT"
    fi
  done

  echo "Pipeline statistics report generated: $REPORT"
}

generate_report

echo "Done. Check taxonomy/ for BLAST results, and markers/ for per-marker consensuses."
