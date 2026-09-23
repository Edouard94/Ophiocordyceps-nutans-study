#!/usr/bin/env bash
set -euo pipefail

#NOTE: RERUN PHYLO WITH FORQUIGNOII (USE AS OG) / MAKE A SEPERATE ITS TREE?

# Multigene phylogeny pipeline driven by a CSV of sequences
# - Extracts per-gene amplicon FASTAs from CSV (including RPB1)
# - Generates consensus sequences for samples with multiple amplicons per gene
# - Uses existing voucher-backed reference sequences from references/ folder
# - Aligns, trims, concatenates (AMAS), and runs IQ-TREE2
# - Optionally computes platform agreement heatmaps

# Usage:
#   bash multigene_phylo.sh -i /path/to/All_seq_data.csv [-o output_dir] [-t THREADS]
#                           [--trim MODE] [--min-bp-per-gene N] [--codon-split]
#                           [--adjust-direction] [--qc] [--qc-max-n FLOAT]
#                           [--qc-max-mismatch FLOAT] [--qc-min-sites INT] [--consensus]
#
# Required:
#   -i, --input FILE       Input CSV of sequences (e.g., All_seq_data.csv)
#
# General options:
#   -o, --outdir DIR       Output directory (default: .)
#   -t, --threads N|AUTO   Thread count for MAFFT/IQ-TREE (default: AUTO)
#
# Trimming and per-sequence length filter:
#   --trim MODE            Trimming mode per gene alignment before concatenation.
#                          MODE in {none,gentle,moderate}. Default: none
#                          gentle  : trimAl -gt 0.9 (drop columns with >90% gaps)
#                          moderate: trimAl -gt 0.8 -resoverlap 0.75 -seqoverlap 70
#   --min-bp-per-gene N    After trimming, drop per-gene sequences with <N non-gap A/C/G/T sites.
#                          Default: 0 (no per-sequence filtering)
#
# Codon partitions (for protein-coding genes):
#   --codon-split          Split TEF, RPB1, RPB2 into 3 codon position partitions using
#                          RAxML step syntax (pos1 start-end\3, pos2 start+1-end\3, pos3 start+2-end\3).
#                          Assumes alignments are in-frame at partition start (no offset).
#
# Orientation correction (MAFFT):
#   --adjust-direction     Use MAFFT --adjustdirection to auto reverse-complement sequences
#                          when needed. Reversed headers are logged to logs/reversed_<GENE>.txt
#                          and MAFFT markers (_R_ / "reverse complemented") are stripped.
#
# Alignment quality control (QC) filtering:
#   --qc                   Enable post-trim QC filtering per gene alignment. Sequences exceeding
#                          thresholds are labeled REMOVE in logs/qc_<GENE>.tsv and excluded from
#                          downstream concatenation for that gene.
#   --qc-max-n FLOAT       Max allowed proportion of N (ambiguous bases) in SitesConsidered.
#                          Default: 0.5
#   --qc-max-mismatch FLOAT
#                          Max allowed proportion of mismatches to majority-rule consensus across
#                          SitesConsidered. Default: 0.30
#   --qc-min-sites INT     Minimum number of sites to consider for QC measurements (non-gap, A/C/G/T/N).
#                          Default: 100
#
# Consensus generation:
#   --consensus            Generate consensus sequences for samples with multiple amplicons per gene.
#                          Uses MAFFT for alignment and EMBOSS cons for consensus. Default: disabled
#
# Notes on QC behavior:
#   - QC is applied to the file actually used for concatenation (aligned or trimmed).
#   - Sequences marked REMOVE are dropped for that specific gene; the taxon can remain via other genes.
#   - In the concatenated NEXUS, dropped partitions are represented as missing data ("?") for those taxa.
#   - QC decisions per sequence are recorded in logs/qc_<GENE>.tsv with columns:
#       Header, SitesConsidered, N_Prop, Mismatch_Prop, Decision
#
# Example runs:
#   # Gentle trimming + codon partitions (default orientation/QC off)
#   nohup conda run -n base bash multigene_phylo.sh \
#     -i All_seq_data_203filt.csv --trim gentle --codon-split \
#     -o results_last -t AUTO > multigene_pipeline.log 2>&1 &
#
#   # Orientation correction and QC enabled, slightly stricter mismatch threshold
#   nohup conda run -n base bash multigene_phylo.sh \
#     -i All_seq_data_203filt.csv --trim gentle --adjust-direction --qc \
#     --qc-max-mismatch 0.25 --min-bp-per-gene 200 \
#     -o results_last -t AUTO > multigene_pipeline.log 2>&1 &
#
#   # Run without trimming but with QC (keeps full alignments, filters only by QC)
#   bash multigene_phylo.sh -i All_seq_data.csv --trim none --qc -o results_last
#
#   # Full pipeline with consensus, orientation, and QC
#   nohup conda run -n base bash multigene_phylo.sh \
#     -i All_seq_data_203filt.csv --consensus --trim gentle --adjust-direction --qc \
#     --qc-max-n 0.5 --qc-max-mismatch 0.30 --qc-min-sites 100 \
#     -o results_last -t AUTO > multigene_pipeline.log 2>&1 &
#
# Monitoring:
#   tail -f multigene_pipeline.log

CSV_INPUT=""
OUTDIR="."
THREADS="AUTO"
TRIM_MODE="none"
MIN_BP_PER_GENE=0
CODON_SPLIT=0
# New options: orientation adjustment and QC filtering
ADJUST_DIRECTION=0
QC_ENABLE=0
QC_MAX_N=0.5
QC_MAX_MISMATCH=0.30
QC_MIN_SITES=100
CONSENSUS_ENABLE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      CSV_INPUT="$2"
      shift 2
      ;;
    -o|--outdir)
      OUTDIR="$2"
      shift 2
      ;;
    -t|--threads)
      THREADS="$2"
      shift 2
      ;;
    --trim)
      TRIM_MODE="$2"
      shift 2
      ;;
    --min-bp-per-gene)
      MIN_BP_PER_GENE="$2"
      shift 2
      ;;
    --codon-split)
      CODON_SPLIT=1
      shift 1
      ;;
    --adjust-direction)
      ADJUST_DIRECTION=1
      shift 1
      ;;
    --qc)
      QC_ENABLE=1
      shift 1
      ;;
    --qc-max-n)
      QC_MAX_N="$2"
      shift 2
      ;;
    --qc-max-mismatch)
      QC_MAX_MISMATCH="$2"
      shift 2
      ;;
    --qc-min-sites)
      QC_MIN_SITES="$2"
      shift 2
      ;;
    --consensus)
      CONSENSUS_ENABLE=1
      shift 1
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$CSV_INPUT" ]]; then
  echo "Error: please provide -i path/to/CSV (e.g., All_seq_data.csv or All_input.csv)" >&2
  exit 1
fi

# Check dependencies
need() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Error: Missing dependency: $1" >&2
    exit 1
  }
}

for bin in awk sed grep cut sort uniq mafft cons; do
  need "$bin"
done
# Note: iqtree2 and AMAS are checked later before use

# If trimming requested, ensure trimAl is available and TRIM_MODE is valid
if [[ "$TRIM_MODE" != "none" ]]; then
  need trimal
  case "$TRIM_MODE" in
    gentle|moderate) : ;; 
    *) echo "Error: --trim must be one of {none,gentle,moderate}" >&2; exit 1;;
  esac
  if ! [[ "$MIN_BP_PER_GENE" =~ ^[0-9]+$ ]]; then
    echo "Error: --min-bp-per-gene must be a non-negative integer" >&2
    exit 1
  fi
fi

# Layout
WORKDIR=$(pwd)
mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)  # Convert to absolute path

ALIGN_DIR="$OUTDIR/phylogeny/alignments"
TRIM_DIR="$OUTDIR/phylogeny/trimmed"
TREE_DIR="$OUTDIR/phylogeny/trees"
PART_DIR="$OUTDIR/phylogeny/partitions"
REF_DIR="$WORKDIR/references_last"
# Fallback to references/ if references_last/ is not present
if [[ ! -d "$REF_DIR" ]]; then
  REF_DIR="$WORKDIR/references"
fi
AMPL_DIR="$OUTDIR/amplicons"
LOG_DIR="$OUTDIR/logs"
FIG_DIR="$OUTDIR/phylogeny/figures"
AMAS_OUT="$OUTDIR/phylogeny/amas_output"

mkdir -p "$ALIGN_DIR" "$TRIM_DIR" "$TREE_DIR" "$PART_DIR" "$AMPL_DIR" "$LOG_DIR" "$FIG_DIR" "$AMAS_OUT"

# Central pipeline log
PIPELOG="$LOG_DIR/pipeline.log"

# Map THREADS to MAFFT threads: AUTO => -1, numeric => that value
MAFFT_THR="-1"
if [[ "$THREADS" =~ ^[0-9]+$ ]]; then MAFFT_THR="$THREADS"; fi

GENES=(SSU LSU TEF RPB1 RPB2)  # include or exclude genes as needed

echo "Genes included in this run: ${GENES[*]}" >> "$PIPELOG"

# Build a regex from GENES for safe header stripping (only strip known gene tags)
# Example: (ITS|SSU|LSU|TEF|RPB1)
GENE_RE="$(printf "%s|" "${GENES[@]}")"
GENE_RE="(${GENE_RE%|})"

# 1) Extract amplicons from CSV
echo "Step 1: Extracting amplicons from CSV..."
python3 "$WORKDIR/scripts/csv_to_amplicons.py" "$CSV_INPUT" "$AMPL_DIR" \
  2>"$LOG_DIR/csv_to_amplicons.log"

# 1.5) Generate consensus amplicons per sample per gene
if [[ $CONSENSUS_ENABLE -eq 1 ]]; then
echo "Step 1.5: Generating consensus amplicons where multiple sequences exist per sample per gene..."
for gene in "${GENES[@]}"; do
  AMP_FILE="$AMPL_DIR/${gene}_amplicons.fasta"
  if [[ -s "$AMP_FILE" ]]; then
    TMP_DIR="$OUTDIR/tmp_consensus_${gene}"
    mkdir -p "$TMP_DIR"
    # Split by sample: headers are >Platform_Sample_GENE, extract Sample as everything between first _ and _GENE
    awk -v gene="$gene" '/^>/ { hdr = substr($0,2); sub(/_'"$gene"'$/, "", hdr); split(hdr, a, "_"); sample = a[2]; for(i=3; i<=length(a); i++) sample = sample "_" a[i]; file = "'$TMP_DIR'/" sample ".fasta"; print > file; next } { print >> file }' "$AMP_FILE"
    CONSENSUS_COUNT=0
    RENAME_COUNT=0
    for sample_fa in "$TMP_DIR"/*.fasta; do
      if [[ ! -f "$sample_fa" ]]; then continue; fi
      nseq=$(grep -c '^>' "$sample_fa")
      sample=$(basename "$sample_fa" .fasta)
      if (( nseq > 1 )); then
        # Align with MAFFT
        mafft --auto --thread "$MAFFT_THR" "$sample_fa" > "$sample_fa.aln" 2>/dev/null
        # Consensus with cons
        cons -sequence "$sample_fa.aln" -outseq "$sample_fa.cons" -name "${sample}_${gene}" 2>/dev/null
        if [[ -s "$sample_fa.cons" ]]; then
          CONSENSUS_COUNT=$((CONSENSUS_COUNT + 1))
        fi
      elif (( nseq == 1 )); then
        # Standardize singleton header to Sample_GENE (drop platform prefix)
        echo ">${sample}_${gene}" > "$sample_fa.cons"
        awk '!/^>/' "$sample_fa" >> "$sample_fa.cons"
        if [[ -s "$sample_fa.cons" ]]; then
          RENAME_COUNT=$((RENAME_COUNT + 1))
        fi
      fi
    done
    # Rebuild AMP_FILE
    > "$AMP_FILE.new"
    for sample_fa in "$TMP_DIR"/*.fasta; do
      if [[ ! -f "$sample_fa" ]]; then continue; fi
      if [[ -s "$sample_fa.cons" ]]; then
        cat "$sample_fa.cons" >> "$AMP_FILE.new"
      else
        # Fallback: include original if no .cons produced (should be rare)
        cat "$sample_fa" >> "$AMP_FILE.new"
      fi
    done
    mv "$AMP_FILE.new" "$AMP_FILE"
    rm -rf "$TMP_DIR"
    if (( CONSENSUS_COUNT > 0 || RENAME_COUNT > 0 )); then
      echo "  $gene: generated $CONSENSUS_COUNT consensuses; standardized $RENAME_COUNT singletons" >> "$PIPELOG"
    fi
  fi
done
echo "  Consensus generation complete."
fi

# 2) Use existing reference sequences from references/ folder
echo "Step 2: Using existing reference sequences from $REF_DIR..."
for gene in "${GENES[@]}"; do
  REF_FILE="$REF_DIR/${gene}_refs.fasta"
  if [[ -s "$REF_FILE" ]]; then
    echo "  Found references for $gene: $REF_FILE"
  else
    echo "  Warning: No reference file found for $gene at $REF_FILE" | tee -a "$PIPELOG"
  fi
done

# 3) Combine amplicons and references, align, trim
echo "Step 3: Combining amplicon and reference sequences..."
PROCESSED=()
CONCAT_FILES=()

for gene in "${GENES[@]}"; do
  REF_FILE="$REF_DIR/${gene}_refs.fasta"
  AMP_FILE="$AMPL_DIR/${gene}_amplicons.fasta"
  COMBINED="$ALIGN_DIR/${gene}_combined.fasta"
  
  > "$COMBINED"
  
  # Add references with standardized headers (Species_VoucherID only)
  if [[ -s "$REF_FILE" ]]; then
    # Remove gene marker and accession from headers: >Species_VoucherID_GENE_ACC -> >Species_VoucherID
    # Use dynamic gene regex and allow optional accession version (e.g., AB123456 or AB123456.1)
    # Only strip when header matches the exact _GENE_ACC pattern at the end
  # Allow underscore in accessions (e.g., NG_064858.1)
  # After stripping, normalize the voucher token (final token) by removing internal underscores
  # Example: >Ophiocordyceps_brunneaperitheciata_BCC_49312_TEF_XXXX -> >Ophiocordyceps_brunneaperitheciata_BCC49312
  sed -E "s/>(.+)_${GENE_RE}_[A-Za-z0-9_]+(\.[0-9]+)?$/>\1/" "$REF_FILE" \
    | awk 'BEGIN{OFS=""}
      /^>/ {
        s = substr($0,2);
        if (match(s, /_[^_]*$/)) {
          prefix = substr(s, 1, RSTART-1);
          suffix = substr(s, RSTART+1);
          gsub(/_/, "", suffix);
          print ">", prefix, "_", suffix;
        } else {
          print ">", s;
        }
        next
      }
      { print }' >> "$COMBINED"
  fi
  
  # Add amplicons with standardized headers (remove gene suffix if present)
  if [[ -s "$AMP_FILE" ]]; then
    # Remove gene marker from headers: >Sample_ID_GENE -> >Sample_ID, but only if the final token is a known gene
    # This avoids stripping sample suffixes like _A / _B
    sed -E "s/>(.+)_${GENE_RE}$/>\1/" "$AMP_FILE" >> "$COMBINED"
  fi

  # Sanitize headers and sequences
  if [[ -s "$COMBINED" ]]; then
    # Normalize headers: allow only A-Za-z0-9_ in IDs (convert everything else to _)
    # This removes characters like ':', '#', '-', '()' etc. except underscores
    awk 'BEGIN{OFS=""} 
      /^>/ { 
        s = substr($0,2); 
        gsub(/[^A-Za-z0-9_]/, "_", s); 
        gsub(/_+/, "_", s);              # collapse multiple underscores
        sub(/^_+/, "", s);               # trim leading underscores (if any)
        sub(/_+$/, "", s);               # trim trailing underscores
        # Fix missing underscore between species and voucher (e.g., Ophiocordyceps_nutans11BidoupNLU -> Ophiocordyceps_nutans_11BidoupNLU)
        # Apply conservatively only when header starts with a binomial (Genus_species)
        if (s ~ /^[A-Z][a-z]+_/) {
          last = s; sub(/.*_/, "", last);                    # last token after the final underscore
          if (match(last, /^([A-Za-z]+)([0-9].*)$/, a)) {      # letters immediately followed by digits
            base = substr(s, 1, length(s) - length(last));
            s = base a[1] "_" a[2];
          }
        }
        # Unify voucher patterns like Species_CODE_DIGITS -> Species_CODEDIGITS (e.g., OSC_71233 -> OSC71233)
        # Only apply when last two underscore-separated tokens are letters then digits (optionally alphanum at end)
        if (match(s, /(.*)_([A-Za-z]+)_([0-9]+[A-Za-z]*)$/, arr)) {
          s = arr[1] "_" arr[2] arr[3];
        }
        print ">", s; 
        next 
      } 
      { print }' "$COMBINED" > "$COMBINED.tmp" && mv "$COMBINED.tmp" "$COMBINED"
    # Replace non-standard nucleotides with N
    sed -i.bak '/^[^>]/s/[^ATCGNatcgn-]/N/g' "$COMBINED" && rm -f "$COMBINED.bak"

    # Log header counts and detect duplicate headers after sanitization
    HDR_COUNT=$(grep -c '^>' "$COMBINED" || true)
    DUP_FILE="$LOG_DIR/dup_headers_${gene}.txt"
    grep '^>' "$COMBINED" | sort | uniq -d > "$DUP_FILE" || true
    if [[ -s "$DUP_FILE" ]]; then
      echo "  Warning: Found $(wc -l < "$DUP_FILE") duplicate headers for $gene (see $DUP_FILE)" | tee -a "$PIPELOG"
    fi
    echo "  Headers for $gene after sanitization: $HDR_COUNT" >> "$PIPELOG"
  else
    echo "  Warning: No sequences to combine for $gene" | tee -a "$PIPELOG"
    continue
  fi

  if [[ ! -s "$COMBINED" ]]; then
    echo "  Warning: No sequences for $gene, skipping alignment" | tee -a "$LOG_DIR/align_${gene}.log"
    continue
  fi

  echo "  Aligning $gene with MAFFT..."
  # MAFFT options
  MAFFT_OPTS=(--auto --thread "$MAFFT_THR")
  if [[ $ADJUST_DIRECTION -eq 1 ]]; then
    MAFFT_OPTS+=(--adjustdirection)
  fi
  mafft "${MAFFT_OPTS[@]}" "$COMBINED" \
    > "$ALIGN_DIR/${gene}_aligned.fasta" \
    2>>"$LOG_DIR/align_${gene}.log"

  # If adjust-direction used, log sequences that were reverse-complemented and strip MAFFT markers from headers
  if [[ $ADJUST_DIRECTION -eq 1 && -s "$ALIGN_DIR/${gene}_aligned.fasta" ]]; then
    REV_LOG="$LOG_DIR/reversed_${gene}.txt"
    # MAFFT marks reversed sequences by appending '_R_' to names; also sometimes adds '(reverse complemented)'
    grep '^>' "$ALIGN_DIR/${gene}_aligned.fasta" | grep -E '(_R_|reverse complemented)' | sed 's/^>//' > "$REV_LOG" || true
    # Strip markers from headers in aligned file
    sed -E 's/^>(.*)_R_/>\1/; s/\s*\(reverse complemented\)//' "$ALIGN_DIR/${gene}_aligned.fasta" > "$ALIGN_DIR/${gene}_aligned.tmp" && mv "$ALIGN_DIR/${gene}_aligned.tmp" "$ALIGN_DIR/${gene}_aligned.fasta"
  fi
  
  USE_FILE="$ALIGN_DIR/${gene}_aligned.fasta"

  # Optional: per-gene trimming with trimAl
  if [[ "$TRIM_MODE" != "none" && -s "$ALIGN_DIR/${gene}_aligned.fasta" ]]; then
    echo "  Trimming $gene alignment with trimAl ($TRIM_MODE)..." | tee -a "$LOG_DIR/align_${gene}.log"
    TRIM_FILE="$TRIM_DIR/${gene}_trimmed.fasta"
    mkdir -p "$TRIM_DIR"
    if [[ "$TRIM_MODE" == "gentle" ]]; then
      trimal -in "$ALIGN_DIR/${gene}_aligned.fasta" -out "$TRIM_FILE" -gt 0.9 \
        2>>"$LOG_DIR/align_${gene}.log" || true
    else
      trimal -in "$ALIGN_DIR/${gene}_aligned.fasta" -out "$TRIM_FILE" -gt 0.8 -resoverlap 0.75 -seqoverlap 70 \
        2>>"$LOG_DIR/align_${gene}.log" || true
    fi

    # Optional: filter per-sequence by retained non-gap A/C/G/T length
    if [[ -s "$TRIM_FILE" && $MIN_BP_PER_GENE -gt 0 ]]; then
      echo "  Filtering $gene sequences with <${MIN_BP_PER_GENE} bp (A/C/G/T) after trimming..." | tee -a "$LOG_DIR/align_${gene}.log"
      awk -v MIN=$MIN_BP_PER_GENE 'BEGIN{FS="\t"}
        function flush(){
          if(hdr!=""){ if(len >= MIN) { print ">" hdr; print seq } }
          hdr=""; seq=""; len=0;
        }
        /^>/ { flush(); hdr=substr($0,2); next }
        { g=$0; gsub(/[^ACGTacgt]/, "", g); len+=length(g); seq=seq $0 }
        END{ flush() }' "$TRIM_FILE" > "$TRIM_FILE.filtered" && mv "$TRIM_FILE.filtered" "$TRIM_FILE"
    fi

    # If trimming produced a usable alignment, prefer it; otherwise fall back to untrimmed
    if [[ -s "$TRIM_FILE" ]]; then
      # Log counts before/after
      NSEQ_ALN=$(grep -c '^>' "$ALIGN_DIR/${gene}_aligned.fasta" || true)
      NSEQ_TRM=$(grep -c '^>' "$TRIM_FILE" || true)
      echo "  $gene: sequences before/after trimming = ${NSEQ_ALN}/${NSEQ_TRM}" >> "$PIPELOG"
      USE_FILE="$TRIM_FILE"
    else
      echo "  Warning: Trimming for $gene removed all data; using untrimmed alignment" | tee -a "$LOG_DIR/align_${gene}.log"
    fi
  fi

  # Optional: alignment QC filtering (consensus mismatch / N proportion)
  if [[ $QC_ENABLE -eq 1 && -s "$USE_FILE" ]]; then
    echo "  QC filtering $gene alignment (max N: $QC_MAX_N, max mismatch: $QC_MAX_MISMATCH, min sites: $QC_MIN_SITES)..." | tee -a "$LOG_DIR/align_${gene}.log"
    python3 "$WORKDIR/scripts/qc_alignments.py" \
      --input "$USE_FILE" \
      --output "$USE_FILE.qc" \
      --log "$LOG_DIR/qc_${gene}.tsv" \
      --max-n "$QC_MAX_N" \
      --max-mismatch "$QC_MAX_MISMATCH" \
      --min-sites "$QC_MIN_SITES" \
      >>"$LOG_DIR/align_${gene}.log" 2>&1 || true
    if [[ -s "$USE_FILE.qc" ]]; then
      NSEQ_BEFORE=$(grep -c '^>' "$USE_FILE" || true)
      NSEQ_AFTER=$(grep -c '^>' "$USE_FILE.qc" || true)
      mv "$USE_FILE.qc" "$USE_FILE"
      echo "  $gene: sequences before/after QC = ${NSEQ_BEFORE}/${NSEQ_AFTER}" >> "$PIPELOG"
    fi
  fi

  if [[ -s "$ALIGN_DIR/${gene}_aligned.fasta" ]]; then
    PROCESSED+=("$gene")
    # collect file for concatenation (trimmed if present)
    if [[ -s "$USE_FILE" ]]; then
      CONCAT_FILES+=("$USE_FILE")
    fi
  else
    echo "  Warning: Alignment failed for $gene, skipping" | tee -a "$LOG_DIR/align_${gene}.log"
  fi
done


  if (( ${#PROCESSED[@]} == 0 )); then
    echo "Error: No gene produced a valid alignment" >&2
    exit 1
  fi

echo "Successfully processed genes: ${PROCESSED[*]}"

# 4) Concatenate with AMAS
echo "Step 4: Concatenating alignments with AMAS..."

AMAS_CMD=""
if command -v AMAS.py >/dev/null 2>&1; then
  AMAS_CMD=("AMAS.py")
elif python3 -c "import amas.AMAS" 2>/dev/null; then
  AMAS_CMD=("python3" "-m" "amas.AMAS")
else
  echo "Error: AMAS not found. Install with 'pip install AMAS' or 'conda install -c bioconda amas'" >&2
  exit 1
fi

   "${AMAS_CMD[@]}" concat \
    -f fasta \
    -d dna \
    -i "${CONCAT_FILES[@]}" \
  -u nexus \
  -t "$AMAS_OUT/concatenated_alignment.nex" \
  -p "$AMAS_OUT/partitions.txt" \
  -e \
  -c 4 \
  > "$LOG_DIR/amas.log" 2>&1 || {
    echo "Error: AMAS concatenation failed" >&2
    exit 1
  }

echo "  Concatenated alignment written to $AMAS_OUT/concatenated_alignment.nex"

# 5) IQ-TREE2
echo "Step 5: Running IQ-TREE2 phylogenetic inference..."

if command -v iqtree >/dev/null 2>&1; then
  if [[ ! -s "$AMAS_OUT/concatenated_alignment.nex" ]]; then
    echo "Error: Missing concatenated alignment: $AMAS_OUT/concatenated_alignment.nex" >&2
    exit 1
  fi
  # Prepare partition source file (optionally codon-split) then convert to RAxML-style
  PART_SRC="$AMAS_OUT/partitions.txt"
  if [[ -s "$PART_SRC" && $CODON_SPLIT -eq 1 ]]; then
    echo "  Generating codon position partitions (TEF, RPB1, RPB2)..." | tee -a "$LOG_DIR/iqtree.log"
    CODON_FILE="$AMAS_OUT/partitions_codon.txt"
    > "$CODON_FILE"
    CODON_CREATED=0
    while IFS= read -r line; do
      [[ -z "$line" ]] && continue
      name_part="${line%%=*}"; range_part="${line#*=}"
      name_part="$(echo "$name_part" | sed -E 's/^[[:space:]]+|[[:space:]]+$//g')"
      range_part="$(echo "$range_part" | sed -E 's/^[[:space:]]+|[[:space:]]+$//g')"
      # Accept *_aligned or *_trimmed
      gene="$(echo "$name_part" | sed -E 's/^p[0-9]+_([A-Za-z0-9]+)_(aligned|trimmed)$/\1/')"
      start="${range_part%-*}"; end="${range_part#*-}"
      if [[ "$gene" =~ ^(TEF|RPB1|RPB2)$ ]]; then
        if [[ "$start" =~ ^[0-9]+$ && "$end" =~ ^[0-9]+$ && $start -le $end ]]; then
          length=$((end-start+1))
          if (( length % 3 != 0 )); then
            echo "# WARNING: $gene partition length $length not divisible by 3; codon split assumes frame at start" >> "$CODON_FILE"
          fi
          pos1_start=$start; pos2_start=$((start+1)); pos3_start=$((start+2))
          echo "${gene}_pos1 = ${pos1_start}-${end}\\3" >> "$CODON_FILE"
          echo "${gene}_pos2 = ${pos2_start}-${end}\\3" >> "$CODON_FILE"
          echo "${gene}_pos3 = ${pos3_start}-${end}\\3" >> "$CODON_FILE"
          CODON_CREATED=$((CODON_CREATED+3))
        else
          echo "# Skipping invalid range for $gene: $range_part" >> "$CODON_FILE"
          echo "$line" >> "$CODON_FILE"
        fi
      else
        echo "$line" >> "$CODON_FILE"
      fi
    done < "$PART_SRC"
    PART_SRC="$CODON_FILE"
    echo "  Codon partitions generated: $CODON_CREATED (check $CODON_FILE)" | tee -a "$LOG_DIR/iqtree.log"
  fi
  # Convert to RAxML-style for IQ-TREE
  if [[ -s "$PART_SRC" ]]; then
    IQ_PART="$AMAS_OUT/iqtree_partitions.txt"
    awk 'NF && $0 !~ /^#/ { gsub(/\r/,"\n"); split($0,a,"="); if (length(a)==2) {name=a[1]; gsub(/^[[:space:]]+|[[:space:]]+$/,"",name); rng=a[2]; gsub(/^[[:space:]]+|[[:space:]]+$/,"",rng); printf("DNA, %s = %s\n", name, rng);} }' "$PART_SRC" > "$IQ_PART"
    if [[ ! -s "$IQ_PART" ]]; then
      echo "Warning: Failed to convert partition file ($PART_SRC); will run unpartitioned." | tee -a "$LOG_DIR/iqtree.log"
    fi
  fi
  if [[ ! -s "$AMAS_OUT/partitions.txt" || ! -s "$IQ_PART" ]]; then
    echo "Warning: Partition file not found at $AMAS_OUT/partitions.txt; running unpartitioned model" | tee -a "$LOG_DIR/iqtree.log"
    set +e
    iqtree \
      -s "$AMAS_OUT/concatenated_alignment.nex" \
      -m MFP \
      -bb 1000 \
      -bnni \
      -alrt 1000 \
      -nt "$THREADS" \
      -pre "$TREE_DIR/multigene_partitioned" \
      > "$LOG_DIR/iqtree.log" 2>&1
    IQEXIT=$?
    set -e
    if [[ $IQEXIT -ne 0 ]]; then
      # If key outputs exist despite non-zero exit, warn and continue
      if [[ -s "$TREE_DIR/multigene_partitioned.treefile" && -s "$TREE_DIR/multigene_partitioned.contree" ]]; then
        echo "Warning: IQ-TREE exited with code $IQEXIT but tree outputs exist; continuing. See $LOG_DIR/iqtree.log and $TREE_DIR/multigene_partitioned.log" | tee -a "$LOG_DIR/iqtree.log"
      else
        echo "Error: IQ-TREE run failed (exit $IQEXIT) and no tree outputs found" >&2
        exit 1
      fi
    fi
  else
    set +e
    iqtree \
      -s "$AMAS_OUT/concatenated_alignment.nex" \
      -q "$IQ_PART" \
      -m MFP \
      -bb 1000 \
      -bnni \
      -alrt 1000 \
      -nt "$THREADS" \
      -pre "$TREE_DIR/multigene_partitioned" \
      > "$LOG_DIR/iqtree.log" 2>&1
    IQEXIT=$?
    set -e
    if [[ $IQEXIT -ne 0 ]]; then
      if [[ -s "$TREE_DIR/multigene_partitioned.treefile" && -s "$TREE_DIR/multigene_partitioned.contree" ]]; then
        echo "Warning: IQ-TREE exited with code $IQEXIT but tree outputs exist; continuing. See $LOG_DIR/iqtree.log and $TREE_DIR/multigene_partitioned.log" | tee -a "$LOG_DIR/iqtree.log"
      else
        echo "Error: IQ-TREE run failed (exit $IQEXIT) and no tree outputs found" >&2
        exit 1
      fi
    fi
  fi
  echo "  Tree inference complete: $TREE_DIR/multigene_partitioned.treefile"
else
  echo "Warning: iqtree not found in PATH. Skipping tree inference" >&2
fi

# 6) Platform agreement matrices and heatmaps
echo "Step 6: Generating platform comparison plots..."
python3 "$WORKDIR/scripts/compare_platforms.py" \
  --align-dir "$TRIM_DIR" \
  --outdir "$FIG_DIR" \
  > "$LOG_DIR/platform_compare.log" 2>&1 || true

# Summary
echo ""
echo "========================================="
echo "Multigene phylogeny pipeline completed!"
echo "========================================="
echo "Outputs:"
echo "  - Per-gene alignments:         $ALIGN_DIR/*_aligned.fasta"
echo "  - AMAS concatenated alignment: $AMAS_OUT/concatenated_alignment.nex"
echo "  - Partition file:              $AMAS_OUT/partitions.txt"
echo "  - IQ-TREE results (if run):    $TREE_DIR/multigene_partitioned.*"
echo "  - Platform comparison plots:   $FIG_DIR/"
echo "  - Logs:                        $LOG_DIR/"
echo ""
echo "Notes:"
echo "  - Reference sequences from: $REF_DIR/"
echo "  - Input CSV: $CSV_INPUT"
echo "  - Genes processed: ${PROCESSED[*]}"
echo ""
