#!/bin/bash

# ONT Amplicon Pipeline Script
# Usage: bash ONT_Amplicon_Pipeline.sh <CSV_FILE> <FASTQ_FILE>

# Check inputs
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <CSV_FILE> <FASTQ_FILE>"
    exit 1
fi

CSV_FILE="$1"
FASTQ_FILE="$2"

# Parameters
TAG_ERROR_RATE=0 # Maximum allowed error rate (if 0 <= E < 1), or absolute number of errors for full-length adapter match (if E is an integer >= 1).
QUALITY_THRESHOLD=7
MIN_SIZE=100
MAX_SIZE=400
SAMPLE_SIZE=500 # As too many reads can add noise to the clustering and polishing, we recommend subsampling the read data to 300–500 reads per sample (https://doi.org/10.1038/s41596-022-00682-x)
CORES=20
TRIM_ONT_BARCODES=false   # set to false if your reads are already barcode-trimmed
ONT_BARCODE_KIT=""       # e.g., "EXP-NBD196" or "EXP-NBD104"; required if TRIM_ONT_BARCODES=true
MIN_CONSENSUS_SUPPORT=10  # Minimum number of supporting reads for consensus to be kept
NTDB=${NTDB:-/data2/nt_NCBI_database/nt}  # Path to NCBI nt database
TAXONKIT_DB=${TAXONKIT_DB:-$HOME/taxdump}  # Path to TaxonKit taxonomy database
# BLAST performance tuning (match multiplex pipeline defaults)
BLAST_TASK=${BLAST_TASK:-megablast}
BLAST_TAXIDS=${BLAST_TAXIDS:-4751}  # fungi

# Output directory
OUT_DIR="ONT_Pipeline_Output"
mkdir -p "$OUT_DIR"

# Function to compute reverse complement (including wobble bases)
reverse_complement() {
    echo "$1" | tr 'ATCG RYWS KMBD HVN atcg ryws kmbd hvn' 'TAGC YRWS MKVH DBH tagc yrws mkv hdb' | rev
}

# Optional: trim/demultiplex ONT barcodes before tag demultiplexing
if $TRIM_ONT_BARCODES; then
    if [ -z "$ONT_BARCODE_KIT" ]; then
        echo "ERROR: TRIM_ONT_BARCODES is true but ONT_BARCODE_KIT is not set."
        exit 1
    fi
    echo "Trimming/demultiplexing ONT kit barcodes with guppy_barcoder ($ONT_BARCODE_KIT)..."
    DEMUX_DIR="$OUT_DIR/ont_demux"
    mkdir -p "$DEMUX_DIR"
    # This will create per-barcode subfolders and trim barcodes from reads
    guppy_barcoder \
        -i "$(dirname "$FASTQ_FILE")" \
        -s "$DEMUX_DIR" \
        --barcode_kits "$ONT_BARCODE_KIT" \
        --trim_barcodes \
        --require_barcodes_both_ends \
        --num_extra_bases_trim 100 \
        --recursive || true

    # Merge all demultiplexed/truncated reads back into a single cleaned FASTQ for tag demux,
    # or adjust your workflow to use per-barcode bins if desired.
    find "$DEMUX_DIR" -type f -name "*.fastq" -print0 | xargs -0 cat > "$OUT_DIR/ont_barcode_trimmed.fastq" 2>/dev/null || true
    FASTQ_FILE="$OUT_DIR/ont_barcode_trimmed.fastq"
fi

# Get unique samples
samples=$(awk -F, 'NR>1 {print $1}' "$CSV_FILE" | sort | uniq)

for sample in $samples; do
    echo "Processing sample: $sample"

    # Get tags for this sample (assume same for all markers in sample)
    fwd_tag=$(awk -F, '$1=="'$sample'" {print $4; exit}' "$CSV_FILE")
    rev_tag=$(awk -F, '$1=="'$sample'" {print $5; exit}' "$CSV_FILE")
    rev_tag_rc=$(reverse_complement "$rev_tag")
    fwd_tag_rc=$(reverse_complement "$fwd_tag")

    # Per-tag min_overlap = full tag length
    #len_fwd_tag=${#fwd_tag}
    #len_rev_tag=${#rev_tag}

    # Create sample directory
    SAMPLE_DIR="$OUT_DIR/$sample"
    mkdir -p "$SAMPLE_DIR/demultiplexed" "$SAMPLE_DIR/stats"

    # Tags are searched as linked adapters (-a ^ADAPTER1...ADAPTER2 -g ADAPTER1...ADAPTER2)
    # with -a, the adapters that are anchored become required
    # Cutadapt versions before 2.0 anchored the 5’ adapter within linked adapters automatically even if the initial ^ was not specified.
    # Require both adapters using ;required. Do not set min_overlap for anchored adapters.
    # The minimum overlap length cannot be set for anchored adapters as these always need to occur at full length.
    # Demultiplex: Forward orientation.
    cutadapt -a "^$fwd_tag...${rev_tag_rc}$" \
        --error-rate "$TAG_ERROR_RATE" --trimmed-only \
        --cores "$CORES" \
        -o "$SAMPLE_DIR/demultiplexed/${sample}_fwd_tags.fastq" \
        "$FASTQ_FILE" > "$SAMPLE_DIR/stats/${sample}_demultiplex_fwd.log" 2>&1 || true

    # Demultiplex: Reverse orientation
    cutadapt -a "^$rev_tag...${fwd_tag_rc}$" \
        --error-rate "$TAG_ERROR_RATE" --trimmed-only \
        --cores "$CORES" \
        -o "$SAMPLE_DIR/demultiplexed/${sample}_rev_tags.fastq" \
        "$FASTQ_FILE" > "$SAMPLE_DIR/stats/${sample}_demultiplex_rev.log" 2>&1 || true

    # Combine demultiplexed reads
    cat "$SAMPLE_DIR/demultiplexed/${sample}_fwd_tags.fastq" "$SAMPLE_DIR/demultiplexed/${sample}_rev_tags.fastq" > "$SAMPLE_DIR/${sample}_demultiplexed.fastq" 2>/dev/null || true

    # Get markers for this sample
    markers=$(awk -F, '$1=="'$sample'" {print $2}' "$CSV_FILE")

    for marker in $markers; do
        echo "  Processing marker: $marker for sample: $sample"

        # Get primers
        fwd_primer=$(awk -F, '$1=="'$sample'" && $2=="'$marker'" {print $6}' "$CSV_FILE")
        rev_primer=$(awk -F, '$1=="'$sample'" && $2=="'$marker'" {print $7}' "$CSV_FILE")

        # Calculate min overlaps (80% of forward primer length, 50% of reverse primer length)
        min_overlap_forward=$(( ${#fwd_primer} * 80 / 100 ))
        min_overlap_reverse=$(( ${#rev_primer} * 50 / 100 ))

        # Create marker directory
        MARKER_DIR="$SAMPLE_DIR/$marker"
        mkdir -p "$MARKER_DIR"

        # Trim primers: Forward orientation
        # -n 2 required to remove FWD and REV from reads. Remove up to 2 adapters from each read. (Default: 1)
        # With -n 2, it is possible that two 5’ or two 3’ adapters are removed from a read. Linked adapters prevent this.
        # If one has a primer mix (with multiple primer pairs), there is no way to specify which 5’ adapter goes with which 3’ adapter. With linked adapters, one can just use multiple -a options.
        cutadapt -g "$fwd_primer;min_overlap=$min_overlap_forward...$(reverse_complement "$rev_primer");min_overlap=$min_overlap_reverse" \
             --error-rate 0.15 --trimmed-only --cores "$CORES" \
             -o "$MARKER_DIR/${sample}_${marker}_fwd.fastq" \
             "$SAMPLE_DIR/${sample}_demultiplexed.fastq" > "$MARKER_DIR/${sample}_${marker}_primer_fwd.log" 2>&1 || true

        # Trim primers: Reverse orientation
        cutadapt -g "$rev_primer;min_overlap=$min_overlap_reverse...$(reverse_complement "$fwd_primer");min_overlap=$min_overlap_forward" \
             --error-rate 0.15 --trimmed-only --cores "$CORES" \
             -o "$MARKER_DIR/${sample}_${marker}_rev.fastq" \
             "$SAMPLE_DIR/${sample}_demultiplexed.fastq" > "$MARKER_DIR/${sample}_${marker}_primer_rev.log" 2>&1 || true

        # Reverse complement the reverse orientation
        seqkit seq -r -p "$MARKER_DIR/${sample}_${marker}_rev.fastq" -o "$MARKER_DIR/${sample}_${marker}_rev_rc.fastq" -t DNA 2>/dev/null || true

        # Combine orientations
        cat "$MARKER_DIR/${sample}_${marker}_fwd.fastq" "$MARKER_DIR/${sample}_${marker}_rev_rc.fastq" > "$MARKER_DIR/${sample}_${marker}_combined.fastq" 2>/dev/null || true

        # Quality trim
        fastq_quality_trimmer -t "$QUALITY_THRESHOLD" -l "$MIN_SIZE" -i "$MARKER_DIR/${sample}_${marker}_combined.fastq" -o "$MARKER_DIR/${sample}_${marker}_filtered.fastq" 2>/dev/null || true

        # Length filter
        seqkit seq -m "$MIN_SIZE" -M "$MAX_SIZE" -t DNA \
                  -i "$MARKER_DIR/${sample}_${marker}_filtered.fastq" \
                  -o "$MARKER_DIR/${sample}_${marker}_final.fastq" 2>/dev/null || true

        # Run NGSpeciesID
        # Initialize conda if needed
        source $(conda info --base 2>/dev/null)/etc/profile.d/conda.sh 2>/dev/null || true
        conda activate NGSpeciesID
        . ~/medaka/bin/activate
        NGSpeciesID --t "$CORES" --ont --consensus --sample_size "$SAMPLE_SIZE" --medaka --debug \
                    --fastq "$MARKER_DIR/${sample}_${marker}_final.fastq" \
                    --outfolder "$MARKER_DIR/NGSpeciesID_output" 2>&1 | tee "$MARKER_DIR/${sample}_${marker}_NGSpeciesID.log"
        conda deactivate

        # Concatenate consensus
        find "$MARKER_DIR/NGSpeciesID_output" -name "consensus_reference_*.fasta" -exec cat {} + > "$MARKER_DIR/${sample}_${marker}_consensus.fasta" 2>/dev/null || true

        # Filter consensuses by minimum supporting reads
        awk -v min_support="$MIN_CONSENSUS_SUPPORT" '
        /^>/ {
            if (match($0, /total_supporting_reads_([0-9]+)/, arr)) {
                support = arr[1] + 0  # convert to number
                if (support >= min_support) {
                    print_flag = 1
                    print
                } else {
                    print_flag = 0
                }
            } else {
                print_flag = 0  # if no match, skip
            }
        }
        !/^>/ {
            if (print_flag) print
        }
        ' "$MARKER_DIR/${sample}_${marker}_consensus.fasta" > "$MARKER_DIR/${sample}_${marker}_consensus_filtered.fasta"

    done
done

# Create consensuses table
CONSENSUS_TABLE="$OUT_DIR/all_consensuses.tsv"
echo -e "Sample\tMarker\tSupporting_Reads\tConsensus_ID\tSequence" > "$CONSENSUS_TABLE"
find "$OUT_DIR" -name "*_consensus_filtered.fasta" | while read fa; do
    sample=$(basename "$(dirname "$(dirname "$fa")")")
    marker=$(basename "$(dirname "$fa")")
    awk -v s="$sample" -v m="$marker" '/^>/{id=$0; if(match(id, /total_supporting_reads_([0-9]+)/, arr)) support=arr[1]; else support=0; getline seq; print s "\t" m "\t" support "\t" id "\t" seq}' "$fa" >> "$CONSENSUS_TABLE"
done

echo "Pipeline completed. Outputs in $OUT_DIR. Consensuses table: $CONSENSUS_TABLE"

# ---------------------------
# BLAST taxonomy assignment for consensuses
# ---------------------------
echo "Assigning taxonomy to consensuses via BLAST..."

# Create FASTA from consensuses table with a unique, deterministic query ID per sequence
# Use QueryID: Sample|Marker|Consensus_ID (Consensus_ID without leading ">") to avoid collisions across samples/markers
CONSENSUS_FASTA="$OUT_DIR/consensuses.fasta"
awk -F'\t' 'NR==1{next} {id=$4; gsub(/^>/, "", id); printf ">%s|%s|%s\n%s\n", $1, $2, id, $5}' "$CONSENSUS_TABLE" > "$CONSENSUS_FASTA"

# BLAST against nt (prefer megablast and limit to fungi taxon when taxonomy files available)
BLAST_OUT="$OUT_DIR/blast_results.tsv"
# Ensure BLAST can locate taxonomy files alongside DB
TAXDB_DIR=$(dirname "$NTDB")
if [[ -n "${BLASTDB:-}" ]]; then
    export BLASTDB="$TAXDB_DIR:$BLASTDB"
else
    export BLASTDB="$TAXDB_DIR"
fi
echo "BLASTDB set to: $BLASTDB"

TAXIDS_OPT=()
if [[ -f "$TAXDB_DIR/taxdb.bti" && -f "$TAXDB_DIR/taxdb.btd" ]] || [[ -f "$TAXDB_DIR/taxonomy4blast.sqlite3" ]]; then
    TAXIDS_OPT=( -taxids "$BLAST_TAXIDS" )
    echo "Running blastn (task=$BLAST_TASK; taxids=$BLAST_TAXIDS)"
else
    echo "Warning: taxonomy files not found in $TAXDB_DIR; running BLAST without -taxids filtering"
fi

# Try with taxids when available; on failure, retry without
if ! blastn -task "$BLAST_TASK" "${TAXIDS_OPT[@]}" -query "$CONSENSUS_FASTA" -db "$NTDB" -out "$BLAST_OUT" \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle' \
        -evalue 1e-20 -perc_identity 80 -qcov_hsp_perc 80 -max_target_seqs 10 -num_threads "$CORES"; then
    if [[ -n "${TAXIDS_OPT[*]}" ]]; then
        echo "Warning: blastn with -taxids failed; retrying without taxonomic filtering"
        blastn -task "$BLAST_TASK" -query "$CONSENSUS_FASTA" -db "$NTDB" -out "$BLAST_OUT" \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle' \
            -evalue 1e-20 -perc_identity 80 -qcov_hsp_perc 80 -max_target_seqs 10 -num_threads "$CORES" || echo "Warning: blastn failed"
    else
        echo "Warning: blastn failed"
    fi
fi

# Get top hit per query
BLAST_TOP="$OUT_DIR/blast_top.tsv"
sort -k1,1 -k12,12nr "$BLAST_OUT" | awk '!seen[$1]++' > "$BLAST_TOP"

# Extract taxids
cut -f13 "$BLAST_TOP" | tr ';' '\n' | grep -v '^$' | sort -u > "$OUT_DIR/taxids.txt"

# Get lineages with taxonkit
taxonkit lineage --data-dir "$TAXONKIT_DB" "$OUT_DIR/taxids.txt" > "$OUT_DIR/taxids_lineage.tsv"
taxonkit reformat -i 2 -f '{k};{p};{c};{o};{f};{g};{s}' --data-dir "$TAXONKIT_DB" "$OUT_DIR/taxids_lineage.tsv" > "$OUT_DIR/taxids_reformat.tsv"

# Create taxid to taxonomy map
awk -F'\t' 'NF>=3 {split($3, ranks, ";"); print $1 "\t" ranks[1] "\t" ranks[2] "\t" ranks[3] "\t" ranks[4] "\t" ranks[5] "\t" ranks[6] "\t" ranks[7]}' "$OUT_DIR/taxids_reformat.tsv" > "$OUT_DIR/taxid_to_tax.tsv"

# Map BLAST to taxonomy (keyed by composite QueryID)
BLAST_TAX="$OUT_DIR/blast_taxonomy.tsv"
awk -F'\t' 'BEGIN {
    OFS="\t";
    while((getline < "'$OUT_DIR'/taxid_to_tax.tsv") > 0) {
        tax[$1] = $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8;
    }
}
{
    split($13, ids, ";");
    tid = ids[1];
    if(tid in tax) print $1, tax[tid];
    else print $1, "", "", "", "", "", "", "";
}' "$BLAST_TOP" > "$BLAST_TAX"

# Merge taxonomy into consensuses table using an exact key-based join (no grep)
FINAL_TABLE="$OUT_DIR/all_consensuses_with_tax.tsv"
{
  echo -e "Sample\tMarker\tSupporting_Reads\tConsensus_ID\tSequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies"
  awk -F'\t' -v BLTAX="$BLAST_TAX" 'BEGIN{
      OFS="\t";
      while((getline < BLTAX) > 0){
          key=$1; $1=""; sub(/^\t/, ""); m[key]=$0; # store taxonomy by QueryID
      }
  }
  NR==1{next}
  NR>1{
      id=$4; gsub(/^>/, "", id);
      key = $1 "|" $2 "|" id;
      tax = (key in m) ? m[key] : "\t\t\t\t\t\t";
      print $1, $2, $3, $4, $5, tax;
  }' "$CONSENSUS_TABLE"
} > "$FINAL_TABLE"

echo "Final table with taxonomy: $FINAL_TABLE"

# ---------------------------
# Create final table (barcode_summary.tsv)
# ---------------------------
echo "Creating table (barcode_summary.tsv)"
mkdir -p "$OUT_DIR/final_results"
BARCODE_SUMMARY="$OUT_DIR/final_results/barcode_summary.tsv"

python3 - "$CONSENSUS_TABLE" "$BLAST_TOP" "$OUT_DIR/taxid_to_tax.tsv" "$BARCODE_SUMMARY" <<'PY'
import csv
import sys
from pathlib import Path
from collections import defaultdict

cons_path = Path(sys.argv[1])
blast_top_path = Path(sys.argv[2])
tax_map_path = Path(sys.argv[3])
out_path = Path(sys.argv[4])

HEADERS = [
    "Sample", "Marker", "Rank", "Size", "SeqID", "QueryID", "Sequence",
    "FoundMarkerBLAST", "SubjectID", "SubjectTitle", "Pident", "AlignLength",
    "Bitscore", "Evalue", "TaxID", "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
]

def read_consensus(path: Path):
    rows = []
    with path.open() as h:
        r = csv.DictReader(h, delimiter='\t')
        for row in r:
            rows.append(row)
    return rows

def read_blast_top(path: Path):
    rows = {}
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with path.open() as h:
        for line in h:
            p = line.rstrip('\n').split('\t')
            if len(p) < 14:
                continue
            q = p[0]
            rows[q] = {
                'SubjectID': p[1],
                'Pident': p[2],
                'AlignLength': p[3],
                'Evalue': p[10],
                'Bitscore': p[11],
                'TaxID': p[12].split(';')[0] if p[12] else '',
                'SubjectTitle': p[13],
            }
    return rows

def read_tax_map(path: Path):
    m = {}
    if not path.exists() or path.stat().st_size == 0:
        return m
    with path.open() as h:
        for row in csv.reader(h, delimiter='\t'):
            if not row:
                continue
            if len(row) >= 8:
                m[row[0]] = row[1:8]
    return m

cons_rows = read_consensus(cons_path)
blast_rows = read_blast_top(blast_top_path)
tax_map = read_tax_map(tax_map_path)

# Rank by Supporting_Reads per (Sample, Marker)
groups = defaultdict(list)
for r in cons_rows:
    try:
        supp = int(r.get('Supporting_Reads','0'))
    except Exception:
        supp = 0
    groups[(r.get('Sample',''), r.get('Marker',''))].append((supp, r))

rank_map = {}
for key, entries in groups.items():
    entries.sort(key=lambda x: (-x[0]))
    for idx, (_, r) in enumerate(entries, start=1):
        rank_map[(r.get('Sample',''), r.get('Marker',''), r.get('Consensus_ID',''))] = idx

out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open('w', newline='') as h:
    w = csv.writer(h, delimiter='\t')
    w.writerow(HEADERS)
    for r in cons_rows:
        sample = r.get('Sample','')
        marker = r.get('Marker','')
        seqid = r.get('Consensus_ID','').lstrip('>')
        # Composite query key to match BLAST qseqid: Sample|Marker|SeqID
        query = f"{sample}|{marker}|{seqid}"
        seq = r.get('Sequence','')
        size = r.get('Supporting_Reads','')
        rank = rank_map.get((sample, marker, r.get('Consensus_ID','')), '')

        b = blast_rows.get(query, {})
        title = b.get('SubjectTitle','')
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

        ranks = tax_map.get(b.get('TaxID',''), ['']*7)

        w.writerow([
            sample,
            marker,
            rank,
            size,
            seqid,
            query,
            seq,
            found,
            b.get('SubjectID',''),
            title,
            b.get('Pident',''),
            b.get('AlignLength',''),
            b.get('Bitscore',''),
            b.get('Evalue',''),
            b.get('TaxID',''),
            ranks[0] if len(ranks)>0 else '',
            ranks[1] if len(ranks)>1 else '',
            ranks[2] if len(ranks)>2 else '',
            ranks[3] if len(ranks)>3 else '',
            ranks[4] if len(ranks)>4 else '',
            ranks[5] if len(ranks)>5 else '',
            ranks[6] if len(ranks)>6 else '',
        ])

print(f"Wrote multiplex-style final table: {out_path}")
PY

echo "Multiplex-style final table generated: $BARCODE_SUMMARY"