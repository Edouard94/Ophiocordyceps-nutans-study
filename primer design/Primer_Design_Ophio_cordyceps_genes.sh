chmod +x gene.sh 
nohup bash gene.sh > gene.out 2> gene.log &

#!/usr/bin/env bash
set -euo pipefail

# Download fungi and hemiptera gene sequences from NCBI

## Create output directories
mkdir -p logs gb uniq_fasta cdhit_out alignments entropy_plots fasta

# Download taxonomy file and set path
obitaxonomy --download-ncbi --out ncbitaxo_20250724.tgz
TAXO="$(pwd)/ncbitaxo_20250724.tgz"

# Genes and organism queries
declare -A FUNGI_QUERIES=(
    [TEF]="Ophiocordyceps[Organism] OR Cordyceps[Organism] (tef1[All Fields] OR tef[All Fields] OR translation elongation factor 1 alpha[All Fields]) AND 0:8000[SLEN]"
    [RPB1]="Ophiocordyceps[Organism] OR Cordyceps[Organism] AND (rpb1[All Fields] OR rpoB[All Fields] OR RNA polymerase II subunit RPB1[All Fields]) AND 0:8000[SLEN]"
    [RPB2]="Ophiocordyceps[Organism] OR Cordyceps[Organism] AND (rpb2[All Fields] OR rpoC[All Fields] OR RNA polymerase II subunit RPB2[All Fields]) AND 0:8000[SLEN]"
    [LSU]="Ophiocordyceps[Organism] OR Cordyceps[Organism] AND (LSU[All Fields] OR large subunit ribosomal RNA[All Fields]) AND 0:8000[SLEN]"
    [SSU]="Ophiocordyceps[Organism] OR Cordyceps[Organism] AND (SSU[All Fields] OR small subunit ribosomal RNA[All Fields]) AND 0:8000[SLEN]"
    [ITS2]="Ophiocordyceps[Organism] OR Cordyceps[Organism] AND (ITS2[All Fields] OR internal transcribed spacer 2[All Fields]) AND 0:8000[SLEN]"
)
declare -A INSECT_QUERIES=(
    [TEF]="Hemiptera[Organism] AND (tef1[All Fields] OR tef[All Fields] OR translation elongation factor 1 alpha[All Fields]) AND 0:8000[SLEN]"
    [RPB1]="Hemiptera[Organism] AND (rpb1[All Fields] OR rpoB[All Fields] OR RNA polymerase II subunit RPB1[All Fields]) AND 0:8000[SLEN]"
    [RPB2]="Hemiptera[Organism] AND (rpb2[All Fields] OR rpoC[All Fields] OR RNA polymerase II subunit RPB2[All Fields]) AND 0:8000[SLEN]"
    [LSU]="Hemiptera[Organism] AND (LSU[All Fields] OR large subunit ribosomal RNA[All Fields]) AND 0:8000[SLEN]"
    [SSU]="Hemiptera[Organism] AND (SSU[All Fields] OR small subunit ribosomal RNA[All Fields]) AND 0:8000[SLEN]"
    [ITS2]="Hemiptera[Organism] AND (ITS2[All Fields] OR internal transcribed spacer 2[All Fields]) AND 0:8000[SLEN]"
)

GENES=(TEF RPB1 RPB2 LSU SSU ITS2)

# Function to robustly filter WGS records from GenBank files using awk for obiconvert
filter_wgs_records() {
    infile="$1"
    outfile="$2"
    awk '
        BEGIN { rec=""; skip=0 }
        /^LOCUS/ { if (rec != "" && !skip) print rec; rec=$0 "\n"; skip=0; next }
        /^WGS/ || /KEYWORDS[ ]+WGS[.]/ || /^TLS/ || /KEYWORDS[ ]+.*TLS/ { skip=1 }
        { if (NR > 1) rec = rec $0 "\n" }
        END { if (rec != "" && !skip) print rec }
    ' "$infile" > "$outfile"
}

mkdir -p stats

for gene in "${GENES[@]}"; do
    # Fungi
    esearch -db nucleotide -query "${FUNGI_QUERIES[$gene]}" | efetch -format gb > "gb/ophio_cordyceps_${gene}.gb" 2> "logs/ophio_${gene,,}.log"
    filter_wgs_records "gb/ophio_cordyceps_${gene}.gb" "gb/ophio_cordyceps_${gene}.filtered.gb"
    obiconvert --skip-empty --update-taxid -t "$TAXO" "gb/ophio_cordyceps_${gene}.filtered.gb" > "fasta/ophio_cordyceps_${gene}.obi.fasta"
    seqkit rmdup -s -i -o "uniq_fasta/uniq${gene}.fasta" "fasta/ophio_cordyceps_${gene}.obi.fasta"
    cd-hit-est -i "uniq_fasta/uniq${gene}.fasta" -o "cdhit_out/${gene}_c99.fasta" -c 0.99 -M 0 -T 0
    seqkit stats "cdhit_out/${gene}_c99.fasta" > "stats/${gene}_fungi.stats"
    mafft --thread -1 --auto "cdhit_out/${gene}_c99.fasta" > "alignments/${gene}_c99.aln" 2> "logs/mafft_${gene,,}.log"
    trimal -in "alignments/${gene}_c99.aln" -out "alignments/${gene}_c99_trimmed.aln" -automated1 -fasta
    plotcon -sequence "alignments/${gene}_c99_trimmed.aln" -winsize 10 -graph png -goutfile "entropy_plots/${gene}_entropy"

    # Hemiptera
    esearch -db nucleotide -query "${INSECT_QUERIES[$gene]}" | efetch -format gb > "gb/hemiptera_${gene}.gb" 2> "logs/hemiptera_${gene,,}.log"
    filter_wgs_records "gb/hemiptera_${gene}.gb" "gb/hemiptera_${gene}.filtered.gb"
    obiconvert --skip-empty --update-taxid -t "$TAXO" "gb/hemiptera_${gene}.filtered.gb" > "fasta/hemiptera_${gene}.obi.fasta"
    seqkit rmdup -s -i -o "uniq_fasta/uniqHemiptera${gene}.fasta" "fasta/hemiptera_${gene}.obi.fasta"
    cd-hit-est -i "uniq_fasta/uniqHemiptera${gene}.fasta" -o "cdhit_out/Hemiptera_${gene}_c99.fasta" -c 0.99 -M 0 -T 0
    seqkit stats "cdhit_out/Hemiptera_${gene}_c99.fasta" > "stats/${gene}_hemiptera.stats"
done

# Combine fungi and hemiptera clustered sequences for each gene
for gene in "${GENES[@]}"; do
    cat "cdhit_out/${gene}_c99.fasta" "cdhit_out/Hemiptera_${gene}_c99.fasta" > "cdhit_out/${gene}_combined.fasta"
done

#!/usr/bin/env bash
set -euo pipefail

mkdir -p entropy_plots1

GENES=(TEF RPB1 RPB2 LSU SSU ITS2)
for gene in "${GENES[@]}"; do
plotcon -sequence "alignments/${gene}_c99.aln" -winsize 10 -graph png -goutfile "entropy_plots1/${gene}_entropy.png"
done

# ──────────────────────────────────────────────────────────────
# MBC-prime analysis of combined Hemiptera and Fungi sequences
# ──────────────────────────────────────────────────────────────
#!/usr/bin/env bash
set -euo pipefail
mkdir -p mbc_alignments mbc_out
eval "$(conda shell.bash hook)"
conda activate mbc-prime

GENES=(TEF RPB1 RPB2 LSU SSU ITS2)

for gene in "${GENES[@]}"; do
    # Combine and align sequences quickly
    cat "cdhit_out/${gene}_c99.fasta" "cdhit_out/Hemiptera_${gene}_c99.fasta" > "mbc_alignments/${gene}_combined.fasta"
    mafft --retree 1 --maxiterate 0 --thread -1 "mbc_alignments/${gene}_combined.fasta" > "mbc_alignments/${gene}_combined.aln"

    # Count number of fungi and hemiptera sequences
    fungi_count=$(grep -c "^>" "cdhit_out/${gene}_c99.fasta")
    hemiptera_count=$(grep -c "^>" "cdhit_out/Hemiptera_${gene}_c99.fasta")

    # Run mbc-prime (adjust -t and other options as needed)
    nohup ~/mbc-prime/mbc-prime -t "$fungi_count" -s 0.5 -v \
        "mbc_alignments/${gene}_combined.aln" >> "mbc_out/primers_${gene}.tsv" 2> "mbc_out/mbc-prime_${gene}.log" &
done

# ──────────────────────────────────────────────────────────────
# ecoPrimers primer design for ~300 bp amplicons (all genes)
# ──────────────────────────────────────────────────────────────
#!/usr/bin/env bash
set -euo pipefail

TAXO="$(pwd)/ncbitaxo_20250724.tgz"
GENES=(TEF RPB1 RPB2 LSU SSU ITS2)

# Requires: ecoPrimers and OBITools4 (install with: mamba install -c bioconda ecoprimers obitools4)
# Output: For each gene, will produce a list of primer pairs for ~300 bp products

# Re-annotation of sequences to species level
mkdir -p obiannotate_out
for gene in "${GENES[@]}"; do
    obiannotate -t "$TAXO" \
                --with-taxon-at-rank species \
                "cdhit_out/${gene}_combined.fasta" | \
    obiannotate -S 'ori_taxid=annotations.taxid' | \
    obiannotate -S 'taxid=annotations.species_taxid' | \
    obiuniq -c taxid > "obiannotate_out/${gene}_obi_ann.fasta"
done

# Formatting data for ecoPrimers

## Unarchiving the taxonomy
## The old OBITools cannot use archived and compressed taxonomies
mkdir ncbitaxo_20250724/
cd ncbitaxo_20250724/
tar zxvf ../ncbitaxo_20250724.tgz
cd ../

# When the -O option is added to a OBITools4 command, the old OBITools format is used instead of the new JSON-based format.
for gene in "${GENES[@]}"; do
    obiconvert -O "obiannotate_out/${gene}_obi_ann.fasta" > "obiannotate_out/${gene}_obi_ann.old.fasta"
done

# Creating gene databases for ecoPrimers
mkdir -p ecoPrimers_db
for gene in "${GENES[@]}"; do
    echo "Creating ecoPrimers database for $gene"
    /home/edouard/ecoPCRFormat -t ncbitaxo_20250724 \
                 -f \
                 -n "ecoPrimers_db/Hypocreales_Hemiptera_${gene}" \
                 "obiannotate_out/${gene}_obi_ann.old.fasta"
done


# Designing primers for Hypocreales taxid 5125 while avoiding Hemipetera taxid 7524
mkdir -p ecoPrimers_out
for gene in "${GENES[@]}"; do
    echo "Designing less strict primers for $gene"
    ecoPrimers -d "ecoPrimers_db/Hypocreales_Hemiptera_${gene}" \
    -e 5 -O 18 \
    -q 0.5 -s 0.5 \
    -l 150 -L 500 \
    -r 5125 -i 7524 \
    > "ecoPrimers_out/Primers_${gene}_moderate.ecoprimers"
done

# ──────────────────────────────────────────────────────────────
# PrimerProspector analysis of the designed primers
# ──────────────────────────────────────────────────────────────

#!/usr/bin/env bash
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate primerprospector

GENES=(TEF RPB1 RPB2 LSU SSU ITS2)

#GENES=(RPB2)

mkdir -p PrimerProspector
for gene in "${GENES[@]}"; do
    seqkit seq -u "cdhit_out/${gene}_c99.fasta" > "PrimerProspector/${gene}_c99.upper.fasta"
done

cd PrimerProspector

# Create primers txt files for each gene
cat > Primers_ITS2.txt << EOF
ITS3tagmixf	gcatcgatgaagaacgcagc
ITS4ngsr	cttttcctccgcttattgatatgc
ITS2f_350bp	gcatcgatgaagaacgcagc
ITS2r_350bp	yttttcctccgcttattgatatgc
ITS2f_230bp	arcaacggatctcttggy
ITS2r_230bp	ccgccactgcatttcggg
EOF
cat > Primers_SSU.txt << EOF
SSUf	ccgtggtaattctagagctaatacatgc
SSUr	rtcgggattgggtaatttgcgc
EOF
cat > Primers_LSU.txt << EOF
LSUf	agcgcacaagtagagtgatcg
LSUr	ggtccgtgtttcaagacgg
LsuF_ophio_poland	AGAGTGATCGAAAKRYG
LsuR_ophio_poland	GGCATAGTTCACCTTCTT
EOF
cat > Primers_TEF.txt << EOF
TEFf	tggtacaagggctgggagaagg
TEFr	gctgctcgtggtgcatttcg
mbcf_tef_s81_poland	GAGGCTGGTATYTCCAARGA
mbcr_tef_s79_poland	AGCTGCTCGTGGTGCATYTC
EOF
cat > Primers_RPB1.txt << EOF
RPB1f	ccgtgtgcaagaagaagcg
RPB1r	cggtgatgatcatccactcgg
mbcf_rpb1_s84_poland	GCCAYAACTGCAGCAARGT
mbcr_rpb1_s40_poland	CRGTGAKRATCATCCAYTCRGG
Rpb1f_opti_poland	CAYAACTGCAGCAARGT
Rpb1r_opti_poland	RTTGAGVCCCATGTTG
EOF
cat > Primers_RPB2.txt << EOF
RPB2f	gtctctcatgtgctacgtcagcg
RPB2r	ccgagaaaatcttgaactcctgatcc
EOF

# Analyze primers for each gene
for gene in "${GENES[@]}"; do
    analyze_primers.py -P Primers_${gene}.txt -f ${gene}_c99.upper.fasta -v -o ${gene}_analyze_primers
done

# Optimize primer pairs for each gene
for gene in "${GENES[@]}"; do
    for primer in $(awk '{print $1}' Primers_${gene}.txt); do
        optimize_primers.py -i ${gene}_analyze_primers/${primer}_${gene}_c99_hits.txt -t 5 -v -s weighted_score
    done
done

# Summarise primer hits for all genes
GENES=(TEF RPB1 RPB2 LSU SSU ITS2)

for gene in "${GENES[@]}"; do
  echo -e "Primer\tHits\tAvgWeightedScore\tAvgNon3Mm\tAvg3Mm" > "${gene}_primer_summary.tsv"
  for hits in ${gene}_analyze_primers/*_${gene}_c99_hits.txt; do
    p=$(basename "$hits" _${gene}_c99_hits.txt)
    n=$(grep -vc '^#' "$hits")
    read w nm t3 < <(
      awk -F, -v N="$n" '
        $1 !~ /^#/ { sum10 += $10; sum5 += $5; sum6 += $6 }
        END {
          if (N>0)
            printf("%.2f\t%.2f\t%.2f\n", sum10/N, sum5/N, sum6/N)
          else
            print "NA\tNA\tNA"
        }
      ' "$hits"
    )
    printf "%s\t%d\t%s\t%s\t%s\n" "$p" "$n" "$w" "$nm" "$t3"
  done >> "${gene}_primer_summary.tsv"
done

# Combine base frequencies for each gene
for gene in "${GENES[@]}"; do
    cat *_${gene}_c99_hits_base_frequencies.txt > ${gene}_base_frequencies.txt
done
# ──────────────────────────────────────────────────────────────
