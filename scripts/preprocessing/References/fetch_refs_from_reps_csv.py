#!/usr/bin/env python3
"""
Fetch reference sequences from NCBI for multigene phylogeny analysis.

This script reads a CSV of representative species with accessions and vouchers,
then fetches reference sequences from NCBI for specified genes (ITS, SSU, LSU, TEF, RPB1, RPB2).
For each species, it uses direct accession numbers when available, or performs
voucher-based searches for missing accessions. It can also optionally include
outgroup taxa, audit vouchers missing ITS accessions, and populate additional
representatives for specific species (including synonyms) by finding vouchers with multiple genes.

Expects CSV format: Species,Voucher,SSU,LSU,TEF,RPB1,RPB2,Host

Usage:

    nohup conda run -n base python3 scripts/fetch_refs_from_reps_csv.py \
        --csv Ophiocordyceps_species_representatives.csv \
        --outdir references_last \
        --genes ITS SSU LSU TEF RPB1 RPB2 \
        --include-outgroups \
        --audit-its \
        --populate-species "Ophiocordyceps nutans" "Cordyceps nutans" "Ophiocordyceps neonutans" \
        --min-genes 1 \
        > fetch_refs.log 2>&1 &

    # Monitor progress:
    tail -f fetch_refs.log

Options:
    --csv PATH                 Path to representatives CSV (required)
    --outdir DIR               Output directory (default: references)
    --genes GENES              Genes to fetch (default: ITS SSU LSU TEF RPB1)
    --strict-voucher           Strict voucher matching
    --include-outgroups        Include outgroup references
    --outgroups SPECIES...     Outgroup species (default: Tolypocladium inflatum Tolypocladium ophioglossoides)
    --audit-its                Audit missing ITS accessions
    --parse-only               Dry-run without NCBI calls
    --populate-species SPECIES...  Populate additional representatives
    --min-genes N              Min genes for additional reps (default: 3)

Output: Creates {GENE}_refs.fasta files with headers formatted as:
    >SpeciesName_VoucherID_Gene_Accession

Notes:
    - Requires NCBI E-utilities: esearch, efetch (unless using --parse-only)
    - Prefer running via conda base env where entrez-direct is installed:
        conda run -n base ...
    - If missing, install EDirect:
        conda install -n base -c bioconda entrez-direct
    - Use --parse-only to validate input without contacting NCBI
"""
import csv
import sys
import subprocess as sp
from pathlib import Path
from collections import defaultdict
import argparse
import re
import xml.etree.ElementTree as ET

# Shared constants
# Fields in INSDQualifier_name that commonly carry voucher-like identifiers.
# Includes strain/isolate/culture_collection/bio_material and also note/clone to catch codes
# that are sometimes only recorded in generic notes.
VOUCHER_FIELDS = [
    "specimen_voucher",
    "strain",
    "isolate",
    "culture_collection",
    "bio_material",
    "note",
    "clone",
]


GENE_MARKERS = {
    "ITS": "(ITS[All Fields] OR ITS2[All Fields] OR internal transcribed spacer[All Fields] OR internal transcribed spacer 2[All Fields])",
    "SSU": "(SSU[All Fields] OR 18S[All Fields] OR small subunit[All Fields])",
    "LSU": "(LSU[All Fields] OR 28S[All Fields] OR large subunit[All Fields])",
    "TEF": "(tef1[All Fields] OR tef[All Fields] OR tef1a[All Fields] OR \"tef1-alpha\"[All Fields] OR translation elongation factor 1 alpha[All Fields])",
    "RPB1": "(RPB1[All Fields] OR \"RNA polymerase II subunit B1\"[All Fields] OR \"DNA-directed RNA polymerase II largest subunit\"[All Fields])",
    "RPB2": "(RPB2[All Fields] OR \"RNA polymerase II subunit B2\"[All Fields] OR \"DNA-directed RNA polymerase II second largest subunit\"[All Fields])",
}


def need(bin_name):
    try:
        sp.run([bin_name, "-h"], stdout=sp.DEVNULL, stderr=sp.DEVNULL, check=False)
    except FileNotFoundError:
        sys.stderr.write(f"Error: required tool not found: {bin_name}\n")
        sys.exit(1)


def norm(s):
    return re.sub(r"[^A-Za-z0-9]", "", s or "").lower()

def canonical_voucher(v: str) -> str:
    """Return a canonical voucher string:
    - Replace whitespace/semicolons with single underscore
    - Collapse multiple underscores
    - Merge pattern CODE_DIGITS where CODE is letters and DIGITS is numeric/alphanum (e.g. OSC_71233 -> OSC71233)
    - Preserve other internal underscores (e.g. complex composite vouchers) but remove trailing/leading underscores.
    """
    if not v:
        return "NoVoucher"
    # whitespace/semicolon to underscore
    v2 = re.sub(r"[\s;]+", "_", v.strip())
    v2 = re.sub(r"_+", "_", v2)
    # Merge final letter+underscore+digits token if present
    # e.g. ABC_12345 or ABC_12345A at END
    v2 = re.sub(r"([A-Za-z])_([0-9][0-9A-Za-z]*)$", r"\1\2", v2)
    # Also handle cases like OSC_71233 -> OSC71233 when preceding part ends with letters
    v2 = re.sub(r"([A-Za-z]+)_([0-9][0-9A-Za-z]*)$", r"\1\2", v2)
    v2 = v2.strip('_')
    return v2 or "NoVoucher"


def retry_on_failure(max_retries=3, delay=1):
    def decorator(func):
        def wrapper(*args, **kwargs):
            import time
            last_err = None
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    last_err = e
                    if attempt < max_retries - 1:
                        sys.stderr.write(f"[WARN] {func.__name__} failed (attempt {attempt+1}), retrying...\n")
                        time.sleep(delay)
            # If we reach here, all retries failed
            raise last_err
        return wrapper
    return decorator

@retry_on_failure()
def esearch_xml(query, save_path=None):
    # Use 'nuccore' instead of 'nucleotide' for compatibility
    # Use 'gbc' format to get INSDSeq XML format (not Bioseq format)
    p1 = sp.Popen(["esearch", "-db", "nuccore", "-query", query], stdout=sp.PIPE, stderr=sp.PIPE)
    p2 = sp.Popen(["efetch", "-db", "nuccore", "-format", "gbc"], stdin=p1.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    p1.stdout.close()
    out, err = p2.communicate()
    if p2.returncode != 0:
        raise RuntimeError(f"NCBI command failed: {err.decode(errors='replace')}")
    xml_text = out.decode("utf-8", errors="replace")

    # Save XML if path provided
    if save_path:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with open(save_path, "w", encoding="utf-8") as fh:
            fh.write(xml_text)

    return xml_text


def esearch_first_fasta(query):
    """Return (accession, sequence) for the first FASTA hit from esearch|efetch.
    Returns (None, None) if no FASTA entry is found.
    """
    p1 = sp.Popen(["esearch", "-db", "nuccore", "-query", query], stdout=sp.PIPE)
    p2 = sp.Popen(["efetch", "-db", "nuccore", "-format", "fasta"], stdin=p1.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    p1.stdout.close()
    out, _ = p2.communicate()
    fasta_text = out.decode("utf-8", errors="replace")
    if ">" not in fasta_text:
        return (None, None)
    acc = None
    seq_lines = []
    first = True
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if first:
                header = line[1:].strip()
                acc = header.split()[0] if header else None
                first = False
            else:
                break
        else:
            if not first:
                seq_lines.append(line.strip())
    seq = "".join(seq_lines)
    if acc and seq:
        return (acc, seq)
    return (None, None)


def efetch_fasta_batch(accessions):
    if not accessions:
        return ""
    ids = ",".join(accessions)
    res = sp.run(["efetch", "-db", "nuccore", "-id", ids, "-format", "fasta"], capture_output=True, text=True)
    return res.stdout


def valid_xml(xml_text: str) -> bool:
    return bool(xml_text and xml_text.strip() and xml_text.strip().startswith("<"))

def parse_xml_for_voucher(xml_text, voucher, log_prefix=""):
    """Return list of tuples (acc, organism, specimen_voucher_value, seq) filtered by voucher substring match.
    Checks voucher-related qualifiers + definition/title for a flexible substring match."""
    hits = []
    if not valid_xml(xml_text):
        if log_prefix:
            sys.stderr.write(f"[WARN] Empty or invalid XML response {log_prefix}\n")
        return hits
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        if log_prefix:
            sys.stderr.write(f"[WARN] XML parse error {log_prefix}: {str(e)[:100]}\n")
        return hits
    voucher_norm = norm(voucher)
    for insd in root.findall(".//INSDSeq"):
        acc = insd.findtext("INSDSeq_accession-version") or insd.findtext("INSDSeq_primary-accession") or ""
        org = insd.findtext("INSDSeq_organism") or ""
        seq = insd.findtext("INSDSeq_sequence") or ""
        definition = insd.findtext("INSDSeq_definition") or ""
        # Preferred: match voucher against common voucher-carrying qualifiers
        vvals = [
            q.findtext("INSDQualifier_value") or ""
            for q in insd.findall(".//INSDQualifier")
            if (q.findtext("INSDQualifier_name") or "").lower() in VOUCHER_FIELDS
        ]
        match = (any(voucher_norm in norm(vv) for vv in vvals) or voucher_norm in norm(definition)) if voucher_norm else True
        # Fallback: if still not matched, scan all qualifier values (some submitters only put codes in unusual fields)
        if (not match) and voucher_norm:
            all_qvals = [q.findtext("INSDQualifier_value") or "" for q in insd.findall(".//INSDQualifier")]
            match = any(voucher_norm in norm(vv) for vv in all_qvals)
        if match:
            hits.append((acc, org, vvals[0] if vvals else "", seq))
    return hits


def extract_qualifiers(insd):
    quals = []
    for feat in insd.findall(".//INSDFeature"):
        # Qualifiers are under INSDFeature_quals, not directly under INSDFeature
        quals_elem = feat.find("INSDFeature_quals")
        if quals_elem is not None:
            for q in quals_elem.findall("INSDQualifier"):
                name = q.findtext("INSDQualifier_name") or ""
                value = q.findtext("INSDQualifier_value") or ""
                quals.append({"name": name, "value": value})
    return quals


def detect_gene_from_text(text):
    """Detect gene marker from free text (gene name, product, note, title).
    Heuristic: if a record mentions mixed rDNA regions (18S/ITS/28S), prefer ITS (often the complete region) over SSU/LSU.
    Returns standardized gene name (ITS, SSU, LSU, TEF, RPB1, RPB2) or empty string."""
    if not text:
        return ""

    text_norm = text.lower().replace("-", " ").replace("_", " ")
    # collapse multiple spaces
    text_norm = re.sub(r"\s+", " ", text_norm)

    # rDNA markers
    has_its2 = ("its2" in text_norm) or ("internal transcribed spacer 2" in text_norm)
    has_its1 = ("its1" in text_norm) or ("internal transcribed spacer 1" in text_norm)
    has_its_any = (" its " in f" {text_norm} ") or ("internal transcribed spacer" in text_norm) or has_its1 or has_its2
    has_ssu = any(p in text_norm for p in [
        "18s", "ssu", "small subunit ribosomal rna", "small subunit rrna", "small subunit"
    ])
    has_lsu = any(p in text_norm for p in [
        "28s", "lsu", "large subunit ribosomal rna", "large subunit rrna", "large subunit"
    ])

    # Prefer ITS if present, especially in mixed rDNA titles
    if has_its2:
        return "ITS"
    if has_its_any:
        return "ITS"
    if has_ssu:
        return "SSU"
    if has_lsu:
        return "LSU"

    # TEF patterns (translation elongation factor 1-alpha)
    if any(pattern in text_norm for pattern in [
        " tef", " tef1", " ef1", " ef1a", " ef1alpha", " tef1alpha",
        "elongation factor 1", "elongation factor 1alpha",
        "translation elongation factor", "translation elongation factor 1alpha"
    ]):
        return "TEF"

    # RPB1 patterns (RNA polymerase II largest subunit)
    if any(pattern in text_norm for pattern in [
        "rpb1", "rna polymerase ii largest subunit",
        "rna polymerase ii subunit b1", "rna polymerase 2 largest subunit",
        "dna directed rna polymerase ii largest subunit"
    ]):
        return "RPB1"

    # RPB2 patterns (RNA polymerase II second largest subunit)
    if any(pattern in text_norm for pattern in [
        "rpb2", "rna polymerase ii second largest subunit",
        "rna polymerase ii subunit b2", "rna polymerase 2 second largest subunit",
        "dna directed rna polymerase ii second largest subunit"
    ]):
        return "RPB2"

    return ""


def parse_xml_entries(xml_text):
    """Yield dicts with acc, org, seq, voucher, gene from INSD XML.
    Handles multiple INSDSet root elements by wrapping them.
    Uses comprehensive gene detection from gene, product, and note fields.
    Extracts voucher info from specimen_voucher, strain, isolate, culture_collection, etc."""
    if not xml_text or not xml_text.strip() or not xml_text.strip().startswith('<'):
        return []
    
    # Handle multiple INSDSet elements by wrapping them
    # NCBI sometimes returns multiple result sets
    if xml_text.count('<INSDSet>') > 1:
        # Wrap multiple INSDSet elements in a single root
        # Remove XML declaration and DOCTYPE from all but first
        parts = xml_text.split('<INSDSet>')
        wrapped = '<INSDSet_wrapper>'
        for i, part in enumerate(parts):
            if i == 0:
                continue  # Skip text before first INSDSet
            # Remove closing tag and any following content before next INSDSet
            part_clean = part.split('</INSDSet>')[0] if '</INSDSet>' in part else part
            wrapped += '<INSDSet>' + part_clean + '</INSDSet>'
        wrapped += '</INSDSet_wrapper>'
        xml_text = wrapped
    
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []
    out = []
    for insd in root.findall(".//INSDSeq"):
        acc = insd.findtext("INSDSeq_accession-version") or insd.findtext("INSDSeq_primary-accession") or ""
        org = insd.findtext("INSDSeq_organism") or ""
        seq = insd.findtext("INSDSeq_sequence") or ""
        definition = insd.findtext("INSDSeq_definition") or ""
        voucher = ""
        gene = ""
        
        quals = extract_qualifiers(insd)
        
        # Extract voucher from multiple possible qualifier fields (in order of preference)
        for field in VOUCHER_FIELDS:
            if not voucher:
                voucher = next((q["value"] for q in quals if q["name"].lower() == field), "")
        
        # Try to detect gene from multiple sources (in order of preference)
        # 1. Try gene qualifier
        gene_qual = next((q["value"] for q in quals if q["name"].lower() == "gene"), "")
        if gene_qual:
            gene = detect_gene_from_text(gene_qual)
        
        # 2. Try product qualifier
        if not gene:
            product_qual = next((q["value"] for q in quals if q["name"].lower() == "product"), "")
            if product_qual:
                gene = detect_gene_from_text(product_qual)
        
        # 3. Try note qualifier
        if not gene:
            note_qual = next((q["value"] for q in quals if q["name"].lower() == "note"), "")
            if note_qual:
                gene = detect_gene_from_text(note_qual)
        
        # 4. Try sequence definition/title
        if definition:
            def_gene = detect_gene_from_text(definition)
            # If title indicates ITS (mixed rDNA record), prefer ITS over SSU/LSU detected from qualifiers
            if def_gene == "ITS":
                gene = "ITS"
            elif not gene:
                gene = def_gene
        
        out.append({"acc": acc, "org": org, "seq": seq, "voucher": voucher, "gene": gene})
    return out


def write_gene_fastas(outdir, gene_to_records):
    outdir.mkdir(parents=True, exist_ok=True)
    for gene, recs in gene_to_records.items():
        fasta_path = outdir / f"{gene}_refs.fasta"
        with open(fasta_path, "w", encoding="utf-8") as fh:
            for (species, acc, voucher, seq) in recs:
                safe_org = species.replace(" ", "_")
                safe_voucher = canonical_voucher(voucher)
                # Format: SpeciesName_VoucherID_Gene_Accession
                fh.write(f">{safe_org}_{safe_voucher}_{gene}_{acc}\n")
                # wrap sequence to 60 cols
                seq = re.sub(r"[^ACGTNatcgn]", "N", seq).upper()
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i+60] + "\n")


# Lightweight logging to both stderr and a file under outdir/logs
LOG_PATH = None

def ensure_log_path(outdir):
    logs_dir = outdir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    return logs_dir / "fetch.log"

def log(msg):
    line = msg if msg.endswith("\n") else msg + "\n"
    try:
        sys.stderr.write(line)
        sys.stderr.flush()
    except Exception:
        pass
    global LOG_PATH
    if LOG_PATH:
        try:
            with open(LOG_PATH, "a", encoding="utf-8") as lf:
                lf.write(line)
        except Exception:
            pass


def populate_additional_reps(species_list, genes, min_genes, rows):
    """Populate vouchers that have at least min_genes among requested genes across provided species names.
    Strategy: single broad query per species (all sequences length filter) then classify sequences by voucher.
    This avoids per-gene queries and reduces network calls dramatically.
    
    Key fix: Only assign ONE accession per voucher+gene combination (the first one found).
    This prevents mismatches where a voucher's gene gets assigned an accession from a different voucher.
    """
    xml_dir = Path("references") / "xml_cache"
    voucher_to_gene_acc = defaultdict(lambda: defaultdict(str))
    voucher_to_species = {}  # Track which species each voucher belongs to
    
    for species in species_list:
        safe_species = re.sub(r"[^\w\s-]", "", species).replace(" ", "_")
        xml_filename = f"populate_{safe_species}.xml"
        xml_path = xml_dir / xml_filename
        log(f"[INFO] Fetching sequences for {species}...")
        xml = esearch_xml(f'"{species}"[ORGN] AND 0:8000[SLEN]', save_path=xml_path)
        entries = parse_xml_entries(xml)
        log(f"[INFO] Found {len(entries)} total sequences for {species}")
        
        # Debug: check what we're finding
        genes_found = defaultdict(int)
        vouchers_found = set()
        voucher_gene_pairs = defaultdict(list)  # Track all acc per voucher+gene
        
        for e in entries:
            vval = e["voucher"]
            gene = e["gene"]
            acc = e["acc"]
            org = e["org"]
            
            if vval:
                vouchers_found.add(vval)
                # Track species for this voucher (prefer exact match)
                if vval not in voucher_to_species or org == species:
                    voucher_to_species[vval] = org
            
            if gene:
                genes_found[gene] += 1
            
            if not vval or not gene or not acc:
                continue
                
            gene = gene.upper()
            if gene in genes:
                voucher_gene_pairs[(vval, gene)].append(acc)
                # Only assign if not already set for this voucher+gene combo
                if not voucher_to_gene_acc[vval][gene]:
                    voucher_to_gene_acc[vval][gene] = acc
                    log(f"[DEBUG] Assigned {gene} acc {acc} to voucher {vval} (from {org})")
                else:
                    # Warn if we're skipping a different accession for same voucher+gene
                    if voucher_to_gene_acc[vval][gene] != acc:
                        log(f"[WARN] Skipping duplicate {gene} acc {acc} for voucher {vval} (already have {voucher_to_gene_acc[vval][gene]})")
        
        log(f"[INFO] Detected genes: {dict(genes_found)}")
        log(f"[INFO] Found {len(vouchers_found)} unique vouchers")
        log(f"[INFO] Vouchers with gene data: {len(voucher_to_gene_acc)}")
        
        # Check for vouchers with multiple accessions per gene
        for (vval, gene), accs in voucher_gene_pairs.items():
            if len(set(accs)) > 1:
                log(f"[WARN] Voucher {vval} has {len(set(accs))} different {gene} accessions: {set(accs)}")
        
    # Get existing vouchers (normalize for comparison)
    existing_vouchers = set()
    for row in rows:
        if row.get('Species') in species_list:
            existing_vouchers.add(row.get('Voucher', '').strip())
    
    log(f"[INFO] Existing vouchers to exclude: {existing_vouchers}")
    
    added = 0
    skipped = 0
    for voucher, gene_dict in voucher_to_gene_acc.items():
        gene_count = len([g for g in genes if gene_dict.get(g)])
        
        if gene_count < min_genes:
            skipped += 1
            log(f"[SKIP] Voucher {voucher} has only {gene_count} genes (need {min_genes}): {[g for g in genes if gene_dict.get(g)]}")
            continue
            
        if voucher in existing_vouchers:
            skipped += 1
            log(f"[SKIP] Voucher {voucher} already exists in CSV")
            continue
        
        # Use the species from NCBI (might be synonym) or first in list
        species_for_row = voucher_to_species.get(voucher, species_list[0])
        new_row = {'Species': species_for_row, 'Voucher': voucher}
        for g in genes:
            new_row[g] = gene_dict.get(g, '')
        new_row['Host'] = ''
        rows.append(new_row)
        
        genes_list = [f"{g}:{gene_dict[g]}" for g in genes if gene_dict.get(g)]
        log(f"[INFO] Added voucher {voucher} ({species_for_row}) with {gene_count} genes: {genes_list}")
        added += 1
        
    log(f"[INFO] Added {added} additional representatives for {species_list} (skipped {skipped})")


def main():
    ap = argparse.ArgumentParser(description="Fetch reference FASTAs per gene from a representatives CSV (species,voucher,genes). Accession preferred; fallback: voucher-based search (by default uses first sequence if no exact voucher match; use --strict-voucher to skip instead). Supports optional outgroup fetching, ITS voucher auditing, and populating additional representatives for specific species (including synonyms).")
    ap.add_argument("--csv", required=True, help="Path to Ophiocordyceps_species_representatives.csv")
    ap.add_argument("--outdir", default="references", help="Output directory for per-gene FASTAs")
    ap.add_argument("--genes", nargs="*", default=["ITS", "SSU", "LSU", "TEF", "RPB1"], help="Genes to fetch (default: ITS SSU LSU TEF RPB1)")
    ap.add_argument("--strict-voucher", action="store_true", help="Require exact voucher matches (no fallback) AND enforce a single primary voucher per species across all genes; sequences from other vouchers for the same species are skipped.")
    ap.add_argument("--include-outgroups", action="store_true", help="Include outgroup references (Tolypocladium spp.) for specified genes.")
    ap.add_argument("--outgroups", nargs="*", default=["Tolypocladium inflatum", "Tolypocladium ophioglossoides"], help="Outgroup species names to include (default: Tolypocladium inflatum, Tolypocladium ophioglossoides)")
    ap.add_argument("--audit-its", action="store_true", help="Audit vouchers missing ITS accession by searching vouchers with ITS/ITS2; write report.")
    ap.add_argument("--parse-only", action="store_true", help="Parse the representatives file and print summary without contacting NCBI.")
    ap.add_argument("--populate-species", nargs="*", help="Species names to populate with additional representatives from NCBI (can include synonyms)")
    ap.add_argument("--min-genes", type=int, default=3, help="Minimum number of genes required for additional representatives")
    args = ap.parse_args()

    if not args.parse_only:
        for bin_name in ("esearch", "efetch"):
            need(bin_name)

    reps_path = Path(args.csv)
    outdir = Path(args.outdir)
    audit_dir = outdir / "audits"
    # initialize log file
    global LOG_PATH
    LOG_PATH = ensure_log_path(outdir)
    log(f"[START] fetch_refs_from_reps_csv | csv={args.csv} | outdir={outdir}")

    genes = [g.upper() for g in args.genes]
    log(f"[INFO] Genes requested: {' '.join(genes)}")
    for g in genes:
        if g not in GENE_MARKERS:
            sys.stderr.write(f"Unsupported gene: {g}. Supported: {sorted(GENE_MARKERS)}\n")
            sys.exit(1)

    # Read representatives CSV
    rows = []
    with open(reps_path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        for row in r:
            rows.append(row)

    if args.populate_species:
        log(f"[INFO] Populating additional reps for: {', '.join(args.populate_species)} (min_genes={args.min_genes})")
        populate_additional_reps(args.populate_species, genes, args.min_genes, rows)

    # Collect accessions per gene (batch fetch) and voucher-search requests where accession missing
    gene_to_accessions = defaultdict(set)
    voucher_requests = []  # tuples: (species, voucher, gene)

    for row in rows:
        species = row.get("Species", "").strip()
        voucher = row.get("Voucher", "").strip()
        for g in genes:
            acc = (row.get(g, "") or "").strip()
            if acc:
                gene_to_accessions[g].add(acc)
            else:
                voucher_requests.append((species, voucher, g))

    if args.parse_only:
        total = len(rows)
        log(f"[PARSE-ONLY] Loaded {total} species/voucher rows")
        for g in genes:
            filled = sum(1 for row in rows if (row.get(g, "") or "").strip())
            log(f"[PARSE-ONLY] Rows with {g} accession: {filled}")
        return

    # Batch fetch all known accessions (concise implementation)
    gene_to_records = defaultdict(list)
    fetch_summary = defaultdict(list)
    for g, accs in gene_to_accessions.items():
        if not accs:
            continue
        log(f"[INFO] Fetching {len(accs)} direct accessions for {g}...")
        fasta = efetch_fasta_batch(sorted(accs))
        acc_to_seq = {}
        current = None
        buf = []
        for line in fasta.splitlines():
            if line.startswith(">"):
                if current is not None:
                    acc_to_seq[current] = "".join(buf)
                current = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if current is not None:
            acc_to_seq[current] = "".join(buf)

        # index accessions for mapping
        acc_index = {}
        for row in rows:
            if row.get(g, ""):
                a = row[g].strip()
                acc_index[a] = (row.get("Species", "").strip(), row.get("Voucher", "").strip())
                acc_index[a.split('.')[0]] = acc_index[a]

        for acc, seq in acc_to_seq.items():
            if not seq:
                fetch_summary[g].append(f"[WARN] Empty seq {acc}")
                continue
            meta = acc_index.get(acc) or acc_index.get(acc.split('.')[0])
            species, voucher = meta if meta else ("unknown", "")
            if not meta:
                fetch_summary[g].append(f"[WARN] No mapping for {acc}")
            gene_to_records[g].append((species, acc, voucher, seq))
            fetch_summary[g].append(f"[OK] {g}: {species} | {voucher} | {acc}")

    # Voucher-based searches for missing accessions (one query per species+gene, then filter vouchers locally)
    xml_dir = outdir / "xml_cache"
    species_gene_xml_cache = {}
    for species, voucher, g in voucher_requests:
        if not species:
            fetch_summary[g].append(f"[SKIP] Blank species for {g} | {voucher}")
            continue
        marker = GENE_MARKERS[g]
        key = (species, g)
        if key not in species_gene_xml_cache:
            safe_species = re.sub(r"[^\w\s-]", "", species).replace(" ", "_")
            xml_filename = f"{g}_{safe_species}.xml"
            xml_path = xml_dir / xml_filename
            query = f'"{species}"[ORGN] AND {marker} AND 0:8000[SLEN]'
            log(f"[INFO] Voucher search XML fetch: species={species} gene={g}")
            species_gene_xml_cache[key] = esearch_xml(query, save_path=xml_path)
        xml = species_gene_xml_cache[key]
        hits = parse_xml_for_voucher(xml, voucher, f"voucher filter {species} {g}")
        if not hits:
            if args.strict_voucher:
                fetch_summary[g].append(f"[FAIL] No voucher match for {species} | {voucher} | {g}")
                continue
            # fallback: first sequence in species+gene XML
            log(f"[WARN] No exact voucher match for {species} | {voucher} | {g}; using first sequence as fallback. Use --strict-voucher to skip instead.")
            if not valid_xml(xml):
                fetch_summary[g].append(f"[FAIL] Empty XML response for {species} | {voucher} | {g}")
                continue
            try:
                root = ET.fromstring(xml)
                insd = root.find(".//INSDSeq")
                if insd is None:
                    fetch_summary[g].append(f"[FAIL] No sequence found for {species} | {voucher} | {g}")
                    continue
                acc = insd.findtext("INSDSeq_accession-version") or insd.findtext("INSDSeq_primary-accession") or ""
                seq = insd.findtext("INSDSeq_sequence") or ""
                if acc and seq:
                    hits = [(acc, species, voucher, seq)]
            except ET.ParseError:
                fetch_summary[g].append(f"[FAIL] XML parse error for {species} | {voucher} | {g}")
                continue
        if hits:
            acc, org, vval, seq = hits[0]
            gene_to_records[g].append((species or org, acc, vval or voucher, seq))
            fetch_summary[g].append(f"[OK] {g}: {species or org} | {vval or voucher} | {acc}")

    # Optional: include outgroup references
    if args.include_outgroups:
        log(f"[INFO] Including outgroups: {', '.join(args.outgroups)}")
        for g in genes:
            marker = GENE_MARKERS[g]
            for org_name in args.outgroups:
                q = f"\"{org_name}\"[ORGN] AND {marker} AND 0:8000[SLEN]"
                acc, seq = esearch_first_fasta(q)
                if acc and seq:
                    gene_to_records[g].append((org_name, acc, "Outgroup", seq))
                    fetch_summary[g].append(f"[OUTGROUP-OK] {g}: {org_name} | Outgroup | {acc}")
                else:
                    fetch_summary[g].append(f"[OUTGROUP-FAIL] No {g} for {org_name}")

    # Optional: audit vouchers with missing ITS accession
    if args.audit_its:
        audit_dir.mkdir(parents=True, exist_ok=True)
        audit_path = audit_dir / "its_voucher_audit.tsv"
        with open(audit_path, "w", encoding="utf-8") as af:
            af.write("Species\tVoucher\tFound\tAccession\n")
            for row in rows:
                species = (row.get("Species", "") or "").strip()
                voucher = (row.get("Voucher", "") or "").strip()
                its_acc = (row.get("ITS", "") or "").strip()
                if species and voucher and not its_acc:
                    marker = GENE_MARKERS["ITS"]
                    vq = re.sub(r"[\s]+", " ", voucher).strip()
                    q_parts = [f"\"{species}\"[ORGN]", marker, "0:8000[SLEN]"]
                    if vq:
                        q_parts.append(f"(\"{vq}\"[All Fields] OR {vq}[All Fields])")
                    query = " AND ".join(q_parts)
                    xml = esearch_xml(query)
                    found_acc = ""
                    try:
                        root = ET.fromstring(xml)
                        for insd in root.findall(".//INSDSeq"):
                            acc = insd.findtext("INSDSeq_accession-version") or insd.findtext("INSDSeq_primary-accession") or ""
                            # optional: verify voucher in qualifiers
                            match_hits = parse_xml_for_voucher(xml, voucher)
                            if match_hits:
                                found_acc = match_hits[0][0]
                                break
                            if acc:
                                found_acc = acc
                                break
                    except ET.ParseError:
                        pass
                    af.write(f"{species}\t{voucher}\t{'YES' if found_acc else 'NO'}\t{found_acc}\n")
        log(f"[INFO] Wrote ITS voucher audit to {audit_path}")

    # Deduplicate per gene by accession AND by voucher (keep longest sequence)
    # This prevents multiple accessions for the same voucher from creating alignment issues
    for g in list(gene_to_records.keys()):
        seen_acc = set()
        voucher_best = {}  # voucher_key -> (species, acc, voucher, seq)
        
        for rec in gene_to_records[g]:
            species, acc, voucher, seq = rec
            voucher_key = f"{norm(species)}_{norm(voucher)}"
            
            # Skip if we've seen this exact accession
            if acc in seen_acc:
                log(f"[DEDUP] Skipping duplicate accession {acc} for {g}")
                continue
            
            seen_acc.add(acc)
            
            # For voucher deduplication, keep the longest sequence
            if voucher_key in voucher_best:
                existing_seq = voucher_best[voucher_key][3]
                if len(seq) > len(existing_seq):
                    old_acc = voucher_best[voucher_key][1]
                    log(f"[DEDUP] Replacing {g} for {species}|{voucher}: {old_acc} ({len(existing_seq)}bp) -> {acc} ({len(seq)}bp) [longer]")
                    voucher_best[voucher_key] = rec
                else:
                    log(f"[DEDUP] Keeping existing {g} for {species}|{voucher}: acc={voucher_best[voucher_key][1]} ({len(existing_seq)}bp) over {acc} ({len(seq)}bp)")
            else:
                voucher_best[voucher_key] = rec
        
        gene_to_records[g] = list(voucher_best.values())
        log(f"[INFO] {g}: {len(gene_to_records[g])} unique sequences after deduplication")

    # Enforce strict voucher consistency: if --strict-voucher is set, ensure one voucher per species
    if args.strict_voucher:
        # Determine the primary voucher per species (first encountered in input CSV with a non-empty voucher)
        species_primary_voucher = {}
        for row in rows:
            species = (row.get("Species", "") or "").strip()
            voucher = (row.get("Voucher", "") or "").strip()
            if species and voucher and species not in species_primary_voucher:
                species_primary_voucher[species] = voucher

        log(f"[STRICT] Enforcing single voucher per species: {len(species_primary_voucher)} species with primary vouchers")

        def is_primary(species, voucher):
            pv = species_primary_voucher.get(species)
            if not pv:
                return True  # if species had no voucher in CSV, allow
            return norm(pv) == norm(voucher)

        removed = 0
        for g in list(gene_to_records.keys()):
            filtered = []
            for (species, acc, voucher, seq) in gene_to_records[g]:
                # Allow outgroups and empty voucher entries
                if voucher == "Outgroup" or voucher == "NoVoucher":
                    filtered.append((species, acc, voucher, seq))
                    continue
                if is_primary(species, voucher):
                    filtered.append((species, acc, voucher, seq))
                else:
                    removed += 1
                    log(f"[STRICT-SKIP] Dropping {g} accession {acc} for species {species} (voucher {voucher}) != primary voucher {species_primary_voucher.get(species)}")
            gene_to_records[g] = filtered
            log(f"[STRICT] {g}: {len(gene_to_records[g])} sequences retained after voucher consistency filter")
        log(f"[STRICT] Total sequences removed due to voucher inconsistency: {removed}")

    # Write audit for rows not properly downloaded
    def write_download_audit(outdir, rows, genes, gene_to_records):
        audit_dir = outdir / "audits"
        audit_dir.mkdir(parents=True, exist_ok=True)
        audit_path = audit_dir / "download_audit.tsv"
        # Build lookup: which (species,voucher) have a record for each gene
        success = set()  # tuples (norm_species, norm_voucher, gene)
        accs_per_gene = {g: set() for g in genes}
        for g in genes:
            for (species, acc, voucher, seq) in gene_to_records.get(g, []):
                ns = norm(species)
                nv = norm(voucher)
                success.add((ns, nv, g))
                accs_per_gene[g].add(acc)
                accs_per_gene[g].add(acc.split('.')[0])
        with open(audit_path, "w", encoding="utf-8") as af:
            af.write("Species\tVoucher\tGene\tReason\n")
            for row in rows:
                species = (row.get("Species", "") or "").strip()
                voucher = (row.get("Voucher", "") or "").strip()
                ns, nv = norm(species), norm(voucher)
                for g in genes:
                    if (ns, nv, g) in success:
                        continue
                    acc = (row.get(g, "") or "").strip()
                    if acc:
                        # accession present in CSV but not returned
                        reason = "accession not returned"
                        if acc in accs_per_gene.get(g, set()) or acc.split('.')[0] in accs_per_gene.get(g, set()):
                            # Should not happen given success check, but safe-guard
                            continue
                    else:
                        reason = "no voucher match"
                    af.write(f"{species}\t{voucher}\t{g}\t{reason}\n")
        log(f"[INFO] Wrote download audit to {audit_path}")
        return audit_path

    write_download_audit(outdir, rows, genes, gene_to_records)
    write_gene_fastas(outdir, gene_to_records)
    log(f"[DONE] Wrote per-gene FASTAs to {outdir}")
    # Print summary
    for g in genes:
        sys.stderr.write(f"\nSummary for {g}:\n")
        for msg in fetch_summary[g]:
            sys.stderr.write(msg + "\n")
        sys.stderr.write(f"Total {g} sequences: {len(gene_to_records[g])}\n")


if __name__ == "__main__":
    main()