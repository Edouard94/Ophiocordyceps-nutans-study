#!/usr/bin/env python3
"""Fetch GenBank sequences by ID from a TSV and split into gene FASTA files."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import sys
import time
from typing import Dict, Iterable, List, Tuple

from Bio import Entrez, SeqIO


GENE_COLUMNS = ["nSSU rDNA", "nLSU rDNA", "TEF", "RPB1", "RPB2"]
MISSING_VALUES = {"", "-", "–", "—", "na", "n/a"}


def normalize_cell(value: str | None) -> str:
	if value is None:
		return ""
	return value.strip()


def is_missing(value: str | None) -> bool:
	if value is None:
		return True
	return value.strip().lower() in MISSING_VALUES


def split_genbank_ids(cell: str | None) -> List[str]:
	if is_missing(cell):
		return []
	raw = cell.strip()
	parts = re.split(r"[;,]\s*|\s+", raw)
	ids: List[str] = []
	for part in parts:
		token = part.strip()
		if not token or token.lower() in MISSING_VALUES:
			continue
		ids.append(token)
	return ids


def safe_header_field(value: str | None) -> str:
	value = normalize_cell(value)
	if not value:
		return "NA"
	value = re.sub(r"\s+", "_", value)
	value = re.sub(r"[^A-Za-z0-9_.-]", "", value)
	return value or "NA"


def host_abbrev(value: str | None) -> str:
	value = normalize_cell(value)
	if not value:
		return "NA"
	value = re.sub(r"[^A-Za-z0-9]", "", value)
	if not value:
		return "NA"
	return value[:3].upper()


def fetch_fasta_record(gb_id: str, retries: int = 3, delay: float = 0.34):
	last_error: Exception | None = None
	for attempt in range(1, retries + 1):
		try:
			with Entrez.efetch(
				db="nucleotide",
				id=gb_id,
				rettype="fasta",
				retmode="text",
			) as handle:
				records = list(SeqIO.parse(handle, "fasta"))
			if not records:
				raise ValueError("No records returned")
			if len(records) > 1:
				# Keep first record but allow tracing in stdout.
				print(
					f"Warning: {gb_id} returned {len(records)} records; using first",
					file=sys.stderr,
				)
			return records[0]
		except Exception as exc:  # pragma: no cover - network
			last_error = exc
			if attempt < retries:
				time.sleep(delay)
			else:
				raise
	if last_error:
		raise last_error


def build_header(row: Dict[str, str]) -> str:
	genus = safe_header_field(row.get("CurrentGenus"))
	species = safe_header_field(row.get("CurrentSpecies"))
	herbcol = safe_header_field(row.get("HerbCol"))
	voucher = safe_header_field(row.get("Voucher"))
	host = host_abbrev(row.get("Host (Higher)"))
	return f"{genus}_{species}_{herbcol}_{voucher}_{host}"


def open_gene_files(output_dir: pathlib.Path) -> Dict[str, pathlib.Path]:
	output_dir.mkdir(parents=True, exist_ok=True)
	gene_paths: Dict[str, pathlib.Path] = {}
	for gene in GENE_COLUMNS:
		safe_gene = re.sub(r"\s+", "_", gene)
		gene_paths[gene] = output_dir / f"{safe_gene}.fasta"
	return gene_paths


def write_records(
	input_path: pathlib.Path,
	output_dir: pathlib.Path,
	email: str,
	delay: float,
	failures_path: pathlib.Path,
) -> None:
	Entrez.email = email

	gene_paths = open_gene_files(output_dir)
	failure_rows: List[Tuple[str, str, str, str]] = []

	with input_path.open("r", encoding="utf-8") as handle:
		reader = csv.DictReader(handle, delimiter="\t")
		missing_headers = [
			col
			for col in [
				"CurrentGenus",
				"CurrentSpecies",
				"HerbCol",
				"Voucher",
				"Host (Higher)",
				*GENE_COLUMNS,
			]
			if col not in reader.fieldnames
		]
		if missing_headers:
			missing = ", ".join(missing_headers)
			raise ValueError(f"Missing required columns: {missing}")

		for row in reader:
			header = build_header(row)
			isolate = normalize_cell(row.get("Isolate/ID"))
			for gene in GENE_COLUMNS:
				ids = split_genbank_ids(row.get(gene))
				if not ids:
					continue
				gene_path = gene_paths[gene]
				with gene_path.open("a", encoding="utf-8") as gene_handle:
					for gb_id in ids:
						try:
							record = fetch_fasta_record(gb_id, delay=delay)
							record.id = header
							record.description = header
							SeqIO.write(record, gene_handle, "fasta")
							time.sleep(delay)
						except Exception as exc:  # pragma: no cover - network
							failure_rows.append(
								(
									isolate or "NA",
									gene,
									gb_id,
									str(exc),
								)
							)

	if failure_rows:
		with failures_path.open("w", encoding="utf-8", newline="") as fail_handle:
			writer = csv.writer(fail_handle, delimiter="\t")
			writer.writerow(["Isolate/ID", "Gene", "GenBank_ID", "Error"])
			writer.writerows(failure_rows)


def build_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(
		description="Extract GenBank sequences by ID and split into gene FASTA files."
	)
	parser.add_argument(
		"-i",
		"--input",
		type=pathlib.Path,
		default=pathlib.Path("List_of_Ophiocordyceps_strains.tsv"),
		help="Input TSV file.",
	)
	parser.add_argument(
		"-o",
		"--output-dir",
		type=pathlib.Path,
		default=pathlib.Path("fasta_by_gene"),
		help="Output directory for FASTA files.",
	)
	parser.add_argument(
		"--email",
		required=True,
		help="Email address required by NCBI Entrez.",
	)
	parser.add_argument(
		"--delay",
		type=float,
		default=0.34,
		help="Delay (seconds) between NCBI requests.",
	)
	parser.add_argument(
		"--failures",
		type=pathlib.Path,
		default=pathlib.Path("genbank_failures.tsv"),
		help="TSV report for GenBank IDs that failed.",
	)
	return parser


def main() -> int:
	parser = build_parser()
	args = parser.parse_args()

	if not args.input.exists():
		print(f"Input file not found: {args.input}", file=sys.stderr)
		return 1

	try:
		write_records(
			input_path=args.input,
			output_dir=args.output_dir,
			email=args.email,
			delay=args.delay,
			failures_path=args.failures,
		)
	except Exception as exc:
		print(f"Error: {exc}", file=sys.stderr)
		return 1

	return 0


if __name__ == "__main__":
	raise SystemExit(main())
