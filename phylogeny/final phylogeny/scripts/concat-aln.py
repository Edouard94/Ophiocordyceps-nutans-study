#! /usr/bin/env python3

import re
import glob
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Nexus import Nexus
# from Bio.SeqRecord import SeqRecord
# from Bio.Align import MultipleSeqAlignment

def usage():
    msg = '''concatenate genewise alignments in fasta format into single PHYLIP and NEXUS formatted alignments
    Usage:
    concat-aln <alignment-input-dir> <output-file-name> <DNA|protein>
          '''

    print(msg)

def process_all_alignments(alignment_output_dir):

    all_genes = [x for x in glob.glob(alignment_output_dir+"/*") if (x.endswith(".fasta") or x.endswith(".fa"))]
    all_alignments = {re.sub("/|_.*", "", gene.replace(alignment_output_dir, "")):AlignIO.read(gene, "fasta") for gene in all_genes}

    return all_alignments

def combine_alignments(alignment_dict):
    iteron = [(k,Nexus.Nexus(format(alignment_dict[k], "nexus"))) for k in alignment_dict]
    combined_nexus = Nexus.combine(iteron)
    return(combined_nexus)

def main():

    try:
        input_dir   = sys.argv[1]
        output_name = sys.argv[2]
        input_type  = sys.argv[3]
    except IndexError:
        usage()
        sys.exit()

    if len(sys.argv) < 4 or input_dir in ('-h', '--help'):
        usage()
        sys.exit()

    all_alignments = process_all_alignments(input_dir)
    [seq.annotations.update({"molecule_type":input_type}) for gene in all_alignments for seq in all_alignments[gene]]
    combined_nexus = combine_alignments(all_alignments)
    combined_nexus.write_nexus_data(output_name+".nex")
    combined_nexus.export_phylip(output_name+".phy")
    print(f"wrote concatenated alignment from directory {input_dir} to files {output_name}.phy and {output_name}.nex")

    return 0

if __name__ == "__main__":
    main()


