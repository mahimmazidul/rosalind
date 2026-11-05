# FRMT (Data Formats)
# https://rosalind.info/problems/frmt/
#
# What data formats are:
# Biological data can be stored in many file formats — FASTA, GenBank,
# NEXUS, PHYLIP, and more. Each format was designed for different purposes:
# FASTA is simple and universal, GenBank carries rich annotations, NEXUS
# and PHYLIP are used for phylogenetic analysis. Converting between formats
# is a routine task in bioinformatics pipelines.
#
# Why format conversion matters:
# Different tools expect different formats. A sequence downloaded from
# GenBank might need to be converted to FASTA before running it through a
# motif finder, or to PHYLIP before building a phylogenetic tree. Being
# able to parse and convert formats is essential for connecting tools in
# a workflow.
#
# What this solution does:
# It takes GenBank accession IDs, fetches each record from NCBI's Nucleotide
# database in FASTA format using Entrez.efetch, finds the shortest sequence,
# and prints it in FASTA format.


from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'mahimmazidul@users.noreply.github.com'


def fetch_fasta_records(access_ids):
    """Fetch GenBank records in FASTA format for the given accession IDs.

    Args:
        access_ids: A list of GenBank accession ID strings.

    Returns:
        A list of SeqRecord objects parsed from the FASTA data.
    """
    id_string = ','.join(access_ids)
    handle = Entrez.efetch(db='nucleotide', id=id_string, rettype='fasta')
    records = list(SeqIO.parse(handle, 'fasta'))
    handle.close()
    return records


def find_shortest_record(records):
    """Return the SeqRecord with the shortest sequence.

    Args:
        records: A list of SeqRecord objects.

    Returns:
        The SeqRecord with the shortest sequence length.
    """
    return min(records, key=lambda r: len(r.seq))


def main():
    import sys
    import os

    input_file = 'frmt_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            access_ids = f.read().split()
    else:
        access_ids = sys.stdin.read().split()

    records = fetch_fasta_records(access_ids)
    shortest = find_shortest_record(records)

    # Print in FASTA format. description already includes the accession ID.
    print(f'>{shortest.description}')
    # Wrap sequence at 60 characters per line (standard FASTA wrapping)
    seq_str = str(shortest.seq)
    for i in range(0, len(seq_str), 60):
        print(seq_str[i:i + 60])


if __name__ == '__main__':
    main()
