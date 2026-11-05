# GBK (GenBank Introduction)
# https://rosalind.info/problems/gbk/
#
# What GenBank is:
# GenBank is the NIH's open-access database of all publicly available DNA
# sequences. Founded by NCBI in 1982, it has grown exponentially — doubling
# roughly every 18 months — and now holds over 143 billion nucleobases.
# Every sequence has a unique identifier, and records are annotated with
# information about the organism, gene, publication date, and more.
#
# Why GenBank exists in real bioinformatics workflows:
# GenBank is the central repository where researchers deposit sequenced
# DNA from every species on Earth. Before designing primers, running a BLAST
# search, or studying a gene, scientists check GenBank to see if the
# sequence has already been published. It's the starting point for almost
# every sequence analysis project.
#
# What this solution does:
# It uses NCBI's Entrez system (via BioPython's Bio.Entrez module) to search
# the Nucleotide database for entries matching a given genus name and
# publication date range, then returns the count of matching records.


from Bio import Entrez

# Required by NCBI: identify yourself so they can contact you if needed.
Entrez.email = 'mahimmazidul@users.noreply.github.com'


def count_genbank_entries(genus, start_date, end_date):
    """Count Nucleotide GenBank entries for a genus within a date range.

    Uses Entrez esearch to query the nucleotide database. The search term
    combines the genus name (restricted to the Organism field) with a
    publication date range (PDAT = publication date).

    Args:
        genus: The genus name (e.g., 'Anthoxanthum').
        start_date: Start date in YYYY/M/D format.
        end_date: End date in YYYY/M/D format.

    Returns:
        The number of matching GenBank entries as an integer.
    """
    term = f'"{genus}"[Organism] AND {start_date}:{end_date}[PDAT]'
    handle = Entrez.esearch(db='nucleotide', term=term)
    record = Entrez.read(handle)
    handle.close()
    return int(record['Count'])


def main():
    import sys
    import os

    input_file = 'gbk_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().splitlines()
    else:
        lines = sys.stdin.read().splitlines()

    # Filter out empty lines
    lines = [line.strip() for line in lines if line.strip()]
    genus = lines[0]
    start_date = lines[1]
    end_date = lines[2]

    count = count_genbank_entries(genus, start_date, end_date)
    print(count)


if __name__ == '__main__':
    main()
