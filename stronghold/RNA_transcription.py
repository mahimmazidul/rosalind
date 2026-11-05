# RNA (Transcribing DNA into RNA)
# https://rosalind.info/problems/rna/
#
# What is RNA, and how is it different from DNA?
# RNA (ribonucleic acid) is a molecule similar to DNA, but with a few key
# differences. DNA stores genetic information long-term, while RNA acts more
# like a working copy — it carries instructions from DNA to the cell's
# protein-building machinery. DNA uses the sugar deoxyribose and the base
# thymine (T); RNA uses ribose and the base uracil (U) instead of thymine.
# DNA is usually double-stranded, while RNA is typically single-stranded.
#
# What is transcription?
# Transcription is the biological process where an enzyme called RNA
# polymerase reads a DNA template strand and builds a complementary mRNA
# (messenger RNA) strand. It unzips a short region of DNA, reads the DNA
# bases one by one, and assembles a matching RNA strand — except instead of
# pairing A with T as in DNA replication, it pairs A with U. So every T in
# the DNA template becomes a U in the RNA transcript. This mRNA then travels
# to the ribosome, where it's translated into protein.
#
# Why thymine becomes uracil:
# It's not just an arbitrary rule. Uracil and thymine are chemically very
# similar, but uracil lacks a methyl group that thymine has. Using uracil in
# RNA is energetically cheaper (RNA is often short-lived and rebuilt
# frequently), while thymine's extra stability makes it better suited for
# DNA, which needs to last a lifetime. So when RNA polymerase builds an RNA
# strand, it incorporates uracil wherever the DNA template calls for thymine.
#
# Why this matters:
# mRNA is the molecule that actually gets translated into protein at the
# ribosome. Transcription is therefore the first step in gene expression —
# it's how the information stored in DNA gets converted into the proteins
# that do the work in cells. Without transcription, DNA's instructions would
# stay locked away and never become functional biology.
#
# This is Rosalind problem #2 in the Stronghold track. The task: given a DNA
# string, transcribe it into RNA by replacing every T with U.


def transcribe(dna):
    """Transcribe a DNA string into RNA by replacing every T with U.

    Args:
        dna: A string consisting of the characters 'A', 'C', 'G', and 'T'.

    Returns:
        The RNA string with every 'T' replaced by 'U'.
    """
    return dna.replace('T', 'U')


def main():
    import sys
    import os

    input_file = 'rna_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            dna = f.read().strip()
    else:
        dna = sys.stdin.read().strip()

    print(transcribe(dna))


if __name__ == '__main__':
    main()
