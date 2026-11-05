# GC (Computing GC Content)
# https://rosalind.info/problems/gc/
#
# What GC content is and why it matters:
# GC content is the percentage of nitrogenous bases in a DNA sequence that are
# either guanine (G) or cytosine (C). It matters because of hydrogen bonding:
# G and C are held together by three hydrogen bonds, while A and T form only
# two. That extra bond means GC-rich DNA is thermally more stable — it takes
# more heat to break the strands apart (denature them). Organisms living in
# hot environments often have higher GC content in their DNA as an adaptation.
#
# How GC content varies across organisms and genomes:
# GC content isn't uniform across a genome. In vertebrates, DNA forms
# isochores — long regions with relatively homogeneous GC content that can
# range from 30% to over 60%. Different species have different baseline GC
# contents too; for example, the bacterium Streptomyces coelicolor has a GC
# content above 70%, while Plasmodium falciparum (the malaria parasite) is
# down around 20%. Even within a single genome, coding regions often have
# higher GC content than non-coding regions.
#
# Practical uses:
# - Melting temperature: GC content is used to estimate the melting
#   temperature (Tm) of a DNA duplex, which is critical for PCR primer
#   design and hybridization experiments.
# - Gene prediction: Coding sequences tend to have higher GC content than
#   surrounding non-coding DNA, so GC content helps identify likely genes.
# - Species identification and phylogenetics: Comparing GC content across
#   species can reveal evolutionary relationships and help identify
#   unknown DNA samples.
#
# Why this was one of the first things measured:
# When DNA sequencing first became possible in the 1970s, GC content was
# among the earliest metrics computed because it requires only base counting,
# not full sequence alignment or annotation. It gave scientists a quick,
# cheap way to compare genomes before the era of whole-genome sequencing.
#
# This is Rosalind problem #4. The input is in FASTA format — multiple
# labeled DNA sequences. The task: find which sequence has the highest GC
# content and report its label and percentage.


def parse_fasta(text):
    """Parse a FASTA-formatted string into a dictionary of {label: sequence}.

    FASTA format: each sequence starts with a header line beginning with '>',
    followed by one or more lines of sequence data. This function joins all
    sequence lines for each label into a single string.

    Args:
        text: A string in FASTA format.

    Returns:
        A dictionary mapping sequence labels (without the '>') to their
        concatenated DNA sequences.
    """
    sequences = {}
    current_label = None
    current_parts = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            # Save the previous sequence before starting a new one
            if current_label is not None:
                sequences[current_label] = ''.join(current_parts)
            current_label = line[1:].strip()
            current_parts = []
        else:
            current_parts.append(line)

    # Don't forget the last sequence
    if current_label is not None:
        sequences[current_label] = ''.join(current_parts)

    return sequences


def gc_content(dna):
    """Calculate the GC content of a DNA string as a percentage.

    Args:
        dna: A string of 'A', 'C', 'G', 'T' characters.

    Returns:
        The GC content as a float percentage (0 to 100).
    """
    if not dna:
        return 0.0
    gc_count = sum(1 for base in dna if base in ('G', 'C'))
    return (gc_count / len(dna)) * 100.0


def main():
    import sys
    import os

    input_file = 'gc_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    sequences = parse_fasta(text)

    best_label = None
    best_gc = -1.0

    for label, sequence in sequences.items():
        gc = gc_content(sequence)
        if gc > best_gc:
            best_gc = gc
            best_label = label

    print(best_label)
    print(f"{best_gc:.6f}%")


if __name__ == '__main__':
    main()
