# BFIL (Base Filtration by Quality)
# https://rosalind.info/problems/bfil/
#
# What base filtration/trimming is:
# Sequencing quality typically declines toward the ends of reads. Rather
# than discarding an entire read (as in read filtration), base trimming
# removes only the low-quality bases from each end, preserving the
# high-quality middle portion. This saves valuable data for downstream
# analysis.
#
# Why trimming matters:
# Low-quality bases at read ends can cause false variant calls and
# misassemblies. Trimming them improves assembly contiguity and variant
# calling accuracy. Tools like Trimmomatic and FASTQ Quality Trimmer use
# LEADING/TRAILING or sliding window approaches to trim low-quality bases.
#
# What this solution does:
# It reads a quality threshold q and a set of FASTQ entries. For each read,
# it trims leading and trailing bases with quality below q, then outputs
# the trimmed reads in FASTQ format.


from Bio import SeqIO
import sys
import os


def trim_read(sequence, qualities, threshold):
    """Trim low-quality bases from both ends of a read.

    Finds the first base from the left with quality >= threshold and
    the last base from the right with quality >= threshold.

    Args:
        sequence: The nucleotide sequence string.
        qualities: A list of Phred quality scores.
        threshold: The quality cutoff (bases below this are trimmed).

    Returns:
        A tuple of (trimmed_sequence, trimmed_qualities).
    """
    # Find left trim position
    left = 0
    while left < len(qualities) and qualities[left] < threshold:
        left += 1

    # Find right trim position
    right = len(qualities) - 1
    while right >= 0 and qualities[right] < threshold:
        right -= 1

    if left > right:
        # All bases are below threshold
        return '', []

    return sequence[left:right + 1], qualities[left:right + 1]


# (process_fastq logic is inlined in main() for clarity)


def main():
    input_file = 'bfil_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().splitlines()
    else:
        lines = sys.stdin.read().splitlines()

    lines = [line.strip() for line in lines if line.strip()]
    threshold = int(lines[0].split()[0])
    fastq_text = '\n'.join(lines[1:])

    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write(fastq_text)
        temp_path = f.name

    try:
        records = list(SeqIO.parse(temp_path, 'fastq'))
        output_lines = []

        for record in records:
            seq = str(record.seq)
            quals = record.letter_annotations['phred_quality']
            trimmed_seq, trimmed_quals = trim_read(seq, quals, threshold)

            if not trimmed_seq:
                continue

            output_lines.append(f'@{record.id}')
            output_lines.append(trimmed_seq)
            output_lines.append('+')
            qual_string = ''.join(chr(q + 33) for q in trimmed_quals)
            output_lines.append(qual_string)

        print('\n'.join(output_lines))
    finally:
        os.unlink(temp_path)


if __name__ == '__main__':
    main()
