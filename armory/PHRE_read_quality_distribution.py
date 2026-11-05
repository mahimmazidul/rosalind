# PHRE (Read Quality Distribution)
# https://rosalind.info/problems/phre/
#
# What Phred quality scores are:
# Sequencing machines make errors, and not all base calls are equally
# reliable. The Phred quality score Q quantifies the confidence in each
# base call: Q = -10 * log10(P), where P is the probability that the base
# is called incorrectly. Q=10 means 90% accuracy (1 in 10 chance of error),
# Q=20 means 99% (1 in 100), Q=30 means 99.9% (1 in 1000). In FASTQ files,
# the score is encoded as an ASCII character by adding an offset of 33.
#
# Why quality score distribution matters:
# A sequencing run produces millions of reads. Rather than examining each
# read individually, researchers look at the distribution of average
# quality scores across all reads. This is a key quality control metric —
# if most reads have low average quality, there may be problems with the
# sequencing run, sample preparation, or the instrument itself. Tools like
# FastQC automate this analysis.
#
# What this solution does:
# It reads a quality threshold and a set of FASTQ entries. For each read,
# it computes the average Phred quality score across all bases (using
# Biopython's letter_annotations['phred_quality']), then counts how many
# reads have an average quality below the threshold.


from Bio import SeqIO
import sys
import os


def count_low_quality_reads(input_path, threshold):
    """Count reads whose average Phred quality is below the threshold.

    Args:
        input_path: Path to the input file. First line is the threshold,
                    followed by FASTQ records.
        threshold: The quality threshold (integer).

    Returns:
        The number of reads with average quality below the threshold.
    """
    with open(input_path, 'r') as f:
        lines = f.read().splitlines()

    # First line is the threshold
    threshold = int(lines[0].strip())

    # The rest is FASTQ data
    fastq_text = '\n'.join(lines[1:])

    # Write to a temp file for SeqIO to parse
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write(fastq_text)
        temp_path = f.name

    try:
        records = SeqIO.parse(temp_path, 'fastq')
        count = 0
        for record in records:
            qualities = record.letter_annotations['phred_quality']
            avg_quality = sum(qualities) / len(qualities)
            if avg_quality < threshold:
                count += 1
        return count
    finally:
        os.unlink(temp_path)


def main():
    input_file = 'phre_input.txt'
    if os.path.exists(input_file):
        # Read threshold from first line, FASTQ from rest
        with open(input_file, 'r') as f:
            content = f.read()

        lines = content.splitlines()
        threshold = int(lines[0].strip())
        fastq_text = '\n'.join(lines[1:])

        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write(fastq_text)
            temp_path = f.name

        try:
            records = SeqIO.parse(temp_path, 'fastq')
            count = 0
            for record in records:
                qualities = record.letter_annotations['phred_quality']
                avg_quality = sum(qualities) / len(qualities)
                if avg_quality < threshold:
                    count += 1
            print(count)
        finally:
            os.unlink(temp_path)
    else:
        # Read from stdin
        lines = sys.stdin.read().splitlines()
        threshold = int(lines[0].strip())
        fastq_text = '\n'.join(lines[1:])

        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write(fastq_text)
            temp_path = f.name

        try:
            records = SeqIO.parse(temp_path, 'fastq')
            count = 0
            for record in records:
                qualities = record.letter_annotations['phred_quality']
                avg_quality = sum(qualities) / len(qualities)
                if avg_quality < threshold:
                    count += 1
            print(count)
        finally:
            os.unlink(temp_path)


if __name__ == '__main__':
    main()
