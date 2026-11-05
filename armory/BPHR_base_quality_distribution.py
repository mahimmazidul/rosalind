# BPHR (Base Quality Distribution)
# https://rosalind.info/problems/bphr/
#
# What per-base quality distribution is:
# Sequencing quality varies along the length of a read — typically it's
# highest at the beginning and declines toward the end. Per-base quality
# distribution calculates the average quality at each position across all
# reads, revealing where the sequencer is least reliable.
#
# Why this matters in real bioinformatics:
# Researchers use this to decide where to trim reads. If quality drops
# below an acceptable threshold at position 25, they might trim all reads
# to 24 bases. This improves downstream analysis by removing error-prone
# bases. FastQC's "Per Base Sequence Quality" module is a standard QC step.
#
# What this solution does:
# It reads a FASTQ file and a quality threshold. For each base position
# across all reads, it computes the mean Phred quality score. It counts
# how many positions have a mean quality below the threshold.


from Bio import SeqIO
import sys
import os


def count_low_quality_positions(input_path, threshold):
    """Count positions where mean base quality falls below threshold.

    Args:
        input_path: Path to the input file. First line is the threshold,
                    followed by FASTQ records.
        threshold: The quality threshold (integer).

    Returns:
        The number of positions where mean quality is below threshold.
    """
    with open(input_path, 'r') as f:
        content = f.read()

    lines = content.splitlines()
    threshold = int(lines[0].strip())
    fastq_text = '\n'.join(lines[1:])

    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write(fastq_text)
        temp_path = f.name

    try:
        records = list(SeqIO.parse(temp_path, 'fastq'))
        if not records:
            return 0

        read_length = len(records[0].seq)
        position_sums = [0] * read_length
        position_counts = [0] * read_length

        for record in records:
            qualities = record.letter_annotations['phred_quality']
            for i, q in enumerate(qualities):
                if i < read_length:
                    position_sums[i] += q
                    position_counts[i] += 1

        low_quality_positions = 0
        for i in range(read_length):
            if position_counts[i] > 0:
                mean_quality = position_sums[i] / position_counts[i]
                if mean_quality < threshold:
                    low_quality_positions += 1

        return low_quality_positions
    finally:
        os.unlink(temp_path)


def main():
    input_file = 'bphr_input.txt'
    if os.path.exists(input_file):
        result = count_low_quality_positions(input_file, 0)
    else:
        lines = sys.stdin.read().splitlines()
        threshold = int(lines[0].strip())
        fastq_text = '\n'.join(lines[1:])

        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write(fastq_text)
            temp_path = f.name

        try:
            records = list(SeqIO.parse(temp_path, 'fastq'))
            if not records:
                print(0)
                return

            read_length = len(records[0].seq)
            position_sums = [0] * read_length
            position_counts = [0] * read_length

            for record in records:
                qualities = record.letter_annotations['phred_quality']
                for i, q in enumerate(qualities):
                    if i < read_length:
                        position_sums[i] += q
                        position_counts[i] += 1

            low_quality_positions = 0
            for i in range(read_length):
                if position_counts[i] > 0:
                    mean_quality = position_sums[i] / position_counts[i]
                    if mean_quality < threshold:
                        low_quality_positions += 1

            print(low_quality_positions)
        finally:
            os.unlink(temp_path)
        return

    # When reading from file, re-parse to get threshold
    with open(input_file, 'r') as f:
        lines = f.read().splitlines()
    threshold = int(lines[0].strip())
    result = count_low_quality_positions(input_file, threshold)
    print(result)


if __name__ == '__main__':
    main()
