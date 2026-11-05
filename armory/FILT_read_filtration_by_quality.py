# FILT (Read Filtration by Quality)
# https://rosalind.info/problems/filt/
#
# What read filtration is:
# Sequencing produces millions of reads, but not all are reliable. Poor
# quality reads introduce noise and errors into downstream analysis like
# assembly or variant calling. Filtering removes reads that don't meet
# quality standards, improving results and reducing computational cost.
#
# What the quality filter does:
# The FASTX Quality Filter evaluates each read based on two criteria:
# a quality threshold q and a percentage p. A read passes if at least p%
# of its bases have a Phred quality score of at least q. For example,
# with q=20 and p=90, a read passes if 90% or more of its bases have
# quality >= 20 (99% accuracy or better).
#
# Why filtration matters:
# Low-quality reads can cause false variant calls, misassemblies, and
# wasted computation. Filtering them early in the pipeline is a standard
# quality control step in NGS analysis workflows.
#
# What this solution does:
# It reads a quality threshold q, a percentage p, and a set of FASTQ
# entries. For each read, it computes the fraction of bases with quality
# >= q and counts how many reads meet the p% threshold.


from Bio import SeqIO
import sys
import os


def count_passing_reads(input_path, q_threshold, p_percent):
    """Count reads where at least p% of bases have quality >= q.

    Args:
        input_path: Path to a FASTQ file.
        q_threshold: Minimum Phred quality score for a base to count.
        p_percent: Minimum percentage of bases that must meet q_threshold.

    Returns:
        The number of reads that pass the filter.
    """
    records = SeqIO.parse(input_path, 'fastq')
    count = 0
    for record in records:
        qualities = record.letter_annotations['phred_quality']
        total = len(qualities)
        if total == 0:
            continue
        passing = sum(1 for q in qualities if q >= q_threshold)
        if (passing / total) * 100 >= p_percent:
            count += 1
    return count


def main():
    input_file = 'filt_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().splitlines()
    else:
        lines = sys.stdin.read().splitlines()

    lines = [line.strip() for line in lines if line.strip()]
    q_threshold = int(lines[0].split()[0])
    p_percent = float(lines[0].split()[1])
    fastq_text = '\n'.join(lines[1:])

    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write(fastq_text)
        temp_path = f.name

    try:
        result = count_passing_reads(temp_path, q_threshold, p_percent)
        print(result)
    finally:
        os.unlink(temp_path)


if __name__ == '__main__':
    main()
