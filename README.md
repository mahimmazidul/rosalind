# rosalind

Solving Rosalind bioinformatics problems in Python, one by one, starting from the basics. Each solution includes a plain explanation of the biology behind the problem.

## What Rosalind is

Rosalind is a platform for learning bioinformatics through problem solving. It presents short biological questions — mostly about DNA, RNA, proteins, and genetics — and asks you to write a program that computes the answer. It's a way to pick up both programming and computational biology at the same time.

## Problems solved

| Problem | File | Status |
|---------|------|--------|
| DNA (Counting DNA Nucleotides) | `stronghold/DNA_counting_nucleotides.py` | Solved |
| RNA (Transcribing DNA into RNA) | `stronghold/RNA_transcription.py` | Solved |
| REVC (Complementing a Strand of DNA) | `stronghold/REVC_complement_strand.py` | Solved |
| GC (Computing GC Content) | `stronghold/GC_content.py` | Solved |
| HAMM (Counting Point Mutations) | `stronghold/HAMM_point_mutations.py` | Solved |

## How to run a solution

Each problem is a standalone Python script. Run it with Python 3:

```bash
python stronghold/DNA_counting_nucleotides.py
```

If the script reads from `input.txt`, place your input there. Otherwise it will read from standard input.

## About the code

Every solution file starts with a comment block explaining the biology behind the problem — what the molecules are, why the computation matters, and a link to the original Rosalind problem.
