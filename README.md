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
| IPRB (Mendel's First Law) | `stronghold/IPRB_mendels_first_law.py` | Solved |
| PROT (Translating RNA into Protein) | `stronghold/PROT_translating_rna.py` | Solved |
| SUBS (Finding a Motif in DNA) | `stronghold/SUBS_finding_a_motif.py` | Solved |
| CONS (Consensus and Profile) | `stronghold/CONS_consensus_profile.py` | Solved |
| FIBD (Mortal Fibonacci Rabbits) | `stronghold/FIBD_mortal_fibonacci_rabbits.py` | Solved |
| GRPH (Overlap Graphs) | `stronghold/GRPH_overlap_graph.py` | Solved |
| IEV (Calculating Expected Offspring) | `stronghold/IEV_calculating_expected_offspring.py` | Solved |
| LCSM (Finding a Shared Motif) | `stronghold/LCSM_finding_a_shared_motif.py` | Solved |
| LIA (Independent Alleles) | `stronghold/LIA_independent_alleles.py` | Solved |
| MPRT (Finding a Protein Motif) | `stronghold/MPRT_finding_a_protein_motif.py` | Solved |
| MRNA (Inferring mRNA from Protein) | `stronghold/MRNA_inferring_mrna_from_protein.py` | Solved |
| ORF (Open Reading Frames) | `stronghold/ORF_finding_genes_with_orfs.py` | Solved |
| PERM (Enumerating Gene Orders) | `stronghold/PERM_enumerating_gene_orders.py` | Solved |
| PRTM (Calculating Protein Mass) | `stronghold/PRTM_calculating_protein_mass.py` | Solved |
| REVP (Locating Restriction Sites) | `stronghold/REVP_locating_restriction_sites.py` | Solved |
| SPLC (RNA Splicing) | `stronghold/SPLC_rna_splicing.py` | Solved |

## Track Progress

| Track | Folder | Solved | Total |
|-------|--------|--------|-------|
| Bioinformatics Stronghold | `stronghold/` | 21 | 21 |
| Bioinformatics Armory | `armory/` | 15 | 15 |
| Algorithmic Heights | `algorithmic_heights/` | 34 | 34 |
| Bioinformatics Textbook | `textbook/` | 20 | 132 |

## How to run a solution

Each problem is a standalone Python script. Run it with Python 3:

```bash
python stronghold/DNA_counting_nucleotides.py
```

If the script reads from `input.txt`, place your input there. Otherwise it will read from standard input.

## About the code

Every solution file starts with a comment block explaining the biology behind the problem — what the molecules are, why the computation matters, and a link to the original Rosalind problem.
