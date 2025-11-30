# Bioinformatics – Sequence Overlap & Assembly (IFT3295)

This repository contains my solution for a programming assignment in the **Bioinformatics** course (IFT3295, Université de Montréal).  
The goal is to reconstruct a long DNA sequence from short reads by computing optimal overlaps and building an overlap graph.

## Project overview

Given a set of short DNA reads (FASTA/FASTQ):

1. **Pairwise overlap computation**
   - Implement a dynamic programming algorithm to compute the optimal overlap score between pairs of reads (similar to edit-distance alignment, but restricted to suffix–prefix overlaps).
   - Support a configurable scoring scheme (match / mismatch / gap).

2. **Score matrix & overlap graph**
   - Build a **score matrix** where entry *(i, j)* is the best overlap score from read *i* to read *j*.
   - Construct a **directed weighted graph** whose nodes are reads and whose edges represent high-scoring overlaps.

3. **Filtering & path construction**
   - Filter edges based on score thresholds.
   - Remove redundant edges (simple transitive reduction–style heuristics).
   - Traverse the resulting graph greedily to obtain an assembly path.

4. **Sequence reconstruction**
   - Reconstruct an approximate target sequence by following the path and merging overlapping reads.

## How to run

Each Python file contains a short usage description at the top of the file.  
If you want to re-run what we did locally, you can follow these steps:

1. **Clone the repository**

    git clone https://github.com/pierre-emery/bioinfo-sequence-assembly.git
    cd bioinfo-sequence-assembly

2. **(Optional) Open in VS Code**

    - Open VS Code.  
    - Go to **File → Open Folder…** and select the cloned folder.

3. **Open a terminal and go to the project directory**

    In VS Code (or any terminal):

    cd base_python_tp1

From there, you can run each script individually.  
Below are example commands; you can adjust the parameters if you want to experiment with different scoring schemes or thresholds.

---

### Overlap alignment: `prefixe_suffixe.py`

Computes the optimal suffix–prefix overlap between two reads.

    python prefixe_suffixe.py two_reads.fastq

    # With custom scoring:
    python prefixe_suffixe.py two_reads.fastq --match 4 --mismatch -4 --gap -8

You can change `--match`, `--mismatch` and `--gap` to try other scoring values.

---

### Overlap score matrix: `matrice.py`

Builds a 20×20 matrix of overlap scores for the first 20 reads of a FASTQ file and saves it as CSV.

    python matrice.py reads.fq

    # Options:
    python matrice.py reads.fq --out matrice_20x20.csv --match 4 --mismatch -4 --gap -8

You can change:

- the `--out` filename  
- the scoring parameters (`--match`, `--mismatch`, `--gap`).

---

### Overlap graphs: `graph.py`

Constructs overlap graphs from the 20×20 score matrix and writes them as DOT files.

    python graph.py matrice_20x20.csv --threshold 80

This produces:

- `graph_1.dot` – filtered graph before reduction  
- `graph_2.dot` – reduced graph (after removing 2-cycles and applying simplification rules)

You can adjust `--threshold` to change which edges are kept in the graph.

---

### Sequence assembly: `sequence_frag.py`

Reconstructs the final sequence from the reduced graph and the reads.

    python sequence_frag.py

    # Options:
    python sequence_frag.py --reads reads.fq --dot graph_2.dot --match 4 --mismatch -4 --gap -8

You can:

- specify another FASTQ file with `--reads`  
- use a different graph file with `--dot`  
- change the scoring scheme with `--match`, `--mismatch`, `--gap`.

---

### Start codon frame: `codon_start.py`

Finds the reading frame of the start codon of a protein (`geneX.fasta`) inside a genomic region (`sequence.fasta`), assuming the genomic sequence is the coding strand.

    python codon_start.py --genome sequence.fasta --protein geneX.fasta

You can swap in your own FASTA files via `--genome` and `--protein` if you want to test other sequences.
