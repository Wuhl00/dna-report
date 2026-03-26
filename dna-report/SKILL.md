---
name: "dna-report"
description: "DNA sequence analysis Skill. Input a DNA FASTA to run basic property analysis (GC, MW for dsDNA/ssDNA/RNA), restriction enzyme scanning, NCBI BLASTN homology search, and generate a PDF/Markdown report with dynamic AI functional prediction. Invoke when user wants to analyze a DNA sequence."
---

# DNA Sequence Deep Analysis Skill (dna-report)

This Skill is designed for DNA sequences and turns a multi-step bioinformatics workflow into a single input action. Provide a DNA FASTA sequence and get a ready-to-read PDF report plus a Markdown report.

## Features

- **Basic properties**: Sequence length, GC content, and approximate Molecular Weight calculations optimized for different molecule types (dsDNA, ssDNA, and RNA/transcripts).
- **Restriction Enzyme Scanning**: Scans for 10 common 6-cutter restriction enzyme sites (EcoRI, BamHI, HindIII, XhoI, NotI, NdeI, NheI, NcoI, BglII, SalI) and reports their cut counts and positions.
- **AI Genomic Foundation Model (Evo 2)**: Introduces the Evo 2 portal with a direct citation link and an illustrative example (Human actB sequence analysis) featuring an ATGC sequence logo.
- **Homology Search (NCBI BLASTN)**: Performs asynchronous BLASTN search against the comprehensive NCBI `nt` (Nucleotide collection) database via REST API. Neatly parses the top 5 hits, keeping descriptions concise and formatted.
- **AI-Assisted Function Prediction**: 
  - Generates an automated Investigation Summary and Functional Prediction based on BLASTN identity scores.
  - Features dynamic PubMed literature search: automatically extracts relevant keywords from the top BLASTN hit title by stripping common stopwords (like "mRNA", "uncharacterized", "transcript") to generate a highly targeted search link.
- **Outputs**:
  - **PDF**: A beautifully formatted one-stop report with sidebar bookmarks, dynamic page layout to avoid sequence truncation, and clickable external links.
  - **Markdown**: A clean, structured text report easy to edit, copy, or share.

## Tech Stack

- **Parsing & Analysis**: `Biopython` (FASTA parsing, GC calculation)
- **Homology Search**: `NCBI BLAST REST API` (async submit/poll/fetch for blastn) + XML parsing via `Bio.Blast.NCBIXML`
- **PDF Reporting**: `fpdf` (PDF generation, dynamic layouts, custom colors) + `PyPDF2` (write sidebar bookmarks/outline)
- **Regex & String Manipulation**: Advanced `re` parsing for title cleaning, PubMed keyword extraction, and overlapping sequence scanning.

## Usage

- **Input**: A DNA sequence in standard FASTA format.
- **How to run**:
  1. Navigate to the `skills/dna-report` directory.
  2. Provide a sequence in an `input_dna.fasta` file (or pass any valid fasta file to the analyzer script).
  3. Run `python dna_analyzer.py`
- **Output location**: One output folder per run to avoid overwriting results:
  - `dna_analysis_runs/<FASTA_ID>_YYYYMMDD_HHMMSS/`

## Example

### End-to-end analysis
User input (`input_dna.fasta`):
```fasta
>Sample_DNA
ATGCGTACGTAGCTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGCTAGCTGATC
```

Run command:
```bash
python dna_analyzer.py
```

Outputs:
- `dna_analysis_runs/Sample_DNA_YYYYMMDD_HHMMSS/Sample_DNA_report.pdf`
- `dna_analysis_runs/Sample_DNA_YYYYMMDD_HHMMSS/Sample_DNA_report.md`