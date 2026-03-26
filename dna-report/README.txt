# DNA Sequence Deep Analysis Skill (dna-report)

A powerful, automated DNA sequence analysis tool built for molecular biologists and researchers. This tool takes a simple DNA FASTA file, performs a deep analysis, and generates a structured, publication-ready PDF and Markdown report.

## Features
- **Basic Properties Calculation**: Computes Length, GC Content, and accurate Molecular Weights for dsDNA, ssDNA, and RNA.
- **Restriction Enzyme Mapping**: Automatically scans the sequence for 10 common 6-bp restriction enzyme cut sites (e.g., EcoRI, BamHI, NdeI) and lists their precise cut positions.
- **Homology Search (NCBI BLASTN)**: Submits the sequence to the NCBI `nt` database via the official REST API. Retrieves and parses the top 5 homology hits.
- **AI Genomic Foundation Model (Evo 2)**: Integrates an action portal for Evo 2 to analyze perplexity and sequence logo conservation.
- **AI-Assisted Functional Prediction**: Synthesizes a structured summary, functional prediction, and dynamically generates a highly relevant PubMed literature search link based on the top BLAST hit.

## Files Structure
- `dna_analyzer.py`: The core execution script.
- `input_dna.fasta`: The input file where you place your target DNA sequence.
- `evo2_actb_example.png`: An example image used in the PDF report for the Evo 2 section.
- `requirements.txt`: Python dependencies required to run the script.
- `README.txt`: This documentation file.

## Setup & Installation
1. Ensure you have Python 3.8+ installed.
2. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage
1. Open the `input_dna.fasta` file and paste your target DNA or transcript sequence in standard FASTA format.
2. Run the analyzer script:
   ```bash
   python dna_analyzer.py
   ```
3. A new timestamped folder will be created inside a `dna_analysis_runs` directory in your current working directory.
4. Open the generated `_report.pdf` or `_report.md` inside that folder to view your deep analysis results.

## Output Structure
The generated report includes 5 key sections:
1. Input Sequence
2. Basic Properties (MW, GC%, Restriction Sites)
3. AI Genomic Foundation Model (Evo 2)
4. Homology Analysis (BLASTN) - Top 5 Hits
5. AI-Assisted Function Prediction
