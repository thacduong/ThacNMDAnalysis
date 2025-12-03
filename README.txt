Thac’s NMD Classification Pipeline
===================================

A unified pipeline for exon extraction, ORF discovery, optional CDS/CPC2 validation, and NMD classification.

----------------------------------------
Contents
----------------------------------------
1. Overview
2. Features
3. Requirements
4. Installation
5. Inputs
6. Outputs
7. NMD Rules & Coordinate Systems
8. Usage Synopsis
9. Typical Workflows & Examples
10. Column Dictionaries

----------------------------------------
1) Overview
----------------------------------------
The pipeline performs three integrated stages:

[1] Exon extraction
    Reads a GTF file and produces an exon structure table with:
    Transcript_ID, Chr, Exon_Start, Exon_End, Strand.

[2] ORF discovery
    Scans each transcript sequence for ATG→STOP open reading frames in all three forward frames.
    - Default: keep the longest ORF per transcript
    - Optional: filter ORFs to match CPC2 putative peptides
    - Optional: replace or validate ORFs using Ensembl CDS FASTA

[3] NMD classification
    For each transcript, computes distances relevant to NMD (exon length, PTC position, junction distance)
    and applies four heuristic rules to classify:
    NMD-evasive or NMD-sensitive.

----------------------------------------
2) Features
----------------------------------------
- Works with both de novo (StringTie/TACO) and reference (Ensembl) annotations.
- Accepts any transcript FASTA (cDNA or assembled).
- Can automatically generate transcript FASTA from a genome FASTA using pybedtools + BEDTools.
- ORF extraction finds all ATG→STOP ORFs in frames 1–3.
- Optional CPC2 filtering keeps ORFs matching CPC2 "putative_peptide".
- Optional CDS validation replaces ORFs with canonical CDS sequences and marks Match_CDS flags.
- NMD classification provides all rule flags, distances, and final status.

----------------------------------------
3) Requirements
----------------------------------------
Mandatory:
- Python 3.8 or newer
- pandas
- biopython

Install (local/HPC):
    pip install --user pandas biopython

If using genome FASTA auto-build (-r):
- pybedtools
- BEDTools (must include fastaFromBed)

Install:
    pip install --user pybedtools
    module load BEDTools

Optional:
- CPC2 output TSV (columns: "#ID", "putative_peptide")
- Ensembl CDS FASTA

----------------------------------------
4) Installation
----------------------------------------
Save the script (e.g., nmd_classify_clean.py)
Run with:

    python nmd_classify_clean.py <arguments>

----------------------------------------
5) Inputs
----------------------------------------
Required:
- -g / --gtf
  GTF annotation containing exon features with transcript_id

Plus one of the following:
- -f / --fasta
  Transcript FASTA (spliced, 5'→3')
OR
- -r / --reference_fasta
  Genome FASTA (pipeline will auto-build transcripts)

Optional:
- -c / --cpc2        CPC2 TSV file
- -cg / --cds_fasta  Ensembl CDS FASTA for validation

Deprecated but accepted:
- -fg / --cdna_fasta
- -gg / --ensembl_gtf

Output files (three required):
    -o ORF.tsv EXON.tsv NMD.csv

----------------------------------------
6) Outputs
----------------------------------------
1) ORF table (TSV)
   Includes Transcript_ID, Frame, Start/Stop, AA_Seq, NT_Seq, optional Match_CDS

2) Exon table (TSV)
   One row per exon: Transcript_ID, Chr, Exon_Start, Exon_End, Strand

3) NMD table (CSV)
   Genomic/transcript ORF coordinates, distances, rule flags, NMD_Status

----------------------------------------
7) NMD Rules & Coordinate Systems
----------------------------------------
Coordinate systems:
- Transcript coordinates: 1-based, inclusive
- Genomic mapping uses strand from GTF
- FASTA inputs assumed spliced and 5'→3'

NMD-evasion rules (any TRUE → NMD-evasive):
1. Rule_PTC_<150nt
   Stop codon <150 nt from start codon
2. Rule_Exon>407nt
   Stop lies in an exon >= 407 nt
3. Rule_<55nt_to_Junction
   Stop <55 nt upstream of last exon–exon junction
4. Rule_Last_Exon
   Stop codon in last exon

Else → NMD-sensitive.

----------------------------------------
8) Usage Synopsis
----------------------------------------
General:

    python nmd_classify_clean.py         -g GTF         (-f TRANSCRIPTS.fa | -r GENOME.fa)         [-c CPC2.tsv]         [-cg CDS.fa]         -o ORF.tsv EXON.tsv NMD.csv

Notes:
- If -f is missing but -r is provided, transcript FASTA is auto-generated.
- -fg and -gg still work but are deprecated.

----------------------------------------
9) Typical Workflows & Examples
----------------------------------------

A) Ensembl reference with CDS validation:

    python nmd_classify_clean.py         -g Homo_sapiens.GRCh38.114.gtf         -f Homo_sapiens.GRCh38.cdna.all.fa         -cg Homo_sapiens.GRCh38.cds.all.fa         -o orf.tsv exon.tsv nmd.csv

B) De novo transcriptome (StringTie/TACO):

    python nmd_classify_clean.py         -g taco_merged.gtf         -f taco_merged.fa         -o orf.tsv exon.tsv nmd.csv

C) Auto-build transcript FASTA from genome:

    python nmd_classify_clean.py         -g Homo_sapiens.GRCh38.114.gtf         -r Homo_sapiens.GRCh38.dna.primary_assembly.fa         -o orf.tsv exon.tsv nmd.csv

D) Using CPC2 filtering:

    python nmd_classify_clean.py         -g transcripts.gtf         -f transcripts.fa         -c cpc2_output.tsv         -o orf.tsv exon.tsv nmd.csv

----------------------------------------
10) Column Annotation
----------------------------------------

Exon table (exon_out):
- Transcript_ID
- Chr
- Exon_Start
- Exon_End
- Strand

ORF table (orf_out):
- Transcript_ID
- Strand (always +)
- Frame (1–3)
- Start (transcript coordinate)
- Stop (transcript coordinate)
- Length
- AA_Seq
- NT_Seq
- Match_CDS (optional)

NMD table (nmd_out):
- Transcript_ID
- Genomic_Start
- Genomic_Stop
- Transcript_Start
- Transcript_Stop
- Chr
- Strand
- Distance_From_Start_Codon
- PTC_Exon_Length
- Distance_To_Last_Exon_Junction
- PTC_in_Last_Exon
- AA_Seq
- NT_Seq
- NMD_Status
- Rule_PTC_<150nt
- Rule_Exon>407nt
- Rule_<55nt_to_Junction
- Rule_Last_Exon

----------------------------------------
End of README
----------------------------------------
