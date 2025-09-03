Thac's NMD Classification Pipeline
==================================

Contents
--------
1) Overview
2) Features
3) Requirements
4) Installation
5) Inputs
6) Outputs
7) NMD rules & coordinate systems
8) Usage synopsis
9) Typical workflows & examples
10) Column dictionaries

----------------------------------
1) Overview
----------------------------------
The pipeline performs three stages:

[1] Exon extraction — parse a GTF and emit an exon table (one row per exon, with transcript_id, chr, start, end, strand).  
[2] ORF discovery — scan each transcript sequence for ATG…STOP ORFs in all three frames (on the forward orientation of the provided sequence). Optionally filter/align to CPC2 putative peptides and/or validate/replace with Ensembl CDS.  
[3] NMD classification — map the PTC (stop codon) to transcript and genomic positions, compute distances to the last exon junction, and apply a 4‑rule heuristic to label each transcript as NMD-evasive or NMD-sensitive.

----------------------------------
2) Features
----------------------------------
- Accepts either de novo or reference inputs; Ensembl files are supported via dedicated flags.
- Exon table is produced from the supplied GTF, requiring standard attributes with `transcript_id`.
- ORFs: finds all candidate ORFs (start at 'M', stops at '*') per frame; by default keeps the **longest** ORF per transcript unless a CPC2 file is provided.
- Optional CPC2 integration: matches predicted AA sequences to `putative_peptide` and keeps only exact matches.
- Optional validation with Ensembl CDS: if a transcript’s CDS is available, the ORF’s NT/AA are replaced by CDS (and a `Match_CDS` flag indicates exact match vs replacement).
- NMD rules computed from exon structure and ORF endpoints; outputs interpretable distances and rule flags.

----------------------------------
3) Requirements
----------------------------------
- Python ≥ 3.8
- Packages: `pandas`, `biopython`
  - Install with: `pip install pandas biopython`

Optional:
- A CPC2 TSV output containing columns `#ID` and `putative_peptide`.
- Ensembl cDNA FASTA, GTF, and CDS FASTA for validation checks.

----------------------------------
4) Installation
----------------------------------
No special installation is required beyond Python and the packages above. Save the script (e.g., `nmd_pipeline.py`) and run it with `python nmd_pipeline.py <args>`.

----------------------------------
5) Inputs
----------------------------------
Exactly one transcript FASTA and one GTF must be provided, either via generic flags or Ensembl-specific flags. Ensembl flags take precedence if both are set.

- `-f/--fasta`           Transcript FASTA (spliced, 5'→3' orientation; typically cDNA-like).  
- `-g/--gtf`             Transcript GTF with `exon` features and `transcript_id` attributes.  
- `-c/--cpc2`            Optional CPC2 output TSV (columns `#ID`, `putative_peptide` at minimum).  
- `-fg/--cdna_fasta`     Ensembl cDNA FASTA. If set, overrides `-f`.  
- `-gg/--ensembl_gtf`    Ensembl GTF. If set, overrides `-g`.  
- `-cg/--cds_fasta`      Ensembl CDS FASTA for ORF validation/replacement. Optional.

Outputs (exactly three paths required via `-o/--outputs`):
- `orf_out`  (TSV)   — per-transcript ORF table (possibly CPC2-filtered and/or CDS-validated).
- `exon_out` (TSV)   — exon table derived from the GTF.  
- `nmd_out`  (CSV)   — final NMD calls and supporting metrics.

----------------------------------
6) Outputs
----------------------------------
[1] ORF table (TSV): One row per transcript after optional filters/validation. Includes frame, transcript-based start/stop, AA/NT sequences, and (if CDS provided) `Match_CDS` boolean.
[2] Exon table (TSV): One row per exon with transcript ID, chromosome, exon start, exon end, and strand.
[3] NMD table (CSV): One row per transcript with genomic and transcript coordinates of the stop codon, distances to last exon junction, rule flags, and final `NMD_Status` label.

----------------------------------
7) NMD rules & coordinate systems
----------------------------------
Coordinate systems:
- Transcript coordinates are **1-based, inclusive**.
- Genomic coordinates are reported based on exon mapping using the transcript’s strand from the GTF.
- The FASTA is assumed to be the spliced transcript in 5'→3'. The ORF finder scans this orientation in frames 1–3.

Rules (applied to the detected STOP position; if any rule is true → **NMD-evasive**, else **NMD-sensitive**):
- **Rule_PTC_<150nt**: Distance from start codon to STOP < 150 nt.
- **Rule_Exon>407nt**: STOP lies in an exon with length ≥ 407 nt.
- **Rule_<55nt_to_Junction**: STOP is < 55 nt upstream of the last exon–exon junction.
- **Rule_Last_Exon**: STOP is in the last exon.

Thresholds (150/407/55) are conventional heuristics and can be edited in the code if needed.

----------------------------------
8) Usage synopsis
----------------------------------
```
python nmd_pipeline.py \
  [-f FASTA | -fg ENSEMBL_cDNA_FASTA] \
  [-g GTF   | -gg ENSEMBL_GTF] \
  [-c CPC2_TSV] \
  [-cg ENSEMBL_CDS_FASTA] \
  -o ORF.tsv EXON.tsv NMD.csv
```

Precedence:
- If `-fg` is supplied, it is used instead of `-f`.
- If `-gg` is supplied, it is used instead of `-g`.

----------------------------------
9) Typical workflows & examples
----------------------------------
A) De novo transcripts (no CPC2, no Ensembl validation)
```
python nmd_pipeline.py \
  -f transcripts.fa \
  -g transcripts.gtf \
  -o orf.tsv exon.tsv nmd.csv
```

B) With CPC2 filtering (keep ORFs matching CPC2 `putative_peptide` exactly)
```
python nmd_pipeline.py \
  -f transcripts.fa \
  -g transcripts.gtf \
  -c cpc2_output.tsv \
  -o orf.tsv exon.tsv nmd.csv
```

C) Ensembl reference, with CDS validation/replacement
```
python nmd_pipeline.py \
  -fg Homo_sapiens.GRCh38.cdna.all.fa.gz \
  -gg Homo_sapiens.GRCh38.111.gtf.gz \
  -cg Homo_sapiens.GRCh38.cds.all.fa.gz \
  -o orf.tsv exon.tsv nmd.csv
```
In (C), if a transcript’s CDS is available:
- If ORF NT equals CDS NT → `Match_CDS = True` (no change).  
- Else the ORF NT/AA are **replaced** by the Ensembl CDS and `Match_CDS = False`.

----------------------------------
10) Column dictionaries
----------------------------------
Exon table (TSV: `exon_out`):
- `Transcript_ID`  — transcript identifier from GTF `transcript_id` attribute (version-stripped).  
- `Chr`            — chromosome name (string).  
- `Exon_Start`     — exon start (1-based, inclusive).  
- `Exon_End`       — exon end (1-based, inclusive).  
- `Strand`         — `+` or `-`.

ORF table (TSV: `orf_out`):
- `Transcript_ID`  — version-stripped FASTA ID (matches GTF transcript_id when possible).  
- `Strand`         — always `+` (the ORF scan is on the provided 5'→3' transcript sequence).  
- `Frame`          — 1, 2, or 3 (offset of the ORF within the transcript).  
- `Start`          — ORF start in transcript coordinates (1-based).  
- `Stop`           — ORF stop in transcript coordinates (1-based; if no STOP was found, this is the sequence end).  
- `Length`         — ORF nucleotide length.  
- `AA_Seq`         — amino acid sequence (no terminal `*`).  
- `NT_Seq`         — nucleotide sequence of the ORF.  
- `Match_CDS`      — (only if `-cg` used) `True` if NT_Seq == CDS; `False` if replaced by CDS; missing otherwise.

NMD table (CSV: `nmd_out`):
- `Transcript_ID`  
- `Genomic_Start`, `Genomic_Stop` — genomic positions of the ORF start/stop mapped through exon structure.  
- `Transcript_Start`, `Transcript_Stop` — transcript positions of ORF start/stop (1-based).  
- `Chr`, `Strand`  
- `Distance_From_Start_Codon` — nt from ATG to STOP.  
- `PTC_Exon_Length` — length of the exon containing STOP.  
- `Distance_To_Last_Exon_Junction` — nt from STOP to the last exon–exon junction upstream in the transcript.  
- `PTC_in_Last_Exon` — boolean.  
- `AA_Seq`, `NT_Seq`  
- `NMD_Status` — `"NMD-evasive"` if any rule is true; else `"NMD-sensitive"`.  
- `Rule_PTC_<150nt`, `Rule_Exon>407nt`, `Rule_<55nt_to_Junction`, `Rule_Last_Exon` — per-rule booleans.
