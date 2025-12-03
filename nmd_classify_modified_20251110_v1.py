#!/usr/bin/env python3

import argparse
import os
import tempfile
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

############################################
# Argument Parser
############################################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Unified NMD classification pipeline (clean version)."
    )

    # Main, unified arguments
    parser.add_argument(
        "-g", "--gtf", required=False,
        help="GTF annotation file (custom or Ensembl)."
    )
    parser.add_argument(
        "-f", "--fasta", required=False,
        help="Transcript/cDNA FASTA file (custom or Ensembl)."
    )
    parser.add_argument(
        "-r", "--reference_fasta", required=False,
        help="Genome FASTA file used to build transcript FASTA from GTF if -f is not provided."
    )

    # Optional extras
    parser.add_argument(
        "-c", "--cpc2", required=False,
        help="Optional CPC2 output file (tab-delimited)."
    )
    parser.add_argument(
        "-cg", "--cds_fasta", required=False,
        help="Optional Ensembl CDS FASTA file for ORF validation."
    )
    parser.add_argument(
        "-o", "--outputs", nargs=3, required=True,
        metavar=("orf_out", "exon_out", "nmd_out"),
        help="Output files: ORF.tsv EXON.tsv NMD.csv"
    )

    # Backward-compatible deprecated flags
    parser.add_argument(
        "-fg", "--cdna_fasta", required=False,
        help="(Deprecated) Ensembl cDNA FASTA file. Use -f instead."
    )
    parser.add_argument(
        "-gg", "--ensembl_gtf", required=False,
        help="(Deprecated) Ensembl GTF file. Use -g instead."
    )

    return parser.parse_args()

############################################
# GTF Attribute Parsing
############################################
def parse_attributes(attr_string):
    attrs = {}
    for field in attr_string.strip().split(";"):
        if field.strip():
            key_value = field.strip().split(" ", 1)
            if len(key_value) == 2:
                key, value = key_value
                attrs[key] = value.strip('"')
    return attrs

############################################
# Step 1: Extract Exons
############################################
def extract_exons(gtf_file, exon_out):
    """
    Extract exon structure from a GTF and write a simple exon table.
    """
    with open(gtf_file) as gtf, open(exon_out, "w") as out:
        out.write("Transcript_ID\tChr\tExon_Start\tExon_End\tStrand\n")
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, _, _, start, end, _, strand, _, attributes = fields
            attrs = parse_attributes(attributes)
            tx_id = attrs.get("transcript_id")
            if tx_id:
                out.write(f"{tx_id}\t{chrom}\t{start}\t{end}\t{strand}\n")
    print(f"[✓] Exon structure saved to {exon_out}")

############################################
# Step 1b: Build transcript FASTA from exon table + genome (optional)
############################################
def build_transcript_fasta_from_exons(exon_tsv, genome_fa, out_fa):
    """
    Build cDNA-like transcript sequences by concatenating exons
    for each transcript, using a genome FASTA and an exon table.

    Requires pybedtools to be installed.
    """
    try:
        from pybedtools import BedTool
    except ImportError:
        raise ImportError(
            "pybedtools is required to build transcript FASTA from GTF. "
            "Install it or provide -f/--fasta directly."
        )

    exon_df = pd.read_csv(exon_tsv, sep="\t", dtype={"Chr": str})
    exon_df["Exon_Start"] = exon_df["Exon_Start"].astype(int)
    exon_df["Exon_End"] = exon_df["Exon_End"].astype(int)

    bed_rows = []
    # One BED row per exon, ordered in transcript 5'→3' direction
    for tx_id, group in exon_df.groupby("Transcript_ID"):
        strand = group.iloc[0]["Strand"]
        if strand == "+":
            group_sorted = group.sort_values("Exon_Start")
        else:
            # For minus strand, 5' end is at higher genomic coordinate
            group_sorted = group.sort_values("Exon_Start", ascending=False)

        for exon_idx, row in enumerate(group_sorted.itertuples()):
            chrom = row.Chr
            start0 = row.Exon_Start - 1  # BED: 0-based start
            end1 = row.Exon_End         # BED: 1-based end
            name = f"{tx_id}|{exon_idx}"
            bed_rows.append([chrom, start0, end1, name, 0, strand])

    bed = BedTool(bed_rows)

    # Extract sequences per exon; s=True respects strand
    tmp_exon_fa = out_fa + ".exons.tmp.fa"
    bed.sequence(fi=genome_fa, fo=tmp_exon_fa, name=True, s=True)

    # Concatenate exons per transcript in the order they appear in tmp_exon_fa
    exon_seqs = {}
    for record in SeqIO.parse(tmp_exon_fa, "fasta"):
        tx_id = record.id.split("|")[0]
        exon_seqs.setdefault(tx_id, []).append(str(record.seq))

    records = []
    for tx_id, exon_list in exon_seqs.items():
        full_seq = "".join(exon_list)
        records.append(SeqIO.SeqRecord(Seq(full_seq), id=tx_id, description=""))

    SeqIO.write(records, out_fa, "fasta")
    os.remove(tmp_exon_fa)
    print(f"[✓] Transcript FASTA built from GTF + genome and saved to {out_fa}")
    return out_fa

############################################
# Step 2: Find ORFs
############################################
def find_orfs_in_sequence(seq, strand, tx_id):
    orfs = []
    for frame in range(3):
        trans = seq[frame:].translate(to_stop=False)
        aa_seq = str(trans)
        for start_idx in range(len(aa_seq)):
            if aa_seq[start_idx] != "M":
                continue
            found_stop = False
            for end_idx in range(start_idx + 1, len(aa_seq)):
                if aa_seq[end_idx] == "*":
                    nt_start = frame + start_idx * 3
                    nt_end = frame + end_idx * 3 + 2
                    orfs.append({
                        "Transcript_ID": tx_id,
                        "Strand": strand,
                        "Frame": frame + 1,
                        "Start": nt_start + 1,
                        "Stop": nt_end + 1,
                        "Length": len(aa_seq[start_idx:end_idx]) * 3,
                        "AA_Seq": aa_seq[start_idx:end_idx],
                        "NT_Seq": str(seq[nt_start:nt_end + 1])
                    })
                    found_stop = True
                    break
            if not found_stop:
                nt_start = frame + start_idx * 3
                nt_end = len(seq) - 1
                orfs.append({
                    "Transcript_ID": tx_id,
                    "Strand": strand,
                    "Frame": frame + 1,
                    "Start": nt_start + 1,
                    "Stop": nt_end + 1,
                    "Length": len(aa_seq[start_idx:]) * 3,
                    "AA_Seq": aa_seq[start_idx:],
                    "NT_Seq": str(seq[nt_start:nt_end + 1])
                })
    return orfs

def find_all_orfs(fasta_path, cpc2_path=None):
    all_orfs = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        tx_id = record.id.split(".")[0]
        seq = record.seq
        orfs = find_orfs_in_sequence(seq, "+", tx_id)
        if not orfs:
            continue
        if cpc2_path:
            # Keep all for CPC2 filtering later
            all_orfs.extend(orfs)
        else:
            # Keep longest ORF per transcript
            longest = max(orfs, key=lambda x: x["Length"])
            all_orfs.append(longest)

    df = pd.DataFrame(all_orfs)

    if cpc2_path and not df.empty:
        cpc2_df = pd.read_csv(cpc2_path, sep="\t")
        cpc2_df["putative_peptide"] = cpc2_df["putative_peptide"].str.replace("*", "", regex=False)
        cpc2_df = cpc2_df.rename(columns={"#ID": "Transcript_ID"})
        df = df.merge(cpc2_df[["Transcript_ID", "putative_peptide"]], on="Transcript_ID", how="left")
        df["Match_CPC2"] = df.apply(
            lambda r: r["AA_Seq"] == r["putative_peptide"] if pd.notna(r["putative_peptide"]) else False,
            axis=1
        )
        df = df[df["Match_CPC2"]].drop(columns=["putative_peptide", "Match_CPC2"])

    return df

############################################
# Step 2b: Optional CDS validation
############################################
def validate_with_cds(orf_df, cds_fasta):
    cds_dict = {
        record.id.split(".")[0]: str(record.seq)
        for record in SeqIO.parse(cds_fasta, "fasta")
    }

    def replace_orf(row):
        tx_id = row["Transcript_ID"]
        cds_seq = cds_dict.get(tx_id)
        if cds_seq:
            if row["NT_Seq"] == cds_seq:
                row["Match_CDS"] = True
            else:
                row["NT_Seq"] = cds_seq
                row["AA_Seq"] = str(Seq(cds_seq).translate(to_stop=True))
                row["Match_CDS"] = False
        else:
            row["Match_CDS"] = False
        return row

    return orf_df.apply(replace_orf, axis=1)

############################################
# Step 3: NMD Classification
############################################
def transcript_to_genomic_coords(start, stop, exon_list, strand):
    tx_pos = 1
    genomic_start = None
    genomic_stop = None
    exon_iter = exon_list.iterrows() if strand == "+" else exon_list[::-1].iterrows()
    for _, exon in exon_iter:
        exon_len = exon["Exon_End"] - exon["Exon_Start"] + 1
        if genomic_start is None and tx_pos + exon_len - 1 >= start:
            offset = start - tx_pos
            genomic_start = (
                exon["Exon_Start"] + offset if strand == "+"
                else exon["Exon_End"] - offset
            )
        if genomic_stop is None and tx_pos + exon_len - 1 >= stop:
            offset = stop - tx_pos
            genomic_stop = (
                exon["Exon_Start"] + offset if strand == "+"
                else exon["Exon_End"] - offset
            )
        tx_pos += exon_len
        if genomic_start is not None and genomic_stop is not None:
            break
    return genomic_start, genomic_stop

def classify_nmd(orf_df, exon_file, nmd_out):
    exon_df = pd.read_csv(exon_file, sep="\t", dtype={"Chr": str})
    exon_df["Exon_Start"] = exon_df["Exon_Start"].astype(int)
    exon_df["Exon_End"] = exon_df["Exon_End"].astype(int)
    grouped = exon_df.groupby("Transcript_ID")

    orf_ids = set(orf_df["Transcript_ID"])
    exon_ids = set(grouped.groups.keys())
    unmatched_ids = orf_ids - exon_ids
    print(f"[!] {len(unmatched_ids)} transcripts in ORF table not found in exon structure.")
    if unmatched_ids:
        print(f"[i] Example unmatched IDs: {list(unmatched_ids)[:5]}")

    results = []
    for _, row in orf_df.iterrows():
        tx_id = row["Transcript_ID"]
        start = int(row["Start"])
        stop = int(row["Stop"])
        aa_seq = row["AA_Seq"]
        nt_seq = row["NT_Seq"]
        if tx_id not in grouped.groups:
            continue

        exon_list = grouped.get_group(tx_id).sort_values("Exon_Start").reset_index(drop=True)
        chr_name = exon_list.iloc[0]["Chr"]
        strand = exon_list.iloc[0]["Strand"]

        tx_pos = 0
        stop_exon_idx = None
        exon_lengths = (exon_list["Exon_End"] - exon_list["Exon_Start"] + 1).tolist()
        exon_ends_in_tx = [sum(exon_lengths[:i + 1]) for i in range(len(exon_lengths))]

        for i, exon in enumerate(exon_list.itertuples()):
            exon_len = exon.Exon_End - exon.Exon_Start + 1
            if tx_pos < stop <= tx_pos + exon_len:
                stop_exon_idx = i
                break
            tx_pos += exon_len

        if stop_exon_idx is None:
            continue

        is_last_exon = stop_exon_idx == (len(exon_list) - 1)
        last_junction_tx_pos = exon_ends_in_tx[-2] if len(exon_ends_in_tx) >= 2 else 0
        dist_to_last_junction = last_junction_tx_pos - stop
        dist_from_start = stop - start
        stop_exon_len = (
            exon_list.loc[stop_exon_idx]["Exon_End"]
            - exon_list.loc[stop_exon_idx]["Exon_Start"]
            + 1
        )
        genomic_start, genomic_stop = transcript_to_genomic_coords(
            start, stop, exon_list, strand
        )

        # NMD rules
        rule1 = dist_from_start < 150
        rule2 = stop_exon_len >= 407
        rule3 = dist_to_last_junction < 55
        rule4 = is_last_exon
        nmd_status = "NMD-evasive" if any([rule1, rule2, rule3, rule4]) else "NMD-sensitive"

        results.append([
            tx_id, genomic_start, genomic_stop, start, stop,
            chr_name, strand,
            dist_from_start, stop_exon_len, dist_to_last_junction, is_last_exon,
            aa_seq, nt_seq, nmd_status,
            rule1, rule2, rule3, rule4
        ])

    out_cols = [
        "Transcript_ID", "Genomic_Start", "Genomic_Stop",
        "Transcript_Start", "Transcript_Stop",
        "Chr", "Strand",
        "Distance_From_Start_Codon", "PTC_Exon_Length",
        "Distance_To_Last_Exon_Junction", "PTC_in_Last_Exon",
        "AA_Seq", "NT_Seq", "NMD_Status",
        "Rule_PTC_<150nt", "Rule_Exon>407nt",
        "Rule_<55nt_to_Junction", "Rule_Last_Exon"
    ]
    pd.DataFrame(results, columns=out_cols).to_csv(nmd_out, index=False)
    print(f"[✓] Final NMD classification saved to {nmd_out}")

############################################
# Main
############################################
def main():
    args = parse_args()

    # Resolve GTF: prefer -g, fallback to deprecated -gg
    gtf = None
    if args.gtf:
        gtf = args.gtf
    elif args.ensembl_gtf:
        print("[!] Warning: -gg/--ensembl_gtf is deprecated. Use -g/--gtf instead.")
        gtf = args.ensembl_gtf

    if not gtf:
        raise ValueError("You must provide a GTF file via -g/--gtf.")

    # Resolve transcript FASTA: prefer -f, fallback to deprecated -fg (cdna_fasta)
    fasta = None
    built_tmp_fasta = False
    if args.fasta:
        fasta = args.fasta
    elif args.cdna_fasta:
        print("[!] Warning: -fg/--cdna_fasta is deprecated. Use -f/--fasta instead.")
        fasta = args.cdna_fasta

    cpc2 = args.cpc2
    cds_fasta = args.cds_fasta
    orf_out, exon_out, nmd_out = args.outputs

    print("[1/3] Extracting exon structure...")
    extract_exons(gtf, exon_out)

    # If no transcript FASTA given, try building it from GTF + genome FASTA
    if not fasta:
        if not args.reference_fasta:
            raise ValueError(
                "No transcript FASTA (-f) provided and no genome FASTA (-r) to build it from.\n"
                "Provide either:\n"
                "  -f/--fasta : transcript/cDNA FASTA (e.g. Ensembl cDNA), or\n"
                "  -r/--reference_fasta : genome FASTA so the script can build transcript FASTA from the GTF."
            )
        print("[1b] No transcript FASTA given; building transcript FASTA from GTF + genome...")
        tmp_fh = tempfile.NamedTemporaryFile(delete=False, suffix=".fa")
        tmp_fh.close()
        fasta = build_transcript_fasta_from_exons(exon_out, args.reference_fasta, tmp_fh.name)
        built_tmp_fasta = True

    print("[2/3] Finding ORFs...")
    orf_df = find_all_orfs(fasta, cpc2)

    if cds_fasta and not orf_df.empty:
        print("[✓] Validating ORFs using CDS FASTA...")
        orf_df = validate_with_cds(orf_df, cds_fasta)

    orf_df.to_csv(orf_out, sep="\t", index=False)
    print(f"[✓] ORF results saved to {orf_out}")

    # If we built a temporary FASTA, we can safely remove it now
    if built_tmp_fasta:
        try:
            os.remove(fasta)
        except OSError:
            pass

    print("[3/3] Classifying NMD status...")
    classify_nmd(orf_df, exon_out, nmd_out)

if __name__ == "__main__":
    main()
