#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

############################################
# Argument Parser
############################################
def parse_args():
    parser = argparse.ArgumentParser(description="Unified NMD classification pipeline.")
    parser.add_argument("-f", "--fasta", required=False, help="Transcript FASTA file (used unless -fg is provided)")
    parser.add_argument("-g", "--gtf", required=False, help="Transcript GTF annotation file (used unless -gg is provided)")
    parser.add_argument("-c", "--cpc2", required=False, help="Optional CPC2 output file")
    parser.add_argument("-fg", "--cdna_fasta", required=False, help="Ensembl cDNA FASTA file")
    parser.add_argument("-gg", "--ensembl_gtf", required=False, help="Ensembl GTF file")
    parser.add_argument("-cg", "--cds_fasta", required=False, help="Ensembl CDS FASTA file for validation")
    parser.add_argument("-o", "--outputs", nargs=3, required=True, metavar=('orf_out', 'exon_out', 'nmd_out'),
                        help="Output files: ORF.tsv EXON.tsv NMD.csv")
    return parser.parse_args()

############################################
# Step 1: Extract Exons
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

def extract_exons(gtf_file, exon_out):
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

def find_all_orfs(fasta_path, orf_out=None, cpc2_path=None):
    all_orfs = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        tx_id = record.id.split(".")[0]
        seq = record.seq
        orfs = find_orfs_in_sequence(seq, "+", tx_id)
        if not orfs:
            continue
        if cpc2_path:
            all_orfs.extend(orfs)
        else:
            longest = max(orfs, key=lambda x: x["Length"])
            all_orfs.append(longest)

    df = pd.DataFrame(all_orfs)

    if cpc2_path:
        cpc2_df = pd.read_csv(cpc2_path, sep="\t")
        cpc2_df["putative_peptide"] = cpc2_df["putative_peptide"].str.replace("*", "", regex=False)
        cpc2_df = cpc2_df.rename(columns={"#ID": "Transcript_ID"})
        df = df.merge(cpc2_df[["Transcript_ID", "putative_peptide"]], on="Transcript_ID", how="left")
        df["Match_CPC2"] = df.apply(
            lambda r: r["AA_Seq"] == r["putative_peptide"] if pd.notna(r["putative_peptide"]) else False, axis=1
        )
        df = df[df["Match_CPC2"]].drop(columns=["putative_peptide", "Match_CPC2"])

    return df

def validate_with_cds(orf_df, cds_fasta):
    cds_dict = {record.id.split(".")[0]: str(record.seq) for record in SeqIO.parse(cds_fasta, "fasta")}
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
            genomic_start = exon["Exon_Start"] + offset if strand == "+" else exon["Exon_End"] - offset
        if genomic_stop is None and tx_pos + exon_len - 1 >= stop:
            offset = stop - tx_pos
            genomic_stop = exon["Exon_Start"] + offset if strand == "+" else exon["Exon_End"] - offset
        tx_pos += exon_len
        if genomic_start is not None and genomic_stop is not None:
            break
    return genomic_start, genomic_stop

def classify_nmd(orf_df, exon_file, fasta_path, nmd_out):
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
        stop_exon_len = exon_list.loc[stop_exon_idx]["Exon_End"] - exon_list.loc[stop_exon_idx]["Exon_Start"] + 1
        genomic_start, genomic_stop = transcript_to_genomic_coords(start, stop, exon_list, strand)

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
        "Rule_PTC_<150nt", "Rule_Exon>407nt", "Rule_<55nt_to_Junction", "Rule_Last_Exon"
    ]
    pd.DataFrame(results, columns=out_cols).to_csv(nmd_out, index=False)
    print(f"[✓] Final NMD classification saved to {nmd_out}")

############################################
# Main
############################################
def main():
    args = parse_args()

    fasta = args.fasta if args.fasta else args.cdna_fasta
    gtf = args.gtf if args.gtf else args.ensembl_gtf
    cpc2 = args.cpc2
    cds_fasta = args.cds_fasta
    orf_out, exon_out, nmd_out = args.outputs

    print("[1/3] Extracting exon structure...")
    extract_exons(gtf, exon_out)

    print("[2/3] Finding ORFs...")
    orf_df = find_all_orfs(fasta, None, cpc2)

    if cds_fasta:
        print("[✓] Validating ORFs using CDS.fa...")
        orf_df = validate_with_cds(orf_df, cds_fasta)

    orf_df.to_csv(orf_out, sep="\t", index=False)
    print(f"[✓] ORF results saved to {orf_out}")

    print("[3/3] Classifying NMD status...")
    classify_nmd(orf_df, exon_out, fasta, nmd_out)

if __name__ == "__main__":
    main()
