#!/usr/bin/env python3
"""
caps_pipeline.py

新增 --filter 参数：
  - 当 filter=True（默认）时，执行 QUAL>=30 & DP>=10 & biallelic SNP 过滤
  - 当 filter=False 时，跳过以上过滤，直接使用所有 VCF 记录
"""
import argparse, gzip
from cyvcf2 import VCF
from Bio import SeqIO


def open_maybe_gzip(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def load_enzymes(enzyme_file):
    enz = {}
    with open(enzyme_file) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue

            name, site, pattern = ln.strip().split()
            recog      = pattern.replace("'", "").upper()
            cut_offset = pattern.index("'")
            enz[name]  = (recog, cut_offset)
    return enz


def extract_flank(ref_dict, chrom, pos, flank):
    seq   = ref_dict[chrom].seq
    start = max(0, pos - flank - 1)
    mid   = pos - 1
    end   = min(len(seq), pos + flank)
    return seq[start:mid], seq[mid:end]


def detect_caps(ref_seq, alt_seq, enzymes, flank_len):
    results   = []
    snp_index = flank_len
    for name, (site, offset) in enzymes.items():
        L = len(site)
        ref_cuts = {
            i+offset for i in range(len(ref_seq)-L+1)
            if ref_seq[i:i+L] == site
               and i <= snp_index < i+L
        }

        alt_cuts = {
            i+offset for i in range(len(alt_seq)-L+1)
            if alt_seq[i:i+L] == site
               and i <= snp_index < i+L
        }

        if ref_cuts and not alt_cuts:
            results.append((name, 'loss'))

        elif alt_cuts and not ref_cuts:
            results.append((name, 'gain'))
    return results


def mark_snp_lowercase(seq, snp_index):
    """
    将 seq 中索引 snp_index 处字符设为小写，其他保持大写。
    """
    ch = seq[snp_index].lower()      # 将 SNP 位置字符小写 :contentReference[oaicite:1]{index=1}
    return seq[:snp_index] + ch + seq[snp_index+1:]


def main():
    parser = argparse.ArgumentParser(description="Design Caps Marker through standard vcf file.")
    parser.add_argument("-v","--vcf",   required=True,
                        help="Input SNP VCF (.vcf/.vcf.gz)")
    parser.add_argument("-r","--ref",   required=True,
                        help="Reference FASTA")
    parser.add_argument("-e","--enz",   required=True,
                        help="3 colnums Enzyme file: Enzyme_name<TAB>RecognitionSeq<TAB>Pattern")
    parser.add_argument("-f","--flank", type=int, default=300,
                        help="flank sequences, default is 300")
    parser.add_argument("-o","--out",   required=True,
                        help="Output TSV")
    parser.add_argument("--filter", action="store_true", default=True,
                        help="Whether to perform filtering on VCF(QUAL>=30 & DP>=10 & biallelic SNP), default is on")  # :contentReference[oaicite:0]{index=0}
    parser.add_argument("--no-filter", dest="filter", action="store_false",
                        help="skip VCF filter")
    args = parser.parse_args()

    # loading reference 
    ref_dict = SeqIO.to_dict(SeqIO.parse(open_maybe_gzip(args.ref), "fasta"))
    
    # loading enzyme file
    enzymes  = load_enzymes(args.enz)
    
    # loading vcf file
    vcf = VCF(args.vcf)
    if args.filter:
        vcf.set_filter("PASS")

    out = open(args.out, "w")
    out.write("MarkerID\tChrom\tPos\tREF\tALT\tEnzyme\tChange\tFlankSeq\n")

    idx = 1
    for rec in vcf:
        # vcf filtering
        if args.filter:
            if not rec.is_snp or len(rec.ALT) != 1:
                continue

            if rec.QUAL < 30 or rec.INFO.get("DP",0) < 10:
                continue

        chrom, pos = rec.CHROM, rec.POS
        ref, alt   = rec.REF, rec.ALT[0]

        # flank sequence
        up_seq, down_seq = extract_flank(ref_dict, chrom, pos, args.flank)
        ref_seq = str(up_seq) + ref + str(down_seq)
        alt_seq = str(up_seq) + alt + str(down_seq)

        # Detect CAPS gain/loss across SNPs
        for enz, change in detect_caps(ref_seq, alt_seq, enzymes, len(up_seq)):
            marker = f"caps_{idx:05d}"
            seq_full = ref_seq if change == 'loss' else alt_seq
            seq_marked = mark_snp_lowercase(seq_full, len(up_seq))
            out.write(f"{marker}\t{chrom}\t{pos}\t{ref}\t{alt}\t"
                      f"{enz}\t{change}\t{seq_marked}\n")
            idx += 1

    out.close()
    print(f"[+] Total {idx-1} CAPS Markers.")


if __name__ == "__main__":
    main()

