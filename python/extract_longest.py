import os
import argparse
from Bio.SeqIO import FastaIO


def main():
    parser = argparse.ArgumentParser(
        description="Extract the N longest sequences from a FASTA file for calibration runs "
                     "(worst-case timing/memory per bin before choosing shard_size)."
    )
    parser.add_argument("--fasta", type=str, required=True,
                         help="Path to the input multi-FASTA file (e.g. one length bin).")
    parser.add_argument("--n", type=int, default=30,
                         help="Number of longest sequences to extract (default: 30).")
    parser.add_argument("--out", type=str, required=True,
                         help="Output FASTA path for the calibration subset.")

    args = parser.parse_args()

    records = []
    with open(args.fasta, "r") as handle:
        for title, seq in FastaIO.SimpleFastaParser(handle):
            records.append((len(seq), title, seq))

    records.sort(key=lambda r: r[0], reverse=True)
    top = records[:args.n]

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.out, "w") as out_fh:
        for length, title, seq in top:
            out_fh.write(f">{title}\n")
            for i in range(0, len(seq), 60):
                out_fh.write(seq[i:i + 60] + "\n")

    if top:
        print(f"Extracted {len(top)} sequences from {len(records)} total.")
        print(f"Length range in subset: {top[-1][0]}-{top[0][0]} aa")
    else:
        print("No sequences found in input file.")
    print(f"Written to: {args.out}")


if __name__ == "__main__":
    main()
