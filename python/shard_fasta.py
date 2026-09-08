import os
import argparse
from Bio import FastaIO


def main():
    parser = argparse.ArgumentParser(
        description="Split a multi-FASTA file into fixed-size shards for Slurm array jobs."
    )
    parser.add_argument("--fasta", type=str, required=True,
                         help="Path to the input multi-FASTA file (e.g. one length bin).")
    parser.add_argument("--shard-size", type=int, required=True,
                         help="Number of sequences per shard.")
    parser.add_argument("--outdir", type=str, required=True,
                         help="Directory to write shard_NNNN.fasta files into.")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    shard_idx = 0
    seq_in_shard = 0
    out_fh = None
    total = 0

    def open_shard(idx):
        path = os.path.join(args.outdir, f"shard_{idx:04d}.fasta")
        return open(path, "w")

    with open(args.fasta, "r") as handle:
        for title, seq in FastaIO.SimpleFastaParser(handle):
            if out_fh is None or seq_in_shard >= args.shard_size:
                if out_fh is not None:
                    out_fh.close()
                out_fh = open_shard(shard_idx)
                shard_idx += 1
                seq_in_shard = 0

            out_fh.write(f">{title}\n")
            for i in range(0, len(seq), 60):
                out_fh.write(seq[i:i + 60] + "\n")

            seq_in_shard += 1
            total += 1

    if out_fh is not None:
        out_fh.close()

    n_shards = shard_idx
    print(f"Wrote {total} sequences into {n_shards} shard(s) of up to {args.shard_size} "
          f"sequences each in: {args.outdir}")
    print(f"Use with sbatch: --array=0-{max(n_shards - 1, 0)}")


if __name__ == "__main__":
    main()
