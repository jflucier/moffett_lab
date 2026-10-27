import os
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser


def clean_sequence(seq):
    """Cleans the raw FASTA sequence string."""
    # Removes whitespace, newlines, and normalizes to uppercase
    return "".join(seq.split()).upper()


def convert_fasta_to_tsv(fasta_path, tsv_path):
    print(f"Reading FASTA from: {fasta_path}")
    print(f"Writing TSV to:     {tsv_path}")

    count = 0
    with open(fasta_path, "r") as handle, open(tsv_path, "w", encoding="utf-8") as tsv_out:
        # Write headers for the TSV file
        tsv_out.write("id\tseq\n")

        for record in SimpleFastaParser(handle):
            header_line = record[0]
            target_seq = clean_sequence(record[1])

            # Extract identifier using your specified logic
            target_id = header_line.split('|')[0].replace("|", "_").replace("/", "_")

            # Write to the output TSV file
            tsv_out.write(f"{target_id}\t{target_seq}\n")
            count += 1

    print(f"Success! Processed {count:,} sequence records.")


def main():
    parser = argparse.ArgumentParser(description="Convert a FASTA sequence file into a 2-column TSV.")
    parser.add_argument("-f", "--fasta", required=True, help="Path to input FASTA file")
    parser.add_argument("-t", "--tsv", required=True, help="Path to output destination TSV file")

    args = parser.parse_args()

    if not os.path.exists(args.fasta):
        print(f"Error: Input file '{args.fasta}' not found.")
        return

    convert_fasta_to_tsv(args.fasta, args.tsv)


if __name__ == "__main__":
    main()
