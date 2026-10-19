import os
import argparse
from Bio import SeqIO


def main():
    # Set up command-line arguments
    parser = argparse.ArgumentParser(description="Find failed ESMFold jobs and output them to a new FASTA file.")
    parser.add_argument("-i", "--fasta_in", required=True, help="Path to the original input FASTA file.")
    parser.add_argument("-d", "--pdb_dir", required=True, help="Directory containing the successful *.pdb files.")
    parser.add_argument("-o", "--fasta_out", required=True,
                        help="Path to save the generated failed sequences FASTA file.")
    parser.add_argument("-b", "--bin_size", type=int, default=0,
                        help="Maximum number of sequences per output file. If 0 or omitted, all go into one file.")

    args = parser.parse_args()

    # 1. Read existing PDB filenames into a set for fast lookup
    if not os.path.isdir(args.pdb_dir):
        print(f"Error: Directory '{args.pdb_dir}' does not exist.")
        return

    successful_folds = {
        os.path.basename(f)[:-4] for f in os.listdir(args.pdb_dir) if f.endswith('.pdb')
    }
    print(f"Found {len(successful_folds)} successful PDB files.")

    # 2. Parse the input FASTA and filter for missing entries
    if not os.path.isfile(args.fasta_in):
        print(f"Error: FASTA file '{args.fasta_in}' does not exist.")
        return

    failed_records = []
    for record in SeqIO.parse(args.fasta_in, "fasta"):
        # Match how ESMFold names files based on the first chunk of the header before a pipe
        header_id = record.id.split('|')[0].replace("|", "_").replace("/", "_")

        if header_id not in successful_folds:
            failed_records.append(record)

    # 3. Write out the failed sequences
    total_failed = len(failed_records)
    if total_failed == 0:
        print("All sequences were successfully folded! No failed entries found.")
        return

    # If bin_size is set to a valid positive integer, split into chunks
    if args.bin_size > 0:
        base, ext = os.path.splitext(args.fasta_out)
        chunk_idx = 1

        for i in range(0, total_failed, args.bin_size):
            chunk = failed_records[i:i + args.bin_size]
            chunk_filename = f"{base}_chunk_{chunk_idx}{ext}"

            SeqIO.write(chunk, chunk_filename, "fasta")
            print(f"Generated '{chunk_filename}' with {len(chunk)} sequences.")
            chunk_idx += 1

        print(f"Success! Total of {total_failed} failed sequences split across {chunk_idx - 1} files.")
    else:
        # Standard behavior: write everything into a single file
        SeqIO.write(failed_records, args.fasta_out, "fasta")
        print(f"Success! Generated '{args.fasta_out}' with {total_failed} failed sequences.")


if __name__ == "__main__":
    main()
