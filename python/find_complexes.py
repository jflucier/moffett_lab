import os
import glob
import argparse
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio.PDB import PDBParser

# Optional: try importing tqdm for an interactive progress bar
try:
    from tqdm import tqdm

    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False


def process_single_pdb(pdb_path, distance_threshold, min_contacts):
    """Worker function executed across multiple CPU processes."""
    filename = os.path.basename(pdb_path)
    try:
        # Create a fresh parser instance per process to avoid race conditions
        pdb_parser = PDBParser(QUIET=True)
        structure = pdb_parser.get_structure("complex", pdb_path)
        model = structure[0]

        chains = list(model.get_chains())
        if len(chains) < 2:
            return None

        chain_A = chains[0]  # Bait
        chain_B = chains[1]  # Target

        # Extract CA coordinates
        coords_A = [atom.get_coord() for atom in chain_A.get_atoms() if atom.get_name() == "CA"]
        coords_B = [atom.get_coord() for atom in chain_B.get_atoms() if atom.get_name() == "CA"]

        if not coords_A or not coords_B:
            return None

        a_arr = np.array(coords_A)
        b_arr = np.array(coords_B)

        # Matrix broadcasting to compute all pairwise physical distances
        distances = np.linalg.norm(a_arr[:, np.newaxis, :] - b_arr[np.newaxis, :, :], axis=2)
        contact_count = np.sum(distances < distance_threshold)

        if contact_count >= min_contacts:
            return (filename, int(contact_count))

    except Exception as e:
        # Log to structural error log if required, or skip broken/interrupted PDB files safely
        return ("ERROR", f"Error parsing {filename}: {str(e)}")

    return None


def main():
    parser = argparse.ArgumentParser(description="High-throughput parallelized screening for ESMFold multimer PDBs.")

    parser.add_argument("-i", "--input_dir",
                        default="arabidopsis_openprot_bins_500aa_v2/bin_1001-1500aa_split_chunks/esmfold_results",
                        help="Path to the directory containing ESMFold .pdb files")
    parser.add_argument("-o", "--output",
                        default="successful_complexes.txt",
                        help="Path to save the output tab-separated results file")
    parser.add_argument("-d", "--distance",
                        type=float,
                        default=5.0,
                        help="Distance threshold in Angstroms to consider an atom-atom contact (default: 5.0)")
    parser.add_argument("-c", "--min_contacts",
                        type=int,
                        default=10,
                        help="Minimum unique CA-CA contacts required to classify as a complex (default: 10)")
    parser.add_argument("-w", "--workers",
                        type=int,
                        default=24,
                        help="Number of parallel CPU processes to spawn. (Default: Auto-detects all available cores)")

    args = parser.parse_args()

    # Locate targets
    pdb_files = glob.glob(os.path.join(args.input_dir, "*.pdb"))
    total_files = len(pdb_files)

    if total_files == 0:
        print(f"Error: No .pdb files found in {args.input_dir}")
        return

    print(f"Found {total_files:,} PDB files to screen.")
    print(f"Spawning worker pool... (Workers: {args.workers if args.workers else 'All available cores'})")

    hits = []
    errors = []

    # Execute structural parsing concurrently across individual system processes
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        # Map each future object to track processing chunks
        future_to_pdb = {
            executor.submit(process_single_pdb, path, args.distance, args.min_contacts): path
            for path in pdb_files
        }

        # Setup real-time pipeline monitoring
        if HAS_TQDM:
            progress_bar = tqdm(as_completed(future_to_pdb), total=total_files, desc="Screening Proteome")
        else:
            progress_bar = as_completed(future_to_pdb)
            print("Processing structures... (Install 'tqdm' for a live progress metric bar)")

        count = 0
        for future in progress_bar:
            count += 1
            if not HAS_TQDM and count % 5000 == 0:
                print(f"Progress: Processed {count:,} / {total_files:,} files...")

            result = future.result()
            if result:
                if result[0] == "ERROR":
                    errors.append(result[1])
                else:
                    hits.append(result)

    # Sort structural data matches by highest structural density (highest unique contacts)
    hits.sort(key=lambda x: x[1], reverse=True)

    # Flush positive complexes directly to structural log
    with open(args.output, "w") as f:
        f.write("PDB_File\tCA_Contacts\n")
        for filename, contacts in hits:
            f.write(f"{filename}\t{contacts}\n")

    # Output analytical wrap-up matrix
    print(f"\n Screening Complete!")
    print(f"Total Structures Processed: {total_files:,}")
    print(f"Successful Hits Identified : {len(hits):,}")
    print(f"Parsing Issues/Fails      : {len(errors):,}")
    print(f"Results successfully written to: {args.output}")

    if errors:
        error_log = "screening_errors.log"
        with open(error_log, "w") as ef:
            ef.write("\n".join(errors))
        print(f"Warning: Encounters {len(errors)} anomalies. Details dumped inside: {error_log}")


if __name__ == "__main__":
    main()
