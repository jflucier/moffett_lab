import os
import glob
import argparse
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio.PDB import PDBParser

try:
    from tqdm import tqdm

    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False


def process_single_pdb(pdb_path, distance_threshold, min_density_pct):
    """Worker function calculating Interface Density per PDB file."""
    filename = os.path.basename(pdb_path)
    try:
        pdb_parser = PDBParser(QUIET=True)
        structure = pdb_parser.get_structure("complex", pdb_path)
        model = structure

        chains = list(model.get_chains())
        if len(chains) < 2:
            return None

        chain_A = chains[0]  # Bait protein
        chain_B = chains[1]  # Target protein from Arabidopsis proteome

        # Extract CA coordinates
        coords_A = [atom.get_coord() for atom in chain_A.get_atoms() if atom.get_name() == "CA"]
        coords_B = [atom.get_coord() for atom in chain_B.get_atoms() if atom.get_name() == "CA"]

        if not coords_A or not coords_B:
            return None

        bait_length = len(coords_A)
        a_arr = np.array(coords_A)
        b_arr = np.array(coords_B)

        # Calculate full distance matrix
        distances = np.linalg.norm(a_arr[:, np.newaxis, :] - b_arr[np.newaxis, :, :], axis=2)

        # STRATEGY 3 LOGIC:
        # Find which specific residues in Chain A (Bait) have AT LEAST ONE contact within distance_threshold
        # distances < distance_threshold creates a boolean matrix (Bait_Residues x Target_Residues)
        # axis=1 checks across the target residues for each bait residue
        interacting_bait_residues = np.any(distances < distance_threshold, axis=1)
        contacting_bait_count = np.sum(interacting_bait_residues)

        # Calculate what percentage of the bait is covered by the target interface
        interface_density = (contacting_bait_count / bait_length) * 100.0

        if interface_density >= min_density_pct:
            return (filename, round(interface_density, 2), int(contacting_bait_count))

    except Exception as e:
        return ("ERROR", f"Error parsing {filename}: {str(e)}")

    return None


def main():
    parser = argparse.ArgumentParser(description="Strategy 3: Interface Packing Density Complex Screener.")

    parser.add_argument("-i", "--input_dir", required=True,
                        help="Path to the directory containing ESMFold .pdb files")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to save the output tab-separated results file")
    parser.add_argument("-d", "--distance", type=float, default=4.5,
                        help="Distance threshold in Angstroms for structural contact (default: 4.5Å)")
    parser.add_argument("-p", "--min_density", type=float, default=10.0,
                        help="Minimum % of Bait residues that must touch the target to pass (default: 10.0%)")
    parser.add_argument("-w", "--workers", type=int, default=None,
                        help="Number of parallel CPU processes. (Default: Auto-detects all cores)")

    args = parser.parse_args()

    pdb_files = glob.glob(os.path.join(args.input_dir, "*.pdb"))
    total_files = len(pdb_files)

    if total_files == 0:
        print(f"Error: No .pdb files found in {args.input_dir}")
        return

    print(
        f"Found {total_files:,} structures. Filtering with Interface Density >= {args.min_density}% at {args.distance}Å...")

    hits = []
    errors = []

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_pdb = {
            executor.submit(process_single_pdb, path, args.distance, args.min_density): path
            for path in pdb_files
        }

        if HAS_TQDM:
            progress_bar = tqdm(as_completed(future_to_pdb), total=total_files, desc="Density Screening")
        else:
            progress_bar = as_completed(future_to_pdb)
            print("Processing structures...")

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

    # Sort results by highest interface packing density (%)
    hits.sort(key=lambda x: x[1], reverse=True)

    with open(args.output, "w") as f:
        f.write("PDB_File\tInterface_Density_Pct\tBait_Contact_Residues\n")
        for filename, density, contacts in hits:
            f.write(f"{filename}\t{density}%\t{contacts}\n")

    print(f"\n Screening Complete!")
    print(f"Total Processed: {total_files:,} | Hits: {len(hits):,} | Errors: {len(errors):,}")
    print(f"Results saved to: {args.output}")


if __name__ == "__main__":
    main()
