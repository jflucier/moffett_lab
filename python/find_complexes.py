import os
import glob
import argparse
import numpy as np
from Bio.PDB import PDBParser


def main():
    # --- CLI Argument Parsing ---
    parser = argparse.ArgumentParser(description="Screen ESMFold multimer PDBs for true bait-prey complexes.")

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

    args = parser.parse_args()

    # --- PDB Screening Logic ---
    pdb_parser = PDBParser(QUIET=True)
    results = []

    print(f"Scanning PDB files in: {args.input_dir}")
    print(f"Criteria: Distance <= {args.distance}Å | Min Contacts >= {args.min_contacts}")

    pdb_files = glob.glob(os.path.join(args.input_dir, "*.pdb"))

    if not pdb_files:
        print(f"Error: No .pdb files found in {args.input_dir}")
        return

    for pdb_path in pdb_files:
        filename = os.path.basename(pdb_path)
        try:
            structure = pdb_parser.get_structure("complex", pdb_path)
            model = structure[0]

            chains = list(model.get_chains())
            if len(chains) < 2:
                continue

            chain_A = chains[0]  # Bait
            chain_B = chains[1]  # Target

            coords_A = [atom.get_coord() for atom in chain_A.get_atoms() if atom.get_name() == "CA"]
            coords_B = [atom.get_coord() for atom in chain_B.get_atoms() if atom.get_name() == "CA"]

            if not coords_A or not coords_B:
                continue

            a_arr = np.array(coords_A)
            b_arr = np.array(coords_B)

            distances = np.linalg.norm(a_arr[:, np.newaxis, :] - b_arr[np.newaxis, :, :], axis=2)
            contact_count = np.sum(distances < args.distance)

            if contact_count >= args.min_contacts:
                results.append((filename, contact_count))
                print(f"Found complex: {filename} ({contact_count} CA contacts)")

        except Exception as e:
            print(f"Error parsing {filename}: {e}")

    # Sort and export
    results.sort(key=lambda x: x[1], reverse=True)

    with open(args.output, "w") as f:
        f.write("PDB_File\tCA_Contacts\n")
        for filename, contacts in results:
            f.write(f"{filename}\t{contacts}\n")

    print(f"\nDone! Results saved to {args.output}")


if __name__ == "__main__":
    main()
