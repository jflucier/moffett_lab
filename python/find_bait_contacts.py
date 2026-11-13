#!/usr/bin/env python3
import os
import json
import csv
import argparse
import numpy as np
from Bio.PDB import MMCIFParser


def count_interchain_contacts(cif_path, cutoff=5.0):
    """
    Parses a .cif file and returns the number of heavy-atom pairs within `cutoff` Angstroms
    between different chains, plus the minimum inter-chain distance.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", cif_path)
    model = structure[0]
    chains = list(model.get_chains())

    if len(chains) < 2:
        return 0, np.nan

    chain_coords = []
    for chain in chains:
        coords = [
            atom.get_coord()
            for res in chain
            for atom in res
            if atom.element != "H"
        ]
        if coords:
            chain_coords.append(np.array(coords))

    if len(chain_coords) < 2:
        return 0, np.nan

    total_contacts = 0
    min_dist = np.inf

    for i in range(len(chain_coords)):
        for j in range(i + 1, len(chain_coords)):
            c1 = chain_coords[i]
            c2 = chain_coords[j]

            for start in range(0, len(c1), 500):
                sub_c1 = c1[start:start + 500]
                dists = np.linalg.norm(sub_c1[:, None, :] - c2[None, :, :], axis=-1)

                current_min = float(np.min(dists))
                if current_min < min_dist:
                    min_dist = current_min

                total_contacts += int(np.sum(dists <= cutoff))

    return total_contacts, round(min_dist, 2) if min_dist != np.inf else np.nan


def analyze_contacts(input_dir, output_tsv, dist_cutoff, min_iptm, min_contacts, max_pae, skip_clash):
    results = []
    folders = sorted(os.listdir(input_dir))
    print(f"Scanning {len(folders)} folders in '{input_dir}'...")
    print(
        f"Thresholds -> Distance <= {dist_cutoff} Å | ipTM >= {min_iptm} | Contacts >= {min_contacts} | Max PAE <= {max_pae}")

    for idx, folder in enumerate(folders, 1):
        folder_path = os.path.join(input_dir, folder)
        if not os.path.isdir(folder_path):
            continue

        cif_path = os.path.join(folder_path, f"{folder}_model.cif")
        summary_path = os.path.join(folder_path, f"{folder}_summary_confidences.json")

        if not os.path.exists(cif_path) or not os.path.exists(summary_path):
            continue

        # 1. Read AlphaFold 3 confidence metrics
        with open(summary_path, "r") as f:
            conf = json.load(f)

        iptm = conf.get("iptm", 0.0) or 0.0
        ptm = conf.get("ptm", 0.0) or 0.0
        ranking_score = conf.get("ranking_score", 0.0) or 0.0
        has_clash = conf.get("has_clash", 0.0)

        if skip_clash and has_clash > 0:
            continue

        if iptm < min_iptm:
            continue

        # Extract minimum off-diagonal PAE across chain pairs
        pae_matrix = conf.get("chain_pair_pae_min", [])
        min_inter_pae = np.nan
        if len(pae_matrix) >= 2:
            off_diag = [
                pae_matrix[i][j]
                for i in range(len(pae_matrix))
                for j in range(len(pae_matrix))
                if i != j and pae_matrix[i][j] is not None
            ]
            if off_diag:
                min_inter_pae = round(min(off_diag), 2)

        if max_pae is not None and not np.isnan(min_inter_pae) and min_inter_pae > max_pae:
            continue

        # 2. Measure 3D physical contacts
        contacts, min_dist = count_interchain_contacts(cif_path, dist_cutoff)

        if contacts < min_contacts:
            continue

        results.append({
            "job_name": folder,
            "iptm": round(iptm, 4),
            "ptm": round(ptm, 4),
            "ranking_score": round(ranking_score, 4),
            "min_inter_pae": min_inter_pae,
            f"atomic_contacts_{dist_cutoff}A": contacts,
            "min_distance_A": min_dist,
            "has_clash": has_clash
        })

        if idx % 200 == 0:
            print(f"Processed {idx}/{len(folders)} folders...")

    # Sort by ipTM descending, then by contact count descending
    contact_col = f"atomic_contacts_{dist_cutoff}A"
    results.sort(key=lambda x: (x["iptm"], x[contact_col]), reverse=True)

    # Write results to TSV
    fieldnames = [
        "job_name", "iptm", "ptm", "ranking_score",
        "min_inter_pae", contact_col, "min_distance_A", "has_clash"
    ]
    with open(output_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print(f"\nDone! Saved {len(results)} matching structures to '{output_tsv}'.")


def main():
    parser = argparse.ArgumentParser(
        description="Screen AlphaFold 3 / AlphaFast completed structures for bait-protein interface contacts."
    )
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Path to the completed outputs directory (e.g., outputs/completed)"
    )
    parser.add_argument(
        "-o", "--output_tsv",
        default="bait_contacts_summary.tsv",
        help="Output TSV file path (default: bait_contacts_summary.tsv)"
    )
    parser.add_argument(
        "-d", "--dist_cutoff",
        type=float,
        default=5.0,
        help="Heavy-atom distance cutoff in Angstroms to define a contact (default: 5.0)"
    )
    parser.add_argument(
        "--min_iptm",
        type=float,
        default=0.0,
        help="Minimum ipTM score required to include a hit (default: 0.0, exports all)"
    )
    parser.add_argument(
        "--min_contacts",
        type=int,
        default=0,
        help="Minimum number of inter-chain atomic contacts required (default: 0)"
    )
    parser.add_argument(
        "--max_pae",
        type=float,
        default=None,
        help="Maximum allowed minimum inter-chain PAE in Angstroms, e.g., 5.0 (default: None)"
    )
    parser.add_argument(
        "--skip_clash",
        action="store_true",
        help="Exclude structures flagged with steric clashes (has_clash > 0)"
    )

    args = parser.parse_args()
    analyze_contacts(
        input_dir=args.input_dir,
        output_tsv=args.output_tsv,
        dist_cutoff=args.dist_cutoff,
        min_iptm=args.min_iptm,
        min_contacts=args.min_contacts,
        max_pae=args.max_pae,
        skip_clash=args.skip_clash
    )


if __name__ == "__main__":
    main()
