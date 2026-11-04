#!/usr/bin/env python3
import os
import json
import shutil
import argparse


def find_incomplete_jobs(input_dir, output_dir, resume_dir):
    if not os.path.exists(resume_dir):
        os.makedirs(resume_dir)

    # 1. Collect all completed PDB file prefixes
    completed_jobs = set()
    if os.path.exists(output_dir):
        for fname in os.listdir(output_dir):
            if fname.endswith(".pdb"):
                base_name = os.path.splitext(fname)[0]
                completed_jobs.add(base_name)
                # Strip out trailing seed numbers (e.g. "jobname_1" -> "jobname") if present
                if "_" in base_name:
                    parts = base_name.rsplit("_", 1)
                    if parts[1].isdigit():
                        completed_jobs.add(parts[0])

    print(f"Detected {len(completed_jobs)} completed structures/seeds in output folder.")

    # 2. Match against input JSONs
    failed_count = 0
    total_count = 0

    for json_file in os.listdir(input_dir):
        if not json_file.endswith(".json"):
            continue
        total_count += 1
        full_json_path = os.path.join(input_dir, json_file)

        try:
            with open(full_json_path, 'r') as f:
                data = json.load(f)

            job_name = data.get("name", "").strip()
            seeds = data.get("modelSeeds", [1])

            # Check if all targeted seeds actually generated their specific outputs
            job_complete = True
            for seed in seeds:
                specific_seed_output = f"{job_name}_{seed}"
                if specific_seed_output not in completed_jobs and job_name not in completed_jobs:
                    job_complete = False
                    break

            if not job_complete:
                shutil.copy(full_json_path, os.path.join(resume_dir, json_file))
                failed_count += 1

        except (json.JSONDecodeError, KeyError) as e:
            print(f"Warning: Could not parse configuration in {json_file}: {e}")
            shutil.copy(full_json_path, os.path.join(resume_dir, json_file))
            failed_count += 1

    print(f"Done! Evaluated {total_count} jobs. Found {failed_count} that are missing or failed.")
    print(f"Ready-to-run JSON files are stored in: '{resume_dir}'")


def main():
    parser = argparse.ArgumentParser(
        description="Filter AlphaFast JSON inputs by checking for missing PDB outputs to build a resume subset."
    )
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Path to the directory containing your original input .json files"
    )
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="Path to the AlphaFast output directory containing .pdb files"
    )
    parser.add_argument(
        "-r", "--resume_dir",
        required=True,
        help="Path to the target directory where incomplete/failed .json files will be copied"
    )

    args = parser.parse_args()
    find_incomplete_jobs(args.input_dir, args.output_dir, args.resume_dir)


if __name__ == "__main__":
    main()
