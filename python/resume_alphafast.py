#!/usr/bin/env python3
import os
import json
import shutil
import argparse


def find_incomplete_jobs(input_dir, output_dir, resume_dir):
    if not os.path.exists(resume_dir):
        os.makedirs(resume_dir)

    # 1. Track which job subdirectories contain a .cif file at their root level
    completed_jobs = set()
    if os.path.exists(output_dir):
        for entry in os.scandir(output_dir):
            if entry.is_dir():
                # Scan ONLY the files directly inside outputs/JOB_NAME/ (no recursion)
                has_cif = False
                try:
                    for f in os.scandir(entry.path):
                        if f.is_file() and f.name.endswith(".cif"):
                            has_cif = True
                            break  # Found a top-level cif, move to next directory
                except PermissionError:
                    continue

                if has_cif:
                    completed_jobs.add(entry.name)

    print(f"Detected {len(completed_jobs)} completed job directories containing structural data at the root level.")

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

            # If the job folder doesn't have a top-level .cif, it's marked to resume
            if job_name not in completed_jobs:
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
        description="Filter AlphaFast JSON inputs by checking top-level job folders for .cif outputs."
    )
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Path to the directory containing your original input .json files"
    )
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="Path to the AlphaFast output directory containing job subfolders"
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
