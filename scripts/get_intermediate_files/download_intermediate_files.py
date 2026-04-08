# Use me to download intermediate files from AWS S3 bucket.

import argparse
import os
import re
import shutil

import pandas as pd

parser = argparse.ArgumentParser(description="Download intermediate files from S3")
parser.add_argument(
    "--base_dir", type=str, required=True, help="Base directory to download files to"
)
args = parser.parse_args()


s3_bucket_files = pd.read_csv(f"{args.base_dir}/task_s3_bucket_map.csv")

for row in s3_bucket_files.itertuples(index=True):
    # skip files that are already downloaded
    file_exists = []
    for split in ["1", "2"]:
        new_filename = (
            f"{args.base_dir}/{row.dataset}/method_out/{row.method}_{split}.h5ad"
        )
        file_exists.append(os.path.exists(new_filename))

    if all(file_exists):
        print(f"Files for {row.method} already exist, skipping...")
        continue

    print(f"Processing row: {row}")
    if pd.isna(row.metric):
        outdir = f"{args.base_dir}/{row.dataset}/method_out/temp"
        os.makedirs(outdir, exist_ok=True)
        os.system(
            f"common/scripts/fetch_task_run --input {row.s3_path} --output {outdir}/"
        )
        for itemname in os.listdir(outdir):
            item_path = os.path.join(outdir, itemname)

            # Skip if it's a directory or not .h5ad file
            if os.path.isdir(item_path) or not itemname.endswith(".h5ad"):
                continue

            elif (
                "no_integration" in row.process_name
                or "perfect_integration" in row.process_name
            ):
                split_name = re.search(r"split\d+", itemname)
                if not split_name:
                    exit(f"Error: cannot find split part in {itemname}")
                split_name = split_name.group()  # Extract the matched string
                new_filename = f"{args.base_dir}/{row.dataset}/method_out/{row.method}_{split_name}.h5ad"
                print(f"Renaming {item_path} to {new_filename}")
                shutil.move(item_path, new_filename)

            else:
                split_name = "2" if "process1" in row.process_name else "1"
                new_filename = f"{args.base_dir}/{row.dataset}/method_out/{row.method}_split{split_name}.h5ad"
                print(f"Renaming {item_path} to {new_filename}")
                shutil.move(item_path, new_filename)
        shutil.rmtree(outdir)

    else:
        new_filename = (
            f"{args.base_dir}/{row.dataset}/metric_out/{row.metric}_{row.method}.h5ad"
        )

        if os.path.exists(new_filename):
            print(f"File {new_filename} already exists, skipping...")
            continue

        outdir = f"{args.base_dir}/{row.dataset}/metric_out/temp"
        os.makedirs(outdir, exist_ok=True)
        os.system(
            f"common/scripts/fetch_task_run --input {row.s3_path} --output {outdir}/"
        )
        for itemname in os.listdir(outdir):
            item_path = os.path.join(outdir, itemname)

            # Skip if it's a directory or the .h5ad file
            if os.path.isdir(item_path) or not itemname.endswith(".h5ad"):
                continue

            else:
                # TODO the split is not needed.. made a mistake before.
                print(f"Renaming {item_path} to {new_filename}")
                shutil.move(item_path, new_filename)

            # clean up temp folder
        shutil.rmtree(outdir)
