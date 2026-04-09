# Using the csv file produced by `find_intermediate_files.py`
# copy the intermediate files out.
# The csv file must be in the same directory as this script, and the files will be copied to the same directory as well.

import os
import re
import shutil
from pathlib import Path

import pandas as pd

base_dir = str(Path(__file__).resolve().parent)
print(f"Downloading to {base_dir}")

s3_bucket_files = pd.read_csv(f"{base_dir}/task_s3_bucket_map.csv")


def process_row(row):
    print(f"Processing {row.method}, dataset {row.dataset}, metric {row.metric}")

    data_dir = row.s3_path

    if pd.isna(row.metric):
        for itemname in os.listdir(data_dir):
            item_path = os.path.join(data_dir, itemname)

            # Skip if it's a directory or not .h5ad file
            if os.path.isdir(item_path) or not itemname.endswith(".h5ad"):
                continue

            if (
                "no_integration" in row.process_name
                or "perfect_integration" in row.process_name
            ):
                split_name = re.search(r"split\d+", itemname)
                if not split_name:
                    exit(f"Error: cannot find split part in {itemname}")
                split_name = split_name.group()  # Extract the matched string
                new_filename = f"{base_dir}/{row.dataset}/method_out/{row.method}_{split_name}.h5ad"

                # create directory if not exists
                os.makedirs(os.path.dirname(new_filename), exist_ok=True)

                print(f"Copying {itemname} to {os.path.basename(new_filename)}")
                shutil.copy2(item_path, new_filename)

            else:
                split_name = "2" if "process1" in row.process_name else "1"
                new_filename = f"{base_dir}/{row.dataset}/method_out/{row.method}_split{split_name}.h5ad"

                # create directory if not exists
                os.makedirs(os.path.dirname(new_filename), exist_ok=True)

                print(f"Copying {itemname} to {os.path.basename(new_filename)}")
                shutil.copy2(item_path, new_filename)

    else:
        new_filename = (
            f"{base_dir}/{row.dataset}/metric_out/{row.metric}_{row.method}.h5ad"
        )

        # create directory if not exists
        os.makedirs(os.path.dirname(new_filename), exist_ok=True)

        for itemname in os.listdir(data_dir):
            item_path = os.path.join(data_dir, itemname)

            # Skip if it's a directory or not .h5ad file
            if os.path.isdir(item_path) or not itemname.endswith(".h5ad"):
                continue

            print(f"Copying {itemname} to {os.path.basename(new_filename)}")
            shutil.copy2(item_path, new_filename)


for row in s3_bucket_files.itertuples(index=True):
    process_row(row)
