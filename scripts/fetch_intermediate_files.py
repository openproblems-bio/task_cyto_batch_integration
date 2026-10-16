#!/usr/bin/env python3
"""
Parse a Nextflow log file to infer task work directories, then copy the
intermediate .h5ad files into an organised output directory.

Supports both SLURM (local paths) and AWS Batch (s3:// paths) logs.
For SLURM, files are copied directly from the local work directory.
For AWS Batch, files are downloaded via `common/scripts/fetch_task_run`.

Usage:
    python fetch_intermediate_files.py <log_file> [--output-dir <dir>] [--save-csv]
"""

import argparse
import os
import re
import shutil
from pathlib import Path

import pandas as pd


# ── Regex patterns ────────────────────────────────────────────────────────────
#
# SLURM line example:
#   [SLURM] submitted process auto:...:harmonypy_process (human_blood.harmonypy)
#   > jobId: 25169450; workDir: /vast/scratch/.../6f/488d302c4e33171f2d69c8107b7a9c
#
# AWS Batch line example:
#   [AWS BATCH] Process `auto:...:harmonypy_process (human_blood.harmonypy)`
#   submitted > job=70af3d53; work-dir=s3://openproblems-work/scratch/.../60/c803767c

AWS_PATTERN = re.compile(
    r"\[AWS BATCH\] Process `.*:([\w\d_]+ \([^)]+\))`.*?work-dir=(s3://[^ ;]+)"
)
SLURM_PATTERN = re.compile(
    r"\[SLURM\] submitted process .*:([\w\d_]+ \([^)]+\)).*?workDir: ([^ ]+)"
)


# ── Log parsing ───────────────────────────────────────────────────────────────

def detect_executor(log_file: Path) -> str:
    """Return 'aws' or 'slurm' based on which executor pattern appears first."""
    with open(log_file) as f:
        for line in f:
            if "[AWS BATCH]" in line:
                return "aws"
            if "[SLURM] submitted process" in line:
                return "slurm"
    raise ValueError(
        f"Could not detect executor type (AWS Batch or SLURM) in {log_file}"
    )


def parse_task_component(task_component: str) -> tuple:
    """
    Split 'process_name (dataset[.method[.metric]])' into its parts.

    Task component examples:
      - "extract_uns_metadata_process (human_blood_mass_cytometry)"
      - "harmonypy_process (human_blood_mass_cytometry.harmonypy)"
      - "flowsom_mapping_similarity_process (mouse_spleen.cytonorm.flowsom_mapping_similarity)"

    Returns (process_name, dataset_name, method_name, metric_name),
    where method_name and metric_name may be None.
    """
    process_name = task_component.split(" ")[0]
    task_info = task_component.split(" ")[1].strip("()").split(".")
    dataset_name = task_info[0]
    method_name = task_info[1] if len(task_info) >= 2 else None
    metric_name = task_info[2] if len(task_info) >= 3 else None
    return process_name, dataset_name, method_name, metric_name


def parse_log(log_file: Path) -> tuple:
    """
    Parse the Nextflow log and return (executor, DataFrame) where executor is
    'aws' or 'slurm' and the DataFrame holds task-to-work-dir mappings.
    """
    executor = detect_executor(log_file)
    pattern = AWS_PATTERN if executor == "aws" else SLURM_PATTERN
    print(f"Detected executor: {executor.upper()}")

    rows = []
    with open(log_file) as f:
        for line_number, line in enumerate(f, start=1):
            match = pattern.search(line)
            if match:
                process_name, dataset, method, metric = parse_task_component(match.group(1))
                work_dir = match.group(2).strip()
                rows.append((line_number, process_name, dataset, method, metric, work_dir))

    df = pd.DataFrame(
        rows,
        columns=["line", "process_name", "dataset", "method", "metric", "work_dir"],
    )
    df.sort_values(by="process_name", inplace=True, ignore_index=True)
    return executor, df


# ── Shared helpers ────────────────────────────────────────────────────────────

def resolve_split_name(process_name: str, filename: str) -> str:
    """
    Determine the split label (e.g. 'split1') for a method output file.

    For no_integration / perfect_integration, the split number is encoded in
    the output filename itself (e.g. 'output_split3.h5ad').
    For all other methods, the split is inferred from the process name:
    'process1' → split2, anything else → split1.
    """
    if "no_integration" in process_name or "perfect_integration" in process_name:
        match = re.search(r"split\d+", filename)
        if not match:
            raise ValueError(
                f"Cannot find 'splitN' in filename '{filename}' "
                f"for process '{process_name}'"
            )
        return match.group()
    else:
        return "split2" if "process1" in process_name else "split1"


# ── SLURM: copy directly from local work directory ────────────────────────────

def copy_method_output(row, output_dir: Path) -> None:
    """
    Copy method output .h5ad to:
        <output_dir>/<dataset>/method_out/<method>_<split>.h5ad
    """
    for itemname in os.listdir(row.work_dir):
        item_path = os.path.join(row.work_dir, itemname)

        if os.path.isdir(item_path) or not itemname.endswith(".h5ad"):
            continue

        split_name = resolve_split_name(row.process_name, itemname)
        dest = output_dir / row.dataset / "method_out" / f"{row.method}_{split_name}.h5ad"
        dest.parent.mkdir(parents=True, exist_ok=True)

        print(f"  {itemname}  ->  {dest.relative_to(output_dir)}")
        shutil.copy2(item_path, dest)


def copy_metric_output(row, output_dir: Path) -> None:
    """
    Copy metric output .h5ad to:
        <output_dir>/<dataset>/metric_out/<metric>_<method>.h5ad
    """
    dest = output_dir / row.dataset / "metric_out" / f"{row.metric}_{row.method}.h5ad"
    dest.parent.mkdir(parents=True, exist_ok=True)

    for itemname in os.listdir(row.work_dir):
        item_path = os.path.join(row.work_dir, itemname)

        if os.path.isdir(item_path) or not itemname.endswith(".h5ad"):
            continue

        print(f"  {itemname}  ->  {dest.relative_to(output_dir)}")
        shutil.copy2(item_path, dest)


def copy_files_slurm(df: pd.DataFrame, output_dir: Path) -> None:
    """Iterate over parsed tasks and copy each output file to the organised tree."""
    for row in df.itertuples(index=False):
        has_method = not pd.isna(row.method)
        has_metric = not pd.isna(row.metric)

        if not has_method and not has_metric:
            print(f"Skipping {row.process_name} ({row.dataset}) — no method or metric")
            continue

        label = f"{row.method}, dataset {row.dataset}"
        if has_metric:
            label += f", metric {row.metric}"
        print(f"Processing {label}")

        if has_metric:
            copy_metric_output(row, output_dir)
        else:
            copy_method_output(row, output_dir)


# ── AWS Batch: download from S3 via fetch_task_run ────────────────────────────

def download_method_output_aws(row, output_dir: Path) -> None:
    """
    Download method output from S3 to a temp dir using fetch_task_run,
    then move each .h5ad to <output_dir>/<dataset>/method_out/<method>_<split>.h5ad.

    Skips if both split1 and split2 files already exist.
    """
    # Check whether both expected split files are already present
    existing = [
        os.path.exists(output_dir / row.dataset / "method_out" / f"{row.method}_{split}.h5ad")
        for split in ["1", "2"]
    ]
    if all(existing):
        print(f"  Files for {row.method} already exist, skipping...")
        return

    temp_dir = output_dir / row.dataset / "method_out" / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    os.system(f"common/scripts/fetch_task_run --input {row.work_dir} --output {temp_dir}/")

    for itemname in os.listdir(temp_dir):
        item_path = os.path.join(temp_dir, itemname)

        if os.path.isdir(item_path) or not itemname.endswith(".h5ad"):
            continue

        split_name = resolve_split_name(row.process_name, itemname)
        dest = output_dir / row.dataset / "method_out" / f"{row.method}_{split_name}.h5ad"

        print(f"  {itemname}  ->  {dest.relative_to(output_dir)}")
        shutil.move(item_path, dest)

    shutil.rmtree(temp_dir)


def download_metric_output_aws(row, output_dir: Path) -> None:
    """
    Download metric output from S3 to a temp dir using fetch_task_run,
    then move the .h5ad to <output_dir>/<dataset>/metric_out/<metric>_<method>.h5ad.

    Skips if the destination file already exists.
    """
    dest = output_dir / row.dataset / "metric_out" / f"{row.metric}_{row.method}.h5ad"

    if dest.exists():
        print(f"  File {dest} already exists, skipping...")
        return

    temp_dir = output_dir / row.dataset / "metric_out" / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    os.system(f"common/scripts/fetch_task_run --input {row.work_dir} --output {temp_dir}/")

    for itemname in os.listdir(temp_dir):
        item_path = os.path.join(temp_dir, itemname)

        if os.path.isdir(item_path) or not itemname.endswith(".h5ad"):
            continue

        print(f"  {itemname}  ->  {dest.relative_to(output_dir)}")
        shutil.move(item_path, dest)

    shutil.rmtree(temp_dir)


def download_files_aws(df: pd.DataFrame, output_dir: Path) -> None:
    """Iterate over parsed tasks and download each output file from S3."""
    for row in df.itertuples(index=False):
        has_method = not pd.isna(row.method)
        has_metric = not pd.isna(row.metric)

        if not has_method and not has_metric:
            print(f"Skipping {row.process_name} ({row.dataset}) — no method or metric")
            continue

        label = f"{row.method}, dataset {row.dataset}"
        if has_metric:
            label += f", metric {row.metric}"
        print(f"Processing {label}")

        if has_metric:
            download_metric_output_aws(row, output_dir)
        else:
            download_method_output_aws(row, output_dir)


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parse a Nextflow log file, infer work directories for each task, "
            "and copy intermediate .h5ad files into an organised directory."
        )
    )
    parser.add_argument(
        "log_file",
        type=Path,
        help="Path to the Nextflow .log file (e.g. nf-14vBbXiMZU8YjS.log)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Directory to write organised files into "
            "(default: 'intermediate_files/' next to the log file)"
        ),
    )
    parser.add_argument(
        "--save-csv",
        action="store_true",
        help="Save the task-to-work-dir mapping as task_work_dir_map.csv in the output directory",
    )
    args = parser.parse_args()

    log_file = args.log_file.resolve()
    if not log_file.exists():
        parser.error(f"Log file not found: {log_file}")

    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else log_file.parent / "intermediate_files"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Log file  : {log_file}")
    print(f"Output dir: {output_dir}")

    executor, df = parse_log(log_file)
    print(f"Found {len(df)} task entries\n")

    if args.save_csv:
        csv_path = output_dir / "task_work_dir_map.csv"
        df.to_csv(csv_path, index=False)
        print(f"Saved mapping CSV to {csv_path}\n")

    if executor == "aws":
        download_files_aws(df, output_dir)
    else:
        copy_files_slurm(df, output_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
