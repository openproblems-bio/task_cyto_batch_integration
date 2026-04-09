# Run me to find the AWS s3 paths for intermediate files for each task.

import argparse
import re

import pandas as pd

parser = argparse.ArgumentParser(
    description="Extract intermediate file paths from log file"
)
parser.add_argument(
    "--log_file_path", type=str, required=True, help="Path to the log file to parse"
)
parser.add_argument(
    "--output_csv_path",
    type=str,
    default="task_s3_bucket_map.csv",
    help="Path to save the output CSV file with task and S3 path mapping",
)
args = parser.parse_args()

# Feb-04 16:09:49.696 [Task submitter] DEBUG nextflow.executor.GridTaskHandler - [SLURM] submitted process auto:run_benchmark:run_wf:runEachWf:harmonypy:processWf:harmonypy_process (human_blood_mass_cytometry.harmonypy) > jobId: 25169450; workDir: /vast/scratch/users/putri.g/nextflow/work/6f/488d302c4e33171f2d69c8107b7a9c


# to pull log lines with:
# Sep-01 00:25:30.723 [AWSBatch-executor-26] DEBUG n.c.aws.batch.AwsBatchTaskHandler - [AWS BATCH] Process `auto:run_benchmark:run_wf:runEachWf:harmonypy:processWf:harmonypy_process (human_blood_mass_cytometry.harmonypy)` submitted > job=70af3d53-8070-45a6-96f4-fd37580b9d00; work-dir=s3://openproblems-work/scratch/2e8BJJovTcoBGj/60/c803767cf237028c908fadf447b1b9
# 1. Match only lines that contain "[AWS BATCH] Process `...`"
# 2. Extract the last part inside the backticks (e.g. "extract_uns_metadata_process (mouse_spleen_flow_cytometry)")
# 3. Extract the work-dir S3 path
# zooming into the ([\w\d_]+ \([^)]+\)) part:
# | Component  | Meaning                                                                                       |
# | ---------- | --------------------------------------------------------------------------------------------- |
# | `[\w\d_]+` | Matches the process name: letters, digits, underscores (e.g., `extract_uns_metadata_process`) |
# | ` `        | Space between process and dataset                                                             |
# | `\(`       | Opening parenthesis (escaped)                                                                 |
# | `[^)]+`    | Match everything inside the parentheses (dataset name), up to the next `)`                    |
# | `\)`       | Closing parenthesis                                                                           |
# that captures harmonypy_process (human_blood_mass_cytometry.harmonypy)
# the work-dir part:
# | Component       | Meaning                                                                       |
# | --------------- | ----------------------------------------------------------------------------- |
# | `.*?`           | Non-greedy match of any characters between the backtick block and `work-dir=` |
# | `work-dir=`     | Literal match                                                                 |
# | `(s3://[^ ;]+)` | The full S3 path. Match any characters except space or semicolon              |


aws_batch_pattern = re.compile(
    r"\[AWS BATCH\] Process `.*:([\w\d_]+ \([^)]+\))`.*?work-dir=(s3://[^ ;]+)"
)

slurm_pattern = re.compile(
    r"\[SLURM\] submitted process .*:([\w\d_]+ \([^)]+\)).*?workDir: ([^ ]+)"
)

is_aws = False

task_s3_bucket_map = []
with open(args.log_file_path, "r") as log_file:
    for line_number, line in enumerate(log_file, start=1):
        if is_aws:
            match = aws_batch_pattern.search(line)
        else:
            match = slurm_pattern.search(line)
        if match:
            # e.g., "extract_uns_metadata_process (mouse_spleen_flow_cytometry)"
            task_component = match.group(1)
            # task can be like this. so can pull the dataset name, task, and whether it is method or metric
            # harmonypy_process (human_blood_mass_cytometry.harmonypy)
            # flowsom_mapping_similarity_process (mouse_spleen_flow_cytometry.cytonorm_no_controls_to_mid.flowsom_mapping_similarity)
            process = task_component.split(" ")[0]
            task_info = task_component.split(" ")[1].strip("()").split(".")
            dataset_name = task_info[0]
            metric_name = None
            method_name = None
            if 2 <= len(task_info) <= 3:
                method_name = task_info[1]
            if len(task_info) == 3:
                metric_name = task_info[2]

            # e.g., "s3://openproblems-work/scratch/2e8"
            s3_path = match.group(2).strip()
            task_s3_bucket_map.append(
                (line_number, process, dataset_name, method_name, metric_name, s3_path)
            )

task_s3_bucket_map = pd.DataFrame(
    task_s3_bucket_map,
    columns=["line", "process_name", "dataset", "method", "metric", "s3_path"],
)
task_s3_bucket_map.sort_values(by="process_name", inplace=True)
task_s3_bucket_map.to_csv(args.output_csv_path, index=False)
