from collections import defaultdict

import anndata as ad
import cytonormpy as cnp
import numpy as np
import pandas as pd

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
    "input_integrated": "resources_test/task_cyto_batch_integration/starter_file/integrated.h5ad",
    "input_unintegrated": "resources_test/task_cyto_batch_integration/starter_file/unintegrated.h5ad",
    "input_validation": "resources_test/task_cyto_batch_integration/starter_file/validation.h5ad",
    "output": "output.h5ad",
}
meta = {"name": "emd_per_samples"}
## VIASH END

print("Reading input files", flush=True)

input_integrated = ad.read_h5ad(par["input_integrated"])
input_unintegrated = ad.read_h5ad(par["input_unintegrated"])
input_validation = ad.read_h5ad(par["input_validation"])

# add all obs columns in the unintegrated data to integrated data
# loc should order the obs based on obs_names
input_integrated.obs = input_unintegrated.obs.loc[input_integrated.obs_names]

# TODO uncomment me if you want to have some samples in validation but not in integrated
# input_validation = input_integrated[
#     input_integrated.obs["sample"] == "Tube1_Batch2_WT"
# ].copy()
# input_validation.layers["preprocessed"] = input_validation.layers["integrated"]
# del input_validation.layers["integrated"]
# input_integrated = input_integrated[
#     input_integrated.obs["sample"].isin(
#         ["Tube1_Batch1_WT", "Tube2_Batch1_KO", "Tube2_Batch2_KO"]
#     )
# ].copy()

markers_to_assess = input_unintegrated.var[
    input_unintegrated.var["to_correct"]
].index.to_numpy()

print("Extracting samples to compute the metric for", flush=True)

# Perhaps overengineering this but the idea is to get all the sample names from
# the integrated data, get the donors name from the combination of input_validation
# and input_unintegrated. Why? because depending on how we structure the input
# for an algorithm we may include matching samples from same donor (which will then
# appear in input_unintegrated) or not (which will then only appear in input_validation).
# then get all samples belong to that donor, double check the batch id,
# and compute the emd across the matching samples

# get all donors-sample map from data
# for some dataset "scheme", matching samples from the same donor may not be given
# as an input to the batch correction algorithm.
# thus they will only exist in validation. we need this for each donor processed
# by the batch correction algorithm.
# Assumption here is that for a given donor, there must be at least 1 sample
# processed by the batch correction algorithm, i.e., in uncorrected_censored (or uncorrected).
sample_donor_map_arr = np.concatenate(
    (
        np.unique(
            input_validation.obs[["sample", "donor"]].to_numpy().astype(str), axis=0
        ),
        np.unique(
            input_integrated.obs[["sample", "donor"]].to_numpy().astype(str), axis=0
        ),
    )
)
# keep just the donors which samples were processed by the batch correction algorithm.
# Why? we only need to compute the metric score for matching samples for donors that
# are processed by the batch correction algorithm.
# TODO can use list instead of set, because validation and unintegrated should not
# share common samples? Keep set for now.
sample_donor_dict = defaultdict(set)

for sample_donor_pair in sample_donor_map_arr:
    sample = sample_donor_pair[0]
    donor = sample_donor_pair[1]
    # print("sample: %s, donor: %s" % (sample, donor))
    sample_donor_dict[donor].add(sample_donor_pair[0])

emd_integrated_per_donor = []
samples_for_emd = []

for donor, samples in sample_donor_dict.items():

    # test 1 donor for now.
    # donor = list(sample_donor_dict.keys())[0]
    # samples = sample_donor_dict[donor]

    # the data for some samples may only be in the validation because they were not
    # included as input for the batch correction algorithm.
    input_integrated_donor = input_integrated[
        input_integrated.obs["sample"].isin(samples)
    ]
    input_validation_donor = input_validation[
        input_validation.obs["sample"].isin(samples)
    ]

    # keep just the markers we need to assess
    input_integrated_donor = input_integrated_donor[:, markers_to_assess]
    input_validation_donor = input_validation_donor[:, markers_to_assess]

    input_for_emd = ad.concat([input_integrated_donor, input_validation_donor])
    input_for_emd.layers["integrated"] = np.concatenate(
        (
            input_integrated_donor.layers["integrated"],
            input_validation_donor.layers["preprocessed"],
        )
    )

    # sanity check.. we need to have at least 2 batches..
    if len(np.unique(input_for_emd.obs["batch"])) < 2:
        raise Exception(
            "The samples to calculate EMD_per_sample for donor %s only belong to 1 batch?"
            % donor
        )

    print("Computing EMD for donor %s" % donor, flush=True)

    # have to change the "sample" column to file_name for emd_comparison_from_anndata to work.
    # Otherwise the _calculate_emd_per_frame used in cytonormpy will error because they
    # harcoded the column file_name and use it in assert.
    # See line 176 of https://github.com/TarikExner/CytoNormPy/blob/main/cytonormpy/_evaluation/_emd_utils.py#L173
    input_for_emd.obs["file_name"] = input_for_emd.obs["sample"]

    emd_integrated = cnp.emd_from_anndata(
        adata=input_for_emd,
        file_list=list(samples),
        channels=markers_to_assess,
        layer="integrated",
        sample_identifier_column="file_name",
    ).reset_index()

    del emd_integrated["label"]

    emd_integrated_per_donor.append(emd_integrated)
    samples_for_emd.append(";".join(samples))

emd_integrated_per_donor = pd.concat(emd_integrated_per_donor, ignore_index=True)

max_emd = np.max(emd_integrated_per_donor)
# mean_emd = np.mean(emd_integrated_per_donor)

print("Assembling output AnnData", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": input_integrated.uns["dataset_id"],
        "method_id": input_integrated.uns["method_id"],
        "metric_ids": ["max_emd"],
        "metric_values": [max_emd],
        "sample_ids": samples_for_emd,
    }
)

print("Write output AnnData to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
